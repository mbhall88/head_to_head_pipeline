import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Union

import click

PathLike = Union[Path, str, os.PathLike]
DEFAULTS = {
    "max-iterations": 10,
    "min-mapq": 10,
    "min-qual": 10,
    "fix": "all",
    "pilon-memory": "8G",
    "final-fasta": "final.fasta",
    "threads": 1,
}
XMX_REGEX = re.compile(r"^[0-9]+[gkm]$", re.IGNORECASE)
XMX_URL = "https://bit.ly/3c3KEot"


def validate_xmx(ctx: click.Context, param, value: str) -> str:
    if not XMX_REGEX.match(value):
        raise click.BadParameter(
            f"The format of {value} is incorrect. "
            "Accepted format is <int>[g|G|m|M|k|K]. e.g. 8g, 1000K, or 980m"
        )
    return value


def remove_pilon_from_fasta_headers(infile: Path, outfile: Path, number_of_pilons: int):
    """Each time you run pilon, it adds "_pilon" to the end of each
    contig name. This makes a new file with them removed, and adds
    ".{number_of_pilons}_pilon_iterations" instead"""
    expected_suffix = "_pilon" * number_of_pilons

    with infile.open() as f_in, outfile.open("w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                line = line.rstrip()
                assert line.endswith(expected_suffix)
                line = line[0 : -len(expected_suffix)]

            print(line.rstrip(), file=f_out)


def get_iteration_files(iteration_num: int, assembly_fasta: str) -> Dict[str, Path]:
    files = {
        "done_file": Path(f"iteration.{iteration_num}.done"),
        "pilon_dir": Path(f"iteration.{iteration_num}.pilon"),
        "mapping_prefix": Path(f"iteration.{iteration_num}.map"),
    }

    if iteration_num == 1:
        files["ref_fasta"] = Path(assembly_fasta)
    else:
        files["ref_fasta"] = (
            Path(f"iteration.{iteration_num - 1}.pilon") / "pilon.fasta"
        )

    files["corrected_fasta"] = files["pilon_dir"] / "pilon.fasta"
    files["changes_file"] = files["pilon_dir"] / "pilon.changes"
    return files


def log_and_run_command(command):
    logging.info(f"Run: {command}")
    completed_process = subprocess.run(
        command,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    if completed_process.returncode == 0:
        logging.info(f"Run ok: {command}")
    else:
        logging.error(f"Error running: {command}")
        logging.error(f"Return code: {completed_process.returncode}")
        logging.error(f"Output from stdout:\n{completed_process.stdout}")
        logging.error(f"Output from stderr:\n{completed_process.stderr}")
        raise Exception("Error in system call. Cannot continue")


def make_pilon_bam(
    reads1: PathLike,
    reads2: PathLike,
    ref_fasta: PathLike,
    outprefix: PathLike,
    threads=1,
) -> Path:
    sam_file = Path(f"{outprefix}.sam")
    sorted_bam = Path(f"{outprefix}.sorted.bam")
    done_file = Path(f"{outprefix}.done")
    if done_file.exists():
        logging.info(f"Found done file {done_file}")
        assert sorted_bam.exists()
        return sorted_bam

    log_and_run_command(
        f"bwa mem -t {threads} -x intractg {ref_fasta} {reads1} {reads2} > {sam_file}"
    )
    log_and_run_command(
        f"samtools sort --threads {threads} --reference {ref_fasta} "
        f"-o {sorted_bam} {sam_file}"
    )
    sam_file.unlink()
    log_and_run_command(f"samtools index {sorted_bam}")
    done_file.touch()
    return sorted_bam


class Pilon:
    def __init__(
        self,
        jarfile: PathLike,
        memory: str,
        threads: int,
        min_mapq: int,
        min_qual: int,
        fixes: str,
    ):
        self.jarfile = jarfile
        self.memory = memory
        self.threads = threads
        self.min_mapq = min_mapq
        self.min_qual = min_qual
        self.fixes = fixes

    @staticmethod
    def cleanup_checkpoints(directory: Path):
        checkpoint_regex = re.compile(r"iteration\.[0-9]+\.(done|map|pilon)(\.done)?")
        for path in directory.iterdir():
            if checkpoint_regex.match(path.name):
                logging.debug(f"Found checkpoint file {path}. Removing file...")
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()

    def generate_params(self) -> str:
        params = [
            f"-Xmx{self.memory}",
            f"-jar {self.jarfile}",
            f"--threads {self.threads}",
            f"--minmq {self.min_mapq}",
            f"--minmq {self.min_mapq}",
            f"--minqual {self.min_qual}",
            "--changes",
        ]
        return " ".join(params)

    def run(
        self, bam: PathLike, assembly: PathLike, outdir: Path,
    ):
        fasta = outdir / "pilon.fasta"
        done_file = outdir / ".done"
        if done_file.exists():
            logging.info(f"Found done file {done_file}")
            assert os.path.exists(f"{fasta}.fai")

        params = self.generate_params()
        log_and_run_command(
            f"java {params} --outdir {outdir} --genome {assembly} --frags {bam}"
        )
        log_and_run_command(f"bwa index {fasta}")
        log_and_run_command(f"samtools faidx {fasta}")
        done_file.touch()

    @staticmethod
    def number_of_changes(changes_file: Path) -> int:
        with changes_file.open() as f:
            return sum(1 for _ in f)


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-1",
    "--reads1",
    help="Name of input reads file 1.",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
)
@click.option(
    "-2",
    "--reads2",
    help="Name of input reads file 2.",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
)
@click.option(
    "-o",
    "--outdir",
    help="Name of output/working directory.",
    required=True,
    type=click.Path(exists=False, file_okay=False, resolve_path=True),
)
@click.option(
    "-f",
    "--final-fasta",
    help="Name of the final, polished fasta file.",
    default=DEFAULTS["final-fasta"],
    show_default=True,
)
@click.option(
    "-a",
    "--assembly",
    help="Name of input assembly fasta file.",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    required=True,
)
@click.option(
    "-j",
    "--pilon-jar",
    help="Name of pilon JAR file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    required=True,
)
@click.option(
    "-i",
    "--max-iterations",
    type=click.IntRange(min=1),
    default=DEFAULTS["max-iterations"],
    help="Max number of iterations.",
    show_default=True,
)
@click.option(
    "-t",
    "--threads",
    type=click.IntRange(min=1),
    default=DEFAULTS["threads"],
    help="Number of threads to use with BWA mem and pilon.",
    show_default=True,
)
@click.option(
    "-m",
    "--pilon-memory",
    default=DEFAULTS["pilon-memory"],
    help=(
        "Maximum memory allocation pool for running Pilon. It is passed as `java "
        f"-Xmx<value>`. Format: <int>[g|G|m|M|k|K]. See {XMX_URL} for more information."
    ),
    show_default=True,
    callback=validate_xmx,
)
@click.option(
    "--fix",
    help=(
        "A comma-separated list of categories of issues to try to fix. Refer to the "
        "pilon documentation for valid values."
    ),
    default=DEFAULTS["fix"],
    show_default=True,
)
@click.option(
    "--min-mapq",
    help="Minimum alignment mapping quality for a read to count in pileups",
    default=DEFAULTS["min-mapq"],
    show_default=True,
)
@click.option(
    "--min-qual",
    help="Minimum base quality to consider for pileups",
    default=DEFAULTS["min-qual"],
    show_default=True,
)
@click.option(
    "--force",
    help="Ignore any previous checkpoints and polish from iteration 1.",
    is_flag=True,
)
@click.option(
    "-l",
    "--log-file",
    help="Name of the log file.",
    default="log.txt",
    show_default=True,
)
def main(
    reads1: str,
    reads2: str,
    outdir: str,
    final_fasta: str,
    assembly: str,
    pilon_jar: str,
    max_iterations: int,
    threads: int,
    pilon_memory: str,
    log_file: str,
    fix: str,
    min_qual: int,
    min_mapq: int,
    force: bool,
):
    """Iteratively run pilon to correct an assembly with Illumina reads. Stops after a
    specified number of iterations, or if no changes were made on the last iteration.
    """
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    os.chdir(outdir)

    log = logging.getLogger()
    log.setLevel(logging.INFO)
    logfile_handle = logging.FileHandler(log_file, mode="w")
    log = logging.getLogger()
    formatter = logging.Formatter(
        "[pilon_iterative %(asctime)s %(levelname)s] %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S",
    )
    logfile_handle.setFormatter(formatter)
    log.addHandler(logfile_handle)

    pilon = Pilon(pilon_jar, pilon_memory, threads, min_mapq, min_qual, fix)
    if force:
        pilon.cleanup_checkpoints(outdir)

    for iteration_num in range(1, max_iterations + 1, 1):
        logging.info(f"Start iteration {iteration_num}")
        files = get_iteration_files(iteration_num, assembly)

        if files["done_file"].exists() and not force:
            logging.info(f'Found iteration done file {files["done_file"]}')
        else:
            bam = make_pilon_bam(
                reads1, reads2, files["ref_fasta"], files["mapping_prefix"], threads
            )
            pilon.run(bam, files["ref_fasta"], files["pilon_dir"])
            bam.unlink()
            Path(f"{bam}.bai").unlink()

        number_of_changes = pilon.number_of_changes(files["changes_file"])
        logging.info(
            f"Number of changes at iteration {iteration_num}: {number_of_changes}"
        )
        files["done_file"].touch()

        if number_of_changes == 0:
            logging.info(f"No changes made in iteration {iteration_num}. Stopping")
            break
        elif iteration_num >= max_iterations:
            logging.info(f"Reached max. iteration number of {max_iterations}. Stopping")
            break
        else:
            logging.info(f"End of iteration {iteration_num}. Continuing")

    logging.info("Making final corrected file final.fasta with renamed contigs")
    remove_pilon_from_fasta_headers(
        files["corrected_fasta"], Path(final_fasta), iteration_num
    )

    logging.info("Finished")


if __name__ == "__main__":
    main()
