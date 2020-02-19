#!/usr/bin/env python3

import argparse
import errno
import logging
import os
from pathlib import Path
import subprocess
from typing import Dict

DEFAULT_MAX_ITERATIONS = 10


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

            print(line, end="", file=f_out)


def get_iteration_files(
    iteration_num: int, options: argparse.Namespace
) -> Dict[str, Path]:
    files = {
        "done_file": Path(f"iteration.{iteration_num}.done"),
        "pilon_dir": Path(f"iteration.{iteration_num}.pilon"),
        "mapping_prefix": Path(f"iteration.{iteration_num}.map"),
    }

    if iteration_num == 1:
        files["ref_fasta"] = options.assembly_fasta
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
    reads1: Path, reads2: Path, ref_fasta: Path, outprefix: Path, threads=1
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
        f"samtools sort --threads {threads} --reference {ref_fasta} -o {sorted_bam} {sam_file}"
    )
    sam_file.unlink()
    log_and_run_command(f"samtools index {sorted_bam}")
    done_file.touch()
    return sorted_bam


def run_pilon(
    bam: Path,
    ref_fasta: Path,
    outdir: Path,
    java_jar: Path,
    java_xmx_opt: str,
    threads=1,
):
    fasta = outdir / "pilon.fasta"
    done_file = outdir / ".done"
    if done_file.exists():
        logging.info(f"Found done file {done_file}")
        assert os.path.exists(f"{fasta}.fai")

    assert java_jar.exists()
    log_and_run_command(
        f"java -Xmx{java_xmx_opt} -jar {java_jar} --outdir {outdir} --genome {ref_fasta} --frags {bam} --changes --threads {threads} --minmq 10 --minqual 10"
    )
    log_and_run_command(f"bwa index {fasta}")
    log_and_run_command(f"samtools faidx {fasta}")
    done_file.touch()


def number_of_pilon_changes(changes_file: Path) -> int:
    with changes_file.open() as f:
        return sum(1 for _ in f)


def check_file_exists(path: Path, filename_for_log: str):
    if path.exists():
        logging.info(f"Found {filename_for_log}: {path}")
    else:
        logging.info(f"Not found: {filename_for_log}: {path}")
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), path)


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Iteratively run pilon to correct assembly. Stops after "
            f"{DEFAULT_MAX_ITERATIONS} (default) iterations, or if no changes made."
        )
    )

    parser.add_argument(
        "--pilon_java_xmx",
        default="8G",
        help="Option used with -Xmx when running javar pilon.jar [%(default)s]",
        metavar="STRING",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use with BWA mem and pilon [%(default)s]",
        metavar="INT",
    )
    parser.add_argument(
        "--max_iterations",
        type=int,
        default=DEFAULT_MAX_ITERATIONS,
        help="Max number of iterations. [%(default)s]",
        metavar="INT",
    )
    parser.add_argument(
        "--pilon_jar", help="Name of pilon JAR file", type=Path, required=True
    )
    parser.add_argument(
        "--assembly_fasta",
        help="Name of input assembly fasta file",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--reads1", help="Name of input reads file 1", type=Path, required=True
    )
    parser.add_argument(
        "--reads2", help="Name of input reads file 2", type=Path, required=True
    )
    parser.add_argument(
        "--outdir",
        help=f"Name of output/working directory. Default [%(default)s]",
        type=Path,
        default=Path(),
    )
    parser.add_argument(
        "--final_fasta",
        help="Name of the final, polished fasta file. [%(default)s]",
        default="final.fasta",
    )
    options = parser.parse_args()

    options.assembly_fasta = options.assembly_fasta.resolve()
    options.reads1 = options.reads1.resolve()
    options.reads2 = options.reads2.resolve()
    options.outdir = options.outdir.resolve()
    options.pilon_jar = options.pilon_jar.resolve()

    options.outdir.mkdir(exist_ok=True)

    return options


def main():
    options = cli()
    os.chdir(options.outdir)

    log = logging.getLogger()
    log.setLevel(logging.INFO)
    logfile_handle = logging.FileHandler("log.txt", mode="w")
    log = logging.getLogger()
    formatter = logging.Formatter(
        "[pilon_iterative %(asctime)s %(levelname)s] %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S",
    )
    logfile_handle.setFormatter(formatter)
    log.addHandler(logfile_handle)

    check_file_exists(options.assembly_fasta, "assembly_fasta")
    check_file_exists(options.reads1, "reads1")
    check_file_exists(options.reads2, "reads2")

    for iteration_num in range(1, options.max_iterations + 1, 1):
        logging.info(f"Start iteration {iteration_num}")
        files = get_iteration_files(iteration_num, options)

        if files["done_file"].exists():
            logging.info(f'Found iteration done file {files["done_file"]}')
        else:
            bam = make_pilon_bam(
                options.reads1,
                options.reads2,
                files["ref_fasta"],
                files["mapping_prefix"],
                threads=options.threads,
            )
            run_pilon(
                bam,
                files["ref_fasta"],
                files["pilon_dir"],
                options.pilon_jar,
                options.pilon_java_xmx,
                threads=options.threads,
            )
            bam.unlink()
            Path(f"{bam}.bai").unlink()

        number_of_changes = number_of_pilon_changes(files["changes_file"])
        logging.info(
            f"Number of changes at iteration {iteration_num}: {number_of_changes}"
        )
        files["done_file"].touch()
        logging.info(f"End iteration {iteration_num}")

        if number_of_changes == 0:
            logging.info(f"No changes made in iteration {iteration_num}. Stopping")
            logging.info("Making final corrected file final.fasta with renamed contigs")
            remove_pilon_from_fasta_headers(
                files["corrected_fasta"], Path(options.final_fasta), iteration_num
            )
            break

    logging.info("Finished")


if __name__ == "__main__":
    main()
