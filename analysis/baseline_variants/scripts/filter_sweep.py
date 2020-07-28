import logging
import subprocess
import json
from itertools import chain, product
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Tuple

import click
import pandas as pd


def interleave(a: Iterable, b: Iterable) -> Iterable:
    return chain(*zip(a, b))


class Filter:
    def __init__(self, script: str, outdir: Path, force: bool = False):
        self.script = script
        self.outdir = outdir
        self.force = force
        self.opts = ["--min-depth", "--max-depth", "--min-qual", "--min-strand-bias"]

    def run(self, infile: Path, args: Tuple[int, int, int, int]) -> Path:
        """Returns the output file path"""
        no_filter = all(x == 0 for x in args)
        if no_filter:
            return infile

        all_filters = all(args)
        if all_filters:
            outname = infile.with_suffix(".all_filters.bcf").name
        else:
            idx_of_filter_param = next((i for i, j in enumerate(args) if j), None)
            filter_name = (
                self.opts[idx_of_filter_param].strip("-")
                + f"={args[idx_of_filter_param]}"
            )
            outname = infile.with_suffix(f".{filter_name}.bcf").name
        outfile = self.outdir / outname

        if (not self.force) and outfile.exists():
            logging.info(f"{outfile} already exists. Will not run filter script again.")
            return outfile

        script_args = list(map(str, interleave(self.opts, args)))
        script_args.extend(["-i", str(infile), "-P", "-o", str(outfile)])
        cmd = " ".join(["python", self.script, *script_args])
        logging.info(f"Running command: {cmd}")
        proc = subprocess.run(cmd.split(), capture_output=True, check=True)
        assert proc.returncode == 0
        assert outfile.exists()
        logging.debug(proc.stderr.decode("utf-8"))
        return outfile


class Concordance:
    def __init__(
        self, script: str, outdir: Path, truth_vcf: str, mask: str, force: bool = False
    ):
        self.script = script
        self.outdir = outdir
        self.truth_vcf = truth_vcf
        self.mask = mask
        self.force = force

    def run(self, infile: Path) -> Path:
        no_filter = not any(s in str(infile) for s in ["filter", "="])
        if no_filter:
            outprefix = infile.with_suffix(".no_filter.concordance").name
        else:
            outprefix = infile.with_suffix(".concordance").name
        outcsv = self.outdir / f"{outprefix}.csv"
        outjson = self.outdir / f"{outprefix}.json"

        if (not self.force) and outjson.exists():
            logging.info(
                f"{outjson} already exists. Will not run concordance script again."
            )
            return outjson

        script_args = [
            "-a",
            self.truth_vcf,
            "-b",
            str(infile),
            "-m",
            self.mask,
            "-c",
            str(outcsv),
            "-j",
            str(outjson),
        ]
        cmd = " ".join(["python", self.script, *script_args])
        logging.info(f"Running command: {cmd}")
        proc = subprocess.run(cmd.split(), capture_output=True, check=True)
        assert proc.returncode == 0
        assert outjson.exists()
        assert outcsv.exists()
        logging.debug(proc.stderr.decode("utf-8"))
        return outjson


@click.command()
@click.help_option("--help", "-h")
@click.argument("vcfs", type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.option(
    "-o",
    "--outdir",
    help="Directory for all output files.",
    type=click.Path(file_okay=False, writable=True),
    default=".",
    show_default=True,
)
@click.option(
    "-M",
    "--markdown-file",
    help=(
        "Filename for output markdown file (within --outdir) containing the table of "
        "results."
    ),
    default="results.md",
    show_default=True,
)
@click.option(
    "-S",
    "--filter-script",
    help="Path to the script for applying filters.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-C",
    "--concordance-script",
    help="Path to the script for checking concordance.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-a",
    "--truth-vcf",
    help="VCF file that is considered 'truth'.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-m",
    "--mask",
    help="BED file containing positions to ignore.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-d",
    "--min-depth",
    help="Minimum depth as a percentage of the expected (median) depth.",
    default=0.2,
    show_default=True,
)
@click.option(
    "-D",
    "--max-depth",
    help="Maximum depth as a fraction of the expected (median) depth.",
    default=2.0,
    show_default=True,
)
@click.option(
    "-s",
    "--min-strand-bias",
    help="Minimum strand bias percentage allowed on the called allele.",
    type=click.IntRange(0, 50),
    metavar="INT",
    default=25,
    show_default=True,
)
@click.option(
    "-q",
    "--min-qual",
    help=f"Filter a variant if QUAL is less than INT.",
    default=0.0,
    show_default=True,
)
@click.option(
    "--force", "-f", help="Overwrite output files if they already exist.", is_flag=True
)
@click.option(
    "-j",
    "--jobs",
    help="Number of subprocesses to use. Set to 0 to use all available.",
    default=1,
    show_default=True,
)
@click.option(
    "--float-precision",
    help="How many positions of precision in results.",
    default=8,
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    vcfs: Iterable[str],
    outdir: str,
    filter_script: str,
    force: bool,
    verbose: bool,
    jobs: int,
    min_qual: float,
    min_depth: float,
    max_depth: float,
    min_strand_bias: int,
    mask: str,
    concordance_script: str,
    truth_vcf: str,
    markdown_file: str,
    float_precision: int,
):
    """Apply filters to VCFs and plot concordance differences.

    VCFS: A list of VCF filters to apply filters to and plot.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    if jobs == 0:
        logging.debug("All available CPUs will be used.")
        jobs = None

    outdir = Path(outdir)
    if not outdir.exists():
        logging.debug("Outdir does not exist. Creating...")
        outdir.mkdir()

    filter_args = [
        (0, 0, 0, 0),
        (min_depth, 0, 0, 0),
        (0, max_depth, 0, 0),
        (0, 0, min_qual, 0),
        (0, 0, 0, min_strand_bias),
        (min_depth, max_depth, min_qual, min_strand_bias),
    ]
    args = product(map(Path, vcfs), filter_args)
    filter_runner = Filter(filter_script, outdir, force)
    concord_runner = Concordance(concordance_script, outdir, truth_vcf, mask, force)

    with Pool(processes=jobs) as pool:
        filtered_vcfs = pool.starmap(filter_runner.run, args)
        concordance_jsons = pool.map(concord_runner.run, filtered_vcfs)

    data = dict()
    for p in concordance_jsons:
        sample, rest = str(p).split(".calls.")
        filter_type = rest.replace(".concordance.json", "")
        with p.open() as fp:
            data[filter_type] = json.load(fp)

    df = pd.DataFrame(data).T
    df["filter"] = df.index

    out_md = outdir / markdown_file
    with pd.option_context("display.float_format", f"{{:0.{float_precision}f}}".format):
        md = df.to_markdown(tablefmt="github", floatfmt=f".{float_precision - 2}%")
    out_md.write_text(md)

    logging.info(f"Results written to {out_md}")
    logging.info("Done!")


if __name__ == "__main__":
    main()
