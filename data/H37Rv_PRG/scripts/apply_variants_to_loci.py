from typing import TextIO, Dict, NamedTuple, Tuple
from pathlib import Path

from intervaltree import IntervalTree
from cyvcf2 import VCF, Writer
import click
import logging


class Record(NamedTuple):
    name: str
    comment: str
    seq: str


def load_loci_info(instream: TextIO) -> IntervalTree:
    intervals = []
    _ = next(instream)  # skip header
    for row in map(str.rstrip, instream):
        fields = row.split(",", maxsplit=6)
        start, end, name = fields[2:5]

        intervals.append((int(start) + 1, int(end), name))

    return IntervalTree.from_tuples(intervals)


def get_record_for_loci(loci_path: Path) -> Record:
    seq = ""
    name = ""
    comment = ""
    for line in loci_path.read_text().splitlines():
        if line.startswith(">"):
            name, comment = line[1:].split(" ", maxsplit=1)
            continue
        seq += line.strip()
    if not name:
        raise KeyError(f"Couldn't parse a header from {loci_path}")
    return Record(name, comment, seq)


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--vcf-path",
    help="VCF file to apply variants from.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "--loci-dir",
    type=click.Path(dir_okay=True, file_okay=False),
    help=(
        "The directory containing the loci template/reference sequences to apply "
        "variants to."
    ),
    required=True,
)
@click.option(
    "--loci-info",
    help=(
        "File containing loci info. Assumed format is VCF with columns: filename, "
        "type, start, end, name, contig."
    ),
    type=click.File(),
    required=True,
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(dir_okay=True, file_okay=False, writable=True),
    help="The directory to write the applied variant fasta files to.",
    default=".",
    show_default=True,
)
@click.option(
    "--chrom", help="Chromosome name in VCF.", default="NC_000962.3", show_default=True
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    vcf_path: str,
    loci_info: TextIO,
    outdir: str,
    verbose: bool,
    chrom: str,
    loci_dir: str,
):
    """Associate information about loci to the relevant VCF records based on position.

    This script will add three new INFO fields to the VCF. The INFO fields are
    loci_name, start, and end. See the VCF header entries for these fields for more
    information."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )
    loci_dir = Path(loci_dir)
    outdir = Path(outdir)

    if not outdir.exists():
        outdir.mkdir()

    logging.info("Loading loci info file...")
    ivtree = load_loci_info(loci_info)
    logging.info(f"Loaded info for {len(ivtree)} loci")

    vcf = VCF(vcf_path)

    logging.info("Applying variants to loci...")
    for iv in ivtree:
        start = iv.begin
        end = iv.end
        loci_name = iv.data
        loci_path = loci_dir / chrom / f"{loci_name}.fa"
        record = get_record_for_loci(loci_path)

        outpath = outdir / loci_path.name
        with outpath.open("w") as outstream:
            # todo: write ref seq to new file
            count = 0
            for variant in vcf(f"{chrom}:{start}-{end}"):
                count += 1
                loci_var_start = variant.POS - start  # 0-based
                loci_var_end = # todo
                # todo: check ref matches loci seq
                # todo: apply variant to seq - https://github.com/mbhall88/pandora_simulations/blob/68e0149099d3a3ba0b56d902017f2216b08b913b/scripts/apply_vcf.py#L45
                # todo: create unique id for variant to use as header id
                # todo: write new seq to file



        if count < 1:
            logging.info(f"No records associated with loci {loci_name}")
        else:
            logging.debug(f"{count} record(s) associated with loci {loci_name}")

    vcf.close()
    logging.info("All done.")


if __name__ == "__main__":
    main()
