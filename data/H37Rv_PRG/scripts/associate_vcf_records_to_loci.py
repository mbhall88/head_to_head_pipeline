from typing import TextIO

from intervaltree import IntervalTree
from cyvcf2 import VCF, Writer
import click
import logging

LOCI_ID = "LOCI"
START_ID = "START"
END_ID = "STOP"


def load_loci_info(instream: TextIO) -> IntervalTree:
    intervals = []
    _ = next(instream)  # skip header
    for row in map(str.rstrip, instream):
        fields = row.split(",", maxsplit=6)
        start, end, name = fields[2:5]

        intervals.append((int(start) + 1, int(end), name))

    return IntervalTree.from_tuples(intervals)


def write_new_info_fields(vcf_writer: Writer):
    vcf_writer.add_info_to_header(
        {
            "ID": LOCI_ID,
            "Description": "name of overlapping loci",
            "Type": "String",
            "Number": "1",
        }
    )
    vcf_writer.add_info_to_header(
        {
            "ID": START_ID,
            "Description": "Loci start position; 1-based inclusive",
            "Type": "Integer",
            "Number": "1",
        }
    )
    vcf_writer.add_info_to_header(
        {
            "ID": END_ID,
            "Description": "Loci end position; 1-based inclusive",
            "Type": "Integer",
            "Number": "1",
        }
    )


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--vcf-path",
    help="VCF file to associate loci to.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "--loci-info",
    help=(
        "File containing loci info. Assumed format is VCF with columns: filename, type,"
        " start, end, name, contig."
    ),
    type=click.File(),
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, writable=True, allow_dash=True),
    help="The filepath to write the updated VCF to.",
    default="-",
    show_default=True,
)
@click.option(
    "--chrom", help="Chromosome name in VCF.", default="NC_000962.3", show_default=True
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(vcf_path: str, loci_info: TextIO, output: str, verbose: bool, chrom: str):
    """Associate information about loci to the relevant VCF records based on position.

    This script will add three new INFO fields to the VCF. The INFO fields are
    loci_name, start, and end. See the VCF header entries for these fields for more
    information."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    logging.info("Loading loci info file...")
    ivtree = load_loci_info(loci_info)
    vcf = VCF(vcf_path)
    vcf_writer = Writer(output, vcf)
    write_new_info_fields(vcf_writer)

    logging.info("Associating loci info to VCF records...")
    for iv in ivtree:
        start = iv.begin
        end = iv.end
        name = iv.data
        count = 0
        for record in vcf(f"{chrom}:{start}-{end}"):
            count += 1
            record.INFO[LOCI_ID] = name
            record.INFO[START_ID] = start
            record.INFO[END_ID] = end
            vcf_writer.write_record(record)

        if count < 1:
            logging.info(f"No records associated with loci {name}")
        else:
            logging.debug(f"{count} record(s) associated with loci {name}")

    vcf_writer.close()
    vcf.close()
    logging.info("All done.")


if __name__ == "__main__":
    main()
