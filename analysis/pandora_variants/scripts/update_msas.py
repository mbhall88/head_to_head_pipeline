import logging
from typing import TextIO

import click


@click.command()
@click.help_option("--help", "-h")
@click.argument()#todo
@click.option(
    "-i",
    "--vcf",
    help="VCF file generate the consensus for.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-f",
    "--ref",
    help="Reference fasta file.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=False),
    required=True,
)
@click.option(
    "-m",
    "--mask",
    "bedfile",
    help="BED file containing positions to mask.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-o",
    "--output",
    help="File path of consensus fasta to output.",
    type=click.File(mode="w"),
    default="-",
    show_default=True,
)
@click.option(
    "-H",
    "--het-default",
    help=(
        "Which heterozygous base should be used? If 'none', an 'N' will be used in the "
        "consensus sequence"
    ),
    type=click.Choice(DEFAULT_HET_CHOICES, case_sensitive=False),
    default=DEFAULT_HET_CHOICES[0],
    show_default=True,
)
@click.option(
    "-I",
    "--ignore",
    help=(
        "Which types of variants should be ignored (i.e. replaced with 'N')? If not "
        "ignoring filters/mask and the position is called ALT, then ALT will be used. "
        "For null/missing the REF base will be used. This option can be specified "
        "multiple times. For example, to ignore null and mask positions use: -I null "
        "-I mask"
    ),
    type=click.Choice(
        ("none", "all", "missing", "null", "filter", "mask"), case_sensitive=False
    ),
    default=["all"],
    show_default=True,
    multiple=True,
)
@click.option(
    "-s",
    "--sample-id",
    help=(
        "Sample identifier to use for the output fasta file. If not specified, the "
        "sample column name in the VCF will be used. If more than one contig is "
        "present, it will be combined with the sample ID (i.e. CONTIG|SAMPLE_ID)"
    ),
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    vcf: str,
    ref: str,
    output: TextIO,
    verbose: bool,
    bedfile: str,
    het_default: str,
    ignore: tuple,
    sample_id: str,
):
    """This script will produce a consensus fasta file from a VCF file.
    The VCF must contain ONLY SNVs; at this stage it cannot handle any other variant
    type. It will also only use the first sample in the VCF currently.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )


if __name__ == "__main__":
    main()
