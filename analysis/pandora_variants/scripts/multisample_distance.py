import click
import pandas as pd
import logging


@click.command()
@click.help_option("--help", "-h")
@click.argument("matrix", type=click.Path(exists=True, dir_okay=False), nargs=1)
@click.option(
    "-o",
    "--output",
    help="Path to write distance matrix to",
    type=click.Path(exists=False, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-d",
    "--delim",
    help="Delimiter to use in the output matrix",
    default=",",
    show_default=True,
)
@click.option(
    "-F",
    "--in-delim",
    help="Delimiter used in the input matrix",
    default=",",
    show_default=True,
)
@click.option("-H", "--no-header", help="The input has no header line", is_flag=True)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    verbose: bool, matrix: str, output: str, delim: str, in_delim: str, no_header: bool
):
    """Produces a distance matrix based on a given genotype matrix.

    Note: The expected input columns are: CHROM,POS,Genotypes
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )


if __name__ == "__main__":
    main()
