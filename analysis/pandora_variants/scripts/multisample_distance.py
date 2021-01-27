import logging
from itertools import combinations
from typing import Dict, Tuple

import click
import pandas as pd


def dist(x: int, y: int) -> int:
    if x == y or x < 0 or y < 0:
        return 0
    return 1


def distance_between(xs: pd.Series, ys: pd.Series) -> int:
    return sum(dist(x, y) for x, y in zip(xs, ys))


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
@click.option(
    "--null", help="Value used for null genotypes", default=-1, show_default=True
)
@click.option(
    "--filtered",
    help="Value used for filtered genotypes",
    default=-2,
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    verbose: bool,
    matrix: str,
    output: str,
    delim: str,
    in_delim: str,
    null: int,
    filtered: int,
):
    """Produces a distance matrix based on a given genotype matrix.

    Note: The expected input columns are: CHROM,POS,Genotypes
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )
    gt_matrix = pd.read_csv(
        matrix,
        delimiter=delim,
        header=0,
        usecols=lambda colname: colname not in ("CHROM", "POS"),
    )

    samples = gt_matrix.columns.values.tolist()
    pairs = combinations(samples, 2)
    pairwise_distances: Dict[Tuple[str, str], int] = dict()

    for x, y in pairs:
        pairwise_distances[(x, y)] = distance_between(gt_matrix[x], gt_matrix[y])

    # todo: turn dict into dist matrix


if __name__ == "__main__":
    main()
