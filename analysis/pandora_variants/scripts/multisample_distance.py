import logging
from itertools import combinations
from typing import Dict, Tuple, Callable

import click
import numpy as np
import pandas as pd


class AsymmetrixMatrixError(Exception):
    pass


def distance_between(xs: pd.Series, ys: pd.Series, dist_func: Callable) -> int:
    return sum(dist_func(x, y) for x, y in zip(xs, ys))


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

    def dist(x: int, y: int) -> int:
        if x == y or any(gt in (null, filtered) for gt in [x, y]):
            return 0
        return 1

    logging.info("Loading genotype matrix...")

    gt_matrix = pd.read_csv(
        matrix,
        delimiter=in_delim,
        header=0,
        usecols=lambda colname: colname not in ("CHROM", "POS"),
    )

    samples = gt_matrix.columns.values.tolist()
    n_samples = len(samples)
    logging.info(f"Loaded genotype matrix for {n_samples} samples")
    pairs = combinations(samples, 2)
    pairwise_distances: Dict[Tuple[str, str], int] = dict()

    logging.info("Calculating distances between all pairs...")
    for a, b in pairs:
        logging.debug(f"Calculating distance for pair ({a}, {b})")
        pairwise_distances[(a, b)] = distance_between(gt_matrix[a], gt_matrix[b], dist)
    logging.info(f"Calculated distances for {len(pairwise_distances)} pairs")

    logging.info("Filling distance matrix...")
    dist_matrix = np.zeros(shape=(n_samples, n_samples), dtype=int)
    # loop through pair indices and fill matrix with distances
    for i, j in combinations(range(n_samples), 2):
        try:
            dist = pairwise_distances[(samples[i], samples[j])]
        except KeyError:
            dist = pairwise_distances[(samples[j], samples[i])]

        dist_matrix[i, j] = dist
        dist_matrix[j, i] = dist

    diagonal_is_zero = all(dist_matrix[i, i] == 0 for i in range(n_samples))
    if not diagonal_is_zero:
        AsymmetrixMatrixError("Diagonal of distance matrix is not all zeros")

    matrix_is_symmetric = np.allclose(dist_matrix, dist_matrix.T)
    if not matrix_is_symmetric:
        raise AsymmetrixMatrixError("Distance matrix is not symmetric")
    logging.info("Distance matrix is symmetrical")

    # write matrix to file
    pd.DataFrame(dist_matrix, index=samples, columns=samples).to_csv(
        output, index=True, header=True, sep=delim
    )
    logging.info(f"Distance matrix written to {output}")


if __name__ == "__main__":
    main()
