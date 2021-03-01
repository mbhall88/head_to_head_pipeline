from pathlib import Path
from typing import List, Collection, Callable
import pandas as pd
import networkx as nx
import numpy as np
import click
import logging

from scipy import stats

WIDTH = 960
HEIGHT = 720
DELIM = ","
PAIR_IDX = ("sample1", "sample2")


class AsymmetrixMatrixError(Exception):
    pass


def parse_threshold(ctx, param, value):
    try:
        return [int(t) for t in value.split(",")]
    except ValueError:
        raise click.BadParameter(
            "threshold values must be integers, and if providing multiple, must be comma-separated"
        )


def parse_samples(ctx, param, value):
    return value.split(",")


def load_matrix(fpath, delim: str = DELIM, name: str = "") -> pd.Series:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        names = np.array(header.split(delim)[1:])
        idx = np.argsort(names)
        sorted_names = names[idx]
        for row in map(str.rstrip, instream):
            # sort row according to the name sorting
            sorted_row = np.array(row.split(delim)[1:], dtype=np.int)[idx]
            matrix.append(sorted_row)

    sorted_matrix = np.array(matrix)[idx]
    n_samples = len(sorted_names)
    diagonal_is_zero = all(sorted_matrix[i, i] == 0 for i in range(n_samples))
    if not diagonal_is_zero:
        raise AsymmetrixMatrixError("Distance matrix diagonal is not all zero")

    matrix_is_symmetric = np.allclose(sorted_matrix, sorted_matrix.T)
    if not matrix_is_symmetric:
        raise AsymmetrixMatrixError("Distance matrix is not symmetric")

    mx = pd.DataFrame(sorted_matrix, columns=sorted_names, index=sorted_names)
    # remove the lower triangle of the matrix and the middle diagonal
    mx = mx.where(np.triu(np.ones(mx.shape), k=1).astype(np.bool))
    mx = mx.stack().rename(name).astype(int)
    mx = mx.rename_axis(PAIR_IDX)

    return mx


def dist_matrix_to_graph(mx: pd.Series, threshold: int) -> nx.Graph:
    edges = [(s1, s2, dist) for (s1, s2), dist in mx.iteritems() if dist <= threshold]
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges)
    return graph


def fit(xs: List[int], ys: List[int]) -> Callable:
    """Determines the line of best fit for the given distances and returns the linear equation
    as a function that takes a threshold, x, and returns the equivalent threshold for the data
    passed to this function.
    Note: xs should be the 'truth' distance values"""
    slope, intercept, rvalue, pvalue, stderr = stats.linregress(xs, ys)

    def threshold_converter(threshold: float) -> float:
        """y = mx + c, where x is the threshold"""
        return slope * threshold + intercept

    return threshold_converter


def cluster_info(graph: nx.Graph, name: str, threshold: int) -> None:
    pad = 40
    logging.info("=" * pad)
    logging.info(f"Cluster information for {name}")
    logging.info("-" * pad)
    logging.info(f"Cluster SNP distance threshold: {threshold}")
    logging.info(f"Number of clusters: {nx.number_connected_components(graph)}\n")
    logging.info("Individual cluster information")
    logging.info("-" * pad)
    for cluster in sorted(nx.connected_components(graph), key=len, reverse=True):
        logging.info(f"Cluster size: {len(cluster)}")
        logging.info(f"Members:")
        logging.info(sorted(cluster))


@click.command()
@click.help_option("--help", "-h")
@click.argument(
    "target",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    nargs=1,
    required=True,
)
@click.argument(
    "queries",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    required=True,
)
@click.option(
    "-S",
    "--samples",
    help=(
        "Names for the matrices. If not given, the stem of the filename(s) will be "
        "used. For multiple names, use a comma-separated list in the same order as "
        "the file paths"
    ),
    callback=parse_samples,
)
@click.option(
    "-o",
    "--outdir",
    help="Path to save output files to.",
    type=click.Path(file_okay=False, writable=True, resolve_path=True),
    default=".",
    show_default=True,
)
@click.option(
    "-d",
    "--delim",
    help="Delimiter used in the matrix.",
    default=DELIM,
    show_default=True,
)
@click.option("--width", help="Plot width in pixels", default=WIDTH, show_default=True)
@click.option(
    "--height", help="Plot height in pixels", default=HEIGHT, show_default=True
)
@click.option(
    "-T",
    "--threshold",
    help=(
        "SNP distance threshold used to define clusters. Samples with distance <= "
        "threshold are considered in the same cluster. This should be a single value "
        "or a comma-separated list of thresholds for each matrix"
    ),
    type=str,
    required=True,
    callback=parse_threshold,
    metavar="INT[,INT...]",
)
@click.option(
    "-A",
    "--adaptive-threshold",
    help=(
        "Calculate the thresholds for the query samples based on the linear "
        "relationship between the target and query pairwise distances. Only those pairs "
        "with a target distance less than or equal to the value passed will be used. "
        "Set to 0 to use all pairs."
    ),
    type=int,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    target: str,
    queries: Collection[str],
    verbose: bool,
    samples: List[str],
    outdir: str,
    delim: str,
    width: int,
    height: int,
    threshold: List[int],
    adaptive_threshold: int,
):
    """Produce visualisations and metrics for the given distance matrices

    TARGET: The distance matrix to compare the others to\n
    QUERIES: Distance matri{x,ces} to compare against the target
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )
    n_matrices = len(queries) + 1

    if not samples:
        samples = [p.name.split(".")[0] for p in map(Path, [target, *queries])]

    if (n_samples := len(samples)) != n_matrices:
        raise ValueError(f"Received {n_matrices}, but only {n_samples}")

    if (n_thresholds := len(threshold)) > 1 and n_thresholds != n_matrices:
        raise ValueError(
            f"Expected either 1 or {n_matrices} threshold values but got {n_thresholds}"
        )

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    target_name = samples[0]
    logging.info("Loading target distance matrix...")
    target_mtx = load_matrix(target, delim=delim, name=target_name)

    logging.info("Loading query matrices...")
    query_mtxs = []
    for name, path in zip(samples[1:], queries):
        query_mtxs.append(load_matrix(path, delim=delim, name=name))

    matrices_have_same_sample_pairs = all(
        target_mtx.index.equals(m.index) for m in query_mtxs
    )
    if not matrices_have_same_sample_pairs:
        raise IndexError("Distance matrices do not have identical sample pairs")

    keep_idx = target_mtx.index
    if adaptive_threshold is not None:
        if adaptive_threshold > 0:
            logging.info(
                f"Restricting to samples pairs with a {target_name} distance of <= "
                f"{adaptive_threshold}..."
            )
            keep_idx = target_mtx <= adaptive_threshold

        expected_dists = target_mtx[keep_idx].to_list()

        for i, mtx in enumerate(query_mtxs, start=1):
            name = samples[i]
            query_dists = mtx[keep_idx].to_list()
            query_threshold_converter = fit(expected_dists, query_dists)
            query_threshold = int(round(query_threshold_converter(threshold[0])))
            if n_thresholds < n_matrices:
                threshold.insert(i, query_threshold)
            else:
                threshold[i] = query_threshold
            logging.info(f"Adaptive threshold used for {name} is {query_threshold}")
    elif n_thresholds == 1:
        threshold = [threshold[0]] * n_matrices

    logging.info("Reducing graphs to on those nodes with edges <= threshold")

    target_graph = dist_matrix_to_graph(target_mtx, threshold=threshold[0])
    cluster_info(target_graph, name=target_name, threshold=threshold[0])
    query_graphs = []
    for mx, t, name in zip(query_mtxs, threshold[1:], samples[1:]):
        query_graph = dist_matrix_to_graph(mx, threshold=t)
        cluster_info(query_graph, name=name, threshold=t)
        query_graphs.append(query_graph)


if __name__ == "__main__":
    main()
