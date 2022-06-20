import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import numpy as np
from itertools import chain
from typing import Tuple, Set

PAIR_IDX = ("sample1", "sample2")


class AsymmetrixMatrixError(Exception):
    pass


def load_matrix(fpath, delim: str = ",", name: str = "") -> pd.Series:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        names = np.array(header.split(delim)[1:])
        idx = np.argsort(names)
        sorted_names = names[idx]
        for row in map(str.rstrip, instream):
            # sort row according to the name sorting
            sorted_row = np.array(row.split(delim)[1:], dtype=int)[idx]
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
    mx = mx.where(np.triu(np.ones(mx.shape), k=1).astype(bool))
    mx = mx.stack().rename(name).astype(int)
    mx = mx.rename_axis(PAIR_IDX)

    return mx


def matrix_to_graph(
    mx: pd.Series, threshold: int, include_singletons: bool = False
) -> nx.Graph:
    edges = [(s1, s2, dist) for (s1, s2), dist in mx.iteritems() if dist <= threshold]
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges)
    if include_singletons:
        samples = set()
        for u in chain.from_iterable(mx.index):
            if u not in samples:
                graph.add_node(u)
                samples.add(u)
    return graph


def tversky_index(
    A: Set[str], B: Set[str], alpha: float = 1.0, beta: float = 1.0
) -> float:
    """If we set alpha and beta to 1 then we get the Jaccard Index.
    If we set alpha to 1 and beta to 0 we get something like recall.
    If we set alpha to 0 and beta to 1 we get something like precision.
    """
    size_of_intersection = len(A & B)
    A_weight = alpha * len(A - B)
    B_weight = beta * len(B - A)
    denominator = size_of_intersection + A_weight + B_weight

    try:
        return size_of_intersection / denominator
    except ZeroDivisionError:
        return 0


def set_precision(A: Set[str], B: Set[str]) -> float:
    return tversky_index(A, B, alpha=0, beta=1)


def set_recall(A: Set[str], B: Set[str]) -> float:
    return tversky_index(A, B, alpha=1, beta=0)


def excess_clustering_rate(A: Set[str], B: Set[str]) -> float:
    """What percentage of true singletons are clustered.
    What percentage of A is not in B
    """
    return len(A - B) / len(A)


def connected_components(G: nx.Graph, node: str) -> Set[str]:
    if node not in G:
        return set()
    return nx.node_connected_component(G, node)


def SACP_AND_SACR(G: nx.Graph, H: nx.Graph) -> Tuple[float, float]:
    G.remove_nodes_from(list(nx.isolates(G)))
    H.remove_nodes_from(list(nx.isolates(H)))

    expected_clusters = list(nx.connected_components(G))
    ppvs = []
    tprs = []
    for i, expected_cluster in enumerate(expected_clusters):
        cluster_ppv = []
        cluster_tpr = []
        for node in expected_cluster:
            actual_cluster = connected_components(H, node)
            tpr = set_recall(expected_cluster, actual_cluster)
            ppv = set_precision(expected_cluster, actual_cluster)
            cluster_ppv.append(ppv)
            cluster_tpr.append(tpr)
        ppvs.extend(cluster_ppv)
        tprs.extend(cluster_tpr)
    return np.mean(ppvs), np.mean(tprs)


# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi


def main():
    compass_mtx = load_matrix(snakemake.input.compass_matrix, name="illumina")
    samples = set(np.unique(np.array(list(chain.from_iterable(compass_mtx.index)))))
    M = compass_mtx
    N = load_matrix(snakemake.input.bcftools_matrix, name="bcftools")

    data = []
    distances = snakemake.params.illumina_thresholds
    ont_dists = list(range(*snakemake.params.nanopore_dist_range))
    for d in distances:
        T = matrix_to_graph(M, threshold=d)
        true_singletons = samples - set(T.nodes)
        for t in ont_dists:
            g = matrix_to_graph(N, threshold=t)
            test_singletons = samples - set(g.nodes)
            xcr = excess_clustering_rate(true_singletons, test_singletons)
            sacp, sacr = SACP_AND_SACR(
                matrix_to_graph(M, threshold=d, include_singletons=False),
                matrix_to_graph(N, threshold=t, include_singletons=False),
            )
            data.append((d, t, sacr, sacp, 1 - xcr))

    df = pd.DataFrame(
        data, columns=["distance", "threshold", "SACR", "SACP", "1-XCR"]
    ).melt(id_vars=["distance", "threshold"], var_name="metric")
    fig, axes = plt.subplots(nrows=len(distances), sharex=True)

    for ax, dist, t in zip(axes.flatten(), distances, ont_dists):
        subdata = df.query("distance == @dist")
        sns.lineplot(
            data=subdata,
            x="threshold",
            y="value",
            hue="metric",
            ax=ax,
            style="metric",
            markers=True,
        )
        ax.set_xticks(ont_dists)
        ax.set_xticklabels(ont_dists)
        ax.set_title(f"Distance = {dist}", fontdict=dict(size=10))

    fig.subplots_adjust(hspace=0.1)

    for fpath in snakemake.output:
        fig.savefig(fpath)


main()
