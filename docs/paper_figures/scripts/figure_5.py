import sys

sys.stderr = open(snakemake.log[0], "w")
from itertools import chain, combinations, product
from typing import Tuple, Set, List, Dict
from fractions import Fraction

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi

seed = snakemake.params.seed  # set to None for no seeding
rng = np.random.default_rng(seed=seed)
PAIR_IDX = ("sample1", "sample2")
DELIM = ","


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
        return round(size_of_intersection / denominator, 4)
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


def dist_matrix_to_graph(mx: pd.Series, threshold: int) -> nx.Graph:
    edges = [(s1, s2, dist) for (s1, s2), dist in mx.iteritems() if dist <= threshold]
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges)
    return graph


def truth_graph(mx: pd.Series, threshold: int) -> nx.Graph:
    names = np.unique(np.array(list(chain.from_iterable(mx.index))))
    mx = mx[combinations(names, 2)]
    return dist_matrix_to_graph(mx, threshold)


def load_matrix(fpath: str, delim: str = DELIM, name: str = "") -> pd.Series:
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

    m = np.array(matrix)[idx]

    df = pd.DataFrame(m, columns=sorted_names, index=sorted_names)
    df = df.stack().rename(name).astype(int)
    df = df.rename_axis(PAIR_IDX)
    # remove the diagonal of the matrix
    ix = [x != y for (x, y) in df.index]
    df = df[ix]
    return df


def connected_components(G: nx.Graph, node: str) -> Set[str]:
    if node not in G:
        return set()
    return nx.node_connected_component(G, node)


def clustered_together(u: str, v: str, G: nx.Graph) -> bool:
    ucc = connected_components(G, u)
    if not ucc:
        return False

    vcc = connected_components(G, v)
    if not vcc:
        return False

    return ucc == vcc


def split(arr: np.ndarray, perc: float) -> Tuple[np.ndarray, np.ndarray]:
    split_idx = int(len(arr) * perc)
    return arr[:split_idx], arr[split_idx:]


def random_split(arr: np.ndarray, perc: float) -> Tuple[np.ndarray, np.ndarray]:
    return split(rng.permuted(arr), perc)


def evaluate_clustering(
    G_true: nx.Graph, G_test: nx.Graph
) -> Tuple[float, float, float]:
    """Returns the ACP, ACR, and cluster size ratio for the clustering"""
    true_clusters = [c for c in nx.connected_components(G_true)]
    sacp_vals = []
    sacr_vals = []
    for cluster in true_clusters:
        for node in cluster:
            test_cluster = connected_components(G_test, node)
            sacr_vals.append(set_recall(cluster, test_cluster))
            sacp_vals.append(set_precision(cluster, test_cluster))

    cluster_num_ratio = nx.number_connected_components(
        G_true
    ) / nx.number_connected_components(G_test)
    return np.mean(sacp_vals), np.mean(sacr_vals), cluster_num_ratio


def run_simulation(ratio: float, threshold: int, G_true: nx.Graph, N: int):
    eval_vals = []
    for _ in range(N):
        ont, ill = random_split(SAMPLES, ratio)
        assert not set(ill) & set(ont)
        edges = [
            (s1, s2, dist)
            for (s1, s2), dist in ILLUMINA_DF[combinations(ill, 2)].iteritems()
            if dist <= threshold
        ]
        ont_threshold = THRESHOLDS[threshold]["ont"]
        edges.extend(
            [
                (s1, s2, dist)
                for (s1, s2), dist in ONT_DF[combinations(ont, 2)].iteritems()
                if dist <= ont_threshold
            ]
        )
        mixed_threshold = THRESHOLDS[threshold]["mixed"]
        edges.extend(
            [
                (s1, s2, dist)
                for (s1, s2), dist in MIXED_DF[product(ill, ont)].iteritems()
                if dist <= mixed_threshold
            ]
        )
        G_test = nx.Graph()
        G_test.add_weighted_edges_from(edges)
        true_singletons = set(SAMPLES) - set(G_true.nodes)
        test_singletons = set(SAMPLES) - set(G_test.nodes)
        xcr = excess_clustering_rate(true_singletons, test_singletons)
        acp, acr, _ = evaluate_clustering(G_true, G_test)
        eval_vals.append((acp, acr, 1 - xcr))
    df = (
        pd.DataFrame(zip(*eval_vals))
        .T.rename(columns={0: "SACP", 1: "SACR", 2: "1-XCR"})
        .melt(var_name="metric")
    )
    df["ratio"] = ratio
    return df


ILLUMINA_DF = load_matrix(snakemake.input.compass_matrix, name="illumina")
ONT_DF = load_matrix(snakemake.input.bcftools_matrix, name="ont")
MIXED_DF = load_matrix(snakemake.input.mixed_matrix, name="mixed")
RATIOS: List[int] = snakemake.params.ratios
NUM_SIMULATIONS: int = snakemake.params.num_simulations
THRESHOLDS: Dict[int, Dict[str, int]] = snakemake.params.thresholds

SAMPLES = np.unique(np.array(list(chain.from_iterable(ILLUMINA_DF.index))))


def main():
    # setup plotting canvas
    fig, axes = plt.subplots(
        nrows=snakemake.params.nrows,
        ncols=snakemake.params.ncols,
        sharex=snakemake.params.sharex,
        sharey=snakemake.params.sharey,
        squeeze=True,
    )
    kwargs = dict(
        x="ratio",
        y="value",
        hue="metric",
        split=False,
        scale="width",
        inner="quartile",
        cut=0,
        linewidth=0.5,
    )

    dfs: List[pd.DataFrame] = []
    for ax, t in zip(axes.flatten(), THRESHOLDS):
        frames = []
        G = truth_graph(ILLUMINA_DF, t)
        for r in RATIOS:
            frame = run_simulation(r, t, G, NUM_SIMULATIONS)
            frame["threshold"] = t
            frames.append(frame)
        data = pd.concat(frames)
        dfs.append(data)
        ax = sns.violinplot(data=data, ax=ax, **kwargs)
        np_t = THRESHOLDS[t]["ont"]
        mix_t = THRESHOLDS[t]["mixed"]
        ax.set_title(
            f"SNP threshold = {'/'.join(map(str, sorted({t, np_t, mix_t})))}",
            fontdict={"fontsize": 14},
        )
        ax.label_outer()
        # we only want one legend for the whole figure
        #  ax.get_legend().remove()
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=snakemake.params.legend_loc, prop=dict(size=16))
        ax.tick_params("y", labelsize=16)
        ax.set_ylabel("")

    ax.tick_params("x", labelsize=16)

    xticklabels = []
    for r in RATIOS:
        f = Fraction(1 / r - 1).limit_denominator()
        illumina_r, nanopore_r = f.as_integer_ratio()
        xticklabels.append(f"{nanopore_r}:{illumina_r}")

    ax.set_xticklabels(xticklabels)
    ax.set_xlabel(snakemake.params.xaxis_label, fontsize=16)
    fig.tight_layout()

    for fpath in snakemake.output:
        fig.savefig(fpath)


main()
