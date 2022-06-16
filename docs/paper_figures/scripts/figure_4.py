import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
import matplotlib.patches as mpatches
from matplotlib import cm
import networkx as nx
from itertools import chain
import pandas as pd
import numpy as np
from typing import Tuple, Set

PAIR_IDX = ("sample1", "sample2")


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


def XCR(G: nx.Graph, H: nx.Graph) -> Tuple[float, int, int]:
    expected_singletons = set(nx.isolates(G))
    actual_singletons = set(nx.isolates(H))
    denom = len(expected_singletons)
    numer = len(expected_singletons - actual_singletons)
    xcr = numer / denom
    return xcr, numer, denom


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
            if u not in samples:
                graph.add_node(v)
                samples.add(v)
    return graph


ILLUMINA_DIST_THRESHOLDS = snakemake.params.illumina_thresholds
NANOPORE_DIST_THRESHOLDS = snakemake.params.nanopore_thresholds
ont_thresholds = {
    ILLUMINA_DIST_THRESHOLDS[0]: NANOPORE_DIST_THRESHOLDS[0],
    ILLUMINA_DIST_THRESHOLDS[1]: NANOPORE_DIST_THRESHOLDS[1],
}
thresholds = ILLUMINA_DIST_THRESHOLDS

# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi
cmap = cm.get_cmap("tab20", 20)
PALETTE = [rgb2hex(cmap(i)) for i in range(cmap.N)]

bcftools_matrix_path = snakemake.input.bcftools_matrix
compass_matrix_path = snakemake.input.compass_matrix
clustering_metrics = []

fig, axes = plt.subplots(ncols=2, nrows=len(ILLUMINA_DIST_THRESHOLDS), squeeze=False)
for r, t in enumerate(thresholds):
    illumina_mtx = load_matrix(compass_matrix_path, name="illumina")
    G_illumina = matrix_to_graph(illumina_mtx, threshold=t, include_singletons=True)

    ont_mtx = load_matrix(bcftools_matrix_path, name="nanopore")
    G_ont = matrix_to_graph(
        ont_mtx, threshold=ont_thresholds[t], include_singletons=True
    )

    assert len(G_illumina.nodes) == len(G_ont.nodes)

    G = matrix_to_graph(illumina_mtx, threshold=t, include_singletons=True)

    node_size = 300
    edge_width = 3
    lw = 1
    font_size = 11
    fw = "bold"
    title_fontdict = dict(fontsize=12, fontweight="bold")
    align = "horizontal"
    black = "#2e3440"
    white = "#eceff4"
    red = "red"
    fc = black
    edge_colour = black
    singleton_label = "S"
    singleton_id = -1
    singleton_node_colour = white
    singleton_line_colour = red
    singleton_lw = 3
    singleton_shape = "s"  # square
    singleton_hatch = "///"
    singleton_ls = "--"
    node_shape = "o"  # circle
    colour_idx = -1
    subset_key = "cluster_id"
    node_colour_key = "node_colour"
    label_key = "label"

    node_colours = dict()
    cluster_id = dict()
    cluster_labels = dict()

    for cluster in nx.connected_components(G):
        if len(cluster) == 1:
            colour = singleton_node_colour
        else:
            colour_idx += 1
            colour = PALETTE[colour_idx]
        for v in cluster:
            node_colours[v] = colour
            if colour == singleton_node_colour:
                cluster_id[v] = singleton_id
                cluster_labels[v] = singleton_label
            else:
                cluster_id[v] = colour_idx
                cluster_labels[v] = str(colour_idx)

    max_cluster_id = colour_idx + 1

    for k, v in cluster_id.items():
        if v == singleton_id:
            cluster_id[k] = max_cluster_id

    nx.set_node_attributes(G, node_colours, node_colour_key)
    nx.set_node_attributes(G, cluster_id, subset_key)
    nx.set_node_attributes(G, cluster_labels, label_key)

    # add Illumina clustering plot
    ax = axes[r][0]

    ax.set_title(f"Illumina | threshold={t}", fontdict=title_fontdict)

    M = ont_mtx
    T = ont_thresholds
    tool = "bcftools"
    ax = axes[r][1]
    tool_t = T[t]
    H = matrix_to_graph(M, threshold=tool_t, include_singletons=True)
    xcr_data = XCR(G_illumina.copy(), H.copy())

    nx.set_node_attributes(H, node_colours, node_colour_key)
    nx.set_node_attributes(H, cluster_id, subset_key)
    nx.set_node_attributes(H, cluster_labels, label_key)

    # remove singletons that were also singletons in Illumina
    singletons_on_both = set()
    for v in list(nx.isolates(H)):
        if node_colours[v] == singleton_node_colour:  # is also singleton in Illumina
            H.remove_node(v)
            singletons_on_both.add(v)
        else:  # is not singleton in Illumina
            H.nodes[v][subset_key] = singleton_id

    # add clustered singletons back to Illumina graph and set their ID
    n = 10
    count = 0
    new_id = singleton_id + 1
    to_change = []
    for v in list(nx.isolates(G)):
        if v in singletons_on_both:
            G.remove_node(v)
            continue
        count += 1
        if count >= n:
            G.nodes[v][subset_key] = singleton_id - 1
        else:
            G.nodes[v][subset_key] = singleton_id

    for v in to_change:
        G.nodes[v][subset_key] = new_id
        for u in connected_components(G, v):
            G.nodes[u][subset_key] = new_id

    # draw Illumina clusters
    ax = axes[r][0]
    pos = nx.multipartite_layout(G, subset_key=subset_key, align=align)
    labs = {v: d[label_key] for v, d in G.nodes(data=True)}

    for shape in [singleton_shape, node_shape]:
        if shape == singleton_shape:
            vs = [k for k, v in labs.items() if v == singleton_label]
        else:
            vs = [k for k, v in labs.items() if v != singleton_label]
        ps = {v: pos[v] for v in vs}
        cols = [G.nodes[v][node_colour_key] for v in vs]

        node_ax = nx.draw_networkx_nodes(
            G,
            pos=ps,
            nodelist=vs,
            node_color=cols,
            node_shape=shape,
            edgecolors=edge_colour,
            ax=ax,
            node_size=node_size,
        )
        if shape == singleton_shape:
            node_ax.set_hatch(singleton_hatch)
            node_ax.set_edgecolor(singleton_line_colour)
            node_ax.set_linestyle(singleton_ls)

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colour, width=edge_width)
    nx.draw_networkx_labels(
        G, pos, labels=labs, ax=ax, font_size=font_size, font_weight=fw, font_color=fc
    )

    # give nodes that were singleton on Illumina, but now clustered
    # in Nanopore the cluster ID of their nanopore cluster
    clusters = nx.connected_components(H)
    for c in clusters:
        ids = {H.nodes[v][subset_key] for v in c}
        if len(ids) == 1:
            continue
        sorted_ids = sorted(ids)
        new_id = sorted_ids[0]

        for v in c:
            H.nodes[v][subset_key] = new_id

    # if there are more than n nodes in the singletons line, split over more lines
    count = 0
    new_id = singleton_id + 1
    to_change = []
    for v, d in H.nodes(data=True):
        if d[subset_key] == singleton_id:
            count += 1
            if count >= n:
                to_change.append(v)

    for v in to_change:
        H.nodes[v][subset_key] = new_id
        for u in connected_components(H, v):
            H.nodes[u][subset_key] = new_id

    counter = 0
    for cc in nx.connected_components(H):
        if all(cluster_labels[v] == singleton_label for v in cc):
            for v in cc:
                H.nodes[v][subset_key] = singleton_id - counter
            counter += 1

    ax = axes[r][1]
    pos = nx.multipartite_layout(H, subset_key=subset_key, align=align)
    labs = {v: d[label_key] for v, d in H.nodes(data=True)}
    sacp, sacr = SACP_AND_SACR(
        matrix_to_graph(illumina_mtx, threshold=t, include_singletons=False),
        matrix_to_graph(M, threshold=tool_t, include_singletons=False),
    )

    for shape in [singleton_shape, node_shape]:
        if shape == singleton_shape:
            vs = [k for k, v in labs.items() if v == singleton_label]
        else:
            vs = [k for k, v in labs.items() if v != singleton_label]
        ps = {v: pos[v] for v in vs}
        cols = [H.nodes[v][node_colour_key] for v in vs]

        node_ax = nx.draw_networkx_nodes(
            H,
            pos=ps,
            nodelist=vs,
            node_color=cols,
            node_shape=shape,
            edgecolors=edge_colour,
            ax=ax,
            node_size=node_size,
        )
        if shape == singleton_shape:
            node_ax.set_hatch(singleton_hatch)
            node_ax.set_edgecolor(singleton_line_colour)
            node_ax.set_linestyle(singleton_ls)

    nx.draw_networkx_edges(H, pos, ax=ax, edge_color=edge_colour, width=edge_width)
    nx.draw_networkx_labels(
        H, pos, labels=labs, ax=ax, font_size=font_size, font_weight=fw, font_color=fc
    )

    fs = font_size

    c = mpatches.Circle((1, 1), 0.0001, color=white)
    metrics = [
        f"SACR={sacr:.3f}",
        f"SACP={sacp:.3f}",
        f"XCR={xcr_data[0]:.3f} ({xcr_data[1]}/{xcr_data[2]})",
    ]

    clustering_metrics.append([t, tool_t, sacr, sacp, *xcr_data])

    location = dict(loc="center right", bbox_to_anchor=(0.99, 0.65))

    ax.legend(
        [c] * len(metrics),
        metrics,
        facecolor=white,
        handlelength=0.1,
        handletextpad=0.1,
        fontsize=fs,
        framealpha=0.5,
        edgecolor=black,
        frameon=True,
        **location,
    )

    ax.set_title(f"Nanopore | threshold={tool_t}", fontdict=title_fontdict)

    for i, ax in enumerate(axes.flatten()):
        fc = ax.get_facecolor()
        if i % 2:
            ax.spines["left"].set_color("black")
        else:
            ax.spines["left"].set_color(fc)

        ax.spines["top"].set_color(fc)
        ax.spines["right"].set_color(fc)
        ax.spines["bottom"].set_color(fc)

plt.tight_layout(h_pad=4.0)

for fpath in snakemake.output:
    fig.savefig(fpath)
