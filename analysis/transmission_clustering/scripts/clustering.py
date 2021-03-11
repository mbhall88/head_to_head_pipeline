import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Collection, Callable, Tuple, Set

import click
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from bokeh.io import output_file, save
from bokeh.models import (
    Range1d,
    Plot,
    Circle,
    MultiLine,
    LinearColorMapper,
    ColorBar,
    BasicTicker,
    HoverTool,
    BoxZoomTool,
    ResetTool,
    PanTool,
    WheelZoomTool,
    UndoTool,
    SaveTool,
    ColumnDataSource,
    LabelSet,
)
from bokeh.palettes import inferno, Category20
from bokeh.plotting import from_networkx
from scipy import stats

WIDTH = 1280
HEIGHT = 720
DELIM = ","
PAIR_IDX = ("sample1", "sample2")
PALETTE = Category20[20] * 10  # effectively 10 cycles of 20 colours


class AsymmetrixMatrixError(Exception):
    pass


@dataclass
class ConfusionMatrix:
    tp: int = 0
    tn: int = 0
    fp: int = 0
    fn: int = 0

    def ravel(self) -> Tuple[int, int, int, int]:
        """Return the matrix as a flattened tuple.
        The order of return is TN, FP, FN, TP
        """
        return self.tn, self.fp, self.fn, self.tp

    def as_matrix(self) -> np.ndarray:
        """Returns a 2x2 matrix [[TN, FP], [FN, TP]]"""
        return np.array([[self.tn, self.fp], [self.fn, self.tp]])

    def precision(self) -> float:
        """Also known as positive predictive value (PPV)"""
        return self.tp / (self.tp + self.fp)

    def recall(self) -> float:
        """Also known as true positive rate (TPR)"""
        return self.tp / (self.tp + self.fn)

    def fowlkes_mallows_index(self) -> float:
        """Geometric mean between precision and recall"""
        return math.sqrt(self.precision() * self.recall())

    def f_score(self, beta: float = 1.0) -> float:
        """Harmonic mean of precision and recall.
        When beta is set to 0, you get precision. When beta is set to 1, you get the
        unweighted F-score which is the harmonic mean of precision and recall. Setting
        beta to 2 weighs recall twice as much as precision. Setting beta to 0.5 weighs
        precision twice as much as recall.
        """
        ppv = self.precision()
        tpr = self.recall()
        beta2 = beta ** 2

        return ((beta2 + 1) * ppv * tpr) / ((beta2 * ppv) + tpr)

    def matthews_correlation_coefficient(self) -> float:
        """A correlation coefficient between the observed and predicted binary
        classifications.
        """
        tn, fp, fn, tp = self.ravel()
        numerator = tp * tn - fp * fn
        denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        return numerator / denominator

    @staticmethod
    def from_predictions(pred: List[bool], truth: List[bool]) -> "ConfusionMatrix":
        assert len(pred) == len(truth)
        mtx = [[0, 0], [0, 0]]
        for y_true, y_pred in zip(truth, pred):
            mtx[y_true][y_pred] += 1
        [tn, fp], [fn, tp] = mtx
        return ConfusionMatrix(tp=tp, tn=tn, fp=fp, fn=fn)


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


def plot_confusion_matrix(
    cm: ConfusionMatrix, name: str, figsize: Tuple[int, int] = (13, 8), dpi: int = 300
):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    sns.heatmap(cm.as_matrix(), annot=True, fmt="d", cmap="Blues")

    # labels, title and ticks
    ax.set_ylabel("True")
    ax.set_xlabel("Predicted")
    ax.set_title(f"{name} Clustering of Pairs - Confusion Matrix")
    ax.xaxis.set_ticklabels(["Negative", "Positive"])
    ax.yaxis.set_ticklabels(["Negative", "Positive"])
    return fig


def plot_graph(
    graph: nx.Graph, outfile: Path, threshold: int, height: int, width: int, name: str
) -> None:
    title = f"{name} clusters"
    # inline effectively allows the plot to work offline
    output_file(str(outfile), title=title, mode="inline")

    edge_attrs = {}
    cmap = inferno(threshold + 1)

    for u, v, d in graph.edges(data=True):
        edge_colour = cmap[d["weight"]]
        edge_attrs[(u, v)] = edge_colour

    nx.set_edge_attributes(graph, edge_attrs, "edge_colour")

    node_attrs = {}
    clusters = (graph.subgraph(c) for c in nx.connected_components(graph))
    for i, cluster in enumerate(clusters):
        colour = "#e5e9f0"
        for v in cluster.nodes:
            node_attrs[v] = colour

    nx.set_node_attributes(graph, node_attrs, "node_colour")

    pos = nx.nx_agraph.graphviz_layout(graph, prog="neato")
    xmin = -5
    xmax = max((x for (x, _) in pos.values())) * 1.05
    xrange = Range1d(xmin, xmax)
    ymin = -5
    ymax = max((y for (_, y) in pos.values())) * 1.05
    yrange = Range1d(ymin, ymax)

    plot = Plot(plot_width=width, plot_height=height, x_range=xrange, y_range=yrange)
    plot.title.align = "center"

    graph_renderer = from_networkx(graph, pos)
    graph_renderer.node_renderer.glyph = Circle(
        size=40, fill_color="node_colour", line_alpha=0.2
    )
    graph_renderer.edge_renderer.glyph = MultiLine(
        line_width=6, line_color="edge_colour", line_alpha=0.9
    )
    plot.renderers.append(graph_renderer)

    cmapper = LinearColorMapper(palette="Inferno256", low=0, high=threshold)
    color_bar = ColorBar(
        color_mapper=cmapper,
        major_label_text_font_size="14px",
        ticker=BasicTicker(desired_num_ticks=threshold),
        label_standoff=6,
        major_label_text_baseline="middle",
        border_line_color=None,
        location=(0, 0),
        title="SNP distance",
        title_text_align="center",
        title_standoff=10,
    )
    plot.add_layout(color_bar, "right")

    node_hover_tool = HoverTool(
        tooltips=[("sample", "@index")],
    )
    plot.add_tools(
        node_hover_tool,
        BoxZoomTool(),
        ResetTool(),
        PanTool(),
        WheelZoomTool(),
        UndoTool(),
        SaveTool(),
    )

    labels = []
    x_vals = []
    y_vals = []
    for lab, (x, y) in pos.items():
        labels.append(lab)
        x_vals.append(x)
        y_vals.append(y)

    d = {"labels": labels, "x_values": x_vals, "y_values": y_vals}
    src = ColumnDataSource(d)

    label_set = LabelSet(
        source=src,
        x="x_values",
        y="y_values",
        text="labels",
        y_offset=-5,
        text_align="center",
        text_font_size="10px",
        text_font_style="bold",
        text_color="black",
        text_font="monospace",
    )
    plot.add_layout(label_set)

    save(plot)


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

    logging.info("Reducing graphs to only those nodes with edges <= threshold")

    target_graph = dist_matrix_to_graph(target_mtx, threshold=threshold[0])
    target_mtx = target_mtx[keep_idx]
    true_labels: List[bool] = [
        clustered_together(u, v, target_graph) for u, v in target_mtx.index
    ]
    cluster_info(target_graph, name=target_name, threshold=threshold[0])
    query_graphs = []
    for mx, t, name in zip(query_mtxs, threshold[1:], samples[1:]):
        query_graph = dist_matrix_to_graph(mx, threshold=t)
        cluster_info(query_graph, name=name, threshold=t)
        query_graphs.append(query_graph)
        predicted_labels: List[bool] = [
            clustered_together(u, v, query_graph) for u, v in target_mtx.index
        ]
        conf_mtx = ConfusionMatrix.from_predictions(predicted_labels, true_labels)
        logging.info(f"Confusion matrix for {name}: {conf_mtx}")
        precision = conf_mtx.precision()
        logging.info(f"Precision: {precision}")
        recall = conf_mtx.recall()
        logging.info(f"Recall: {recall}")
        fmi = conf_mtx.fowlkes_mallows_index()
        logging.info(f"FMI: {fmi}")
        mcc = conf_mtx.matthews_correlation_coefficient()
        logging.info(f"MCC: {mcc}")
        f_score = conf_mtx.f_score()
        logging.info(f"F-score: {f_score}")
        f_beta = conf_mtx.f_score(2.0)
        logging.info(f"F-beta: {f_beta} weighing recall twice as much as precision")
        f_beta = conf_mtx.f_score(0.5)
        logging.info(f"F-beta: {f_beta} weighing precision twice as much as recall\n")
        fig = plot_confusion_matrix(conf_mtx, name=name)
        fig.savefig(outdir / f"{name}.confmatrix.png")

    for i, graph in enumerate([target_graph, *query_graphs]):
        name = samples[i]
        logging.info(f"Making cluster visualisation for {name}")
        plot_path = outdir / f"{name}.clusters.html"
        plot_graph(
            graph,
            outfile=plot_path,
            threshold=threshold[i],
            height=height,
            width=width,
            name=name,
        )
        logging.info(f"Cluster plot saved to {plot_path}")

    logging.info("Done!")


if __name__ == "__main__":
    main()
