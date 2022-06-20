import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np


# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi

ggplot_cm = plt.rcParams["axes.prop_cycle"].by_key()["color"]

DELIM = ","
PAIR_IDX = ("sample1", "sample2")
BLACK = "#4c566a"
BLUE = ggplot_cm[1]
PURPLE = ggplot_cm[2]
YELLOW = ggplot_cm[4]
RED = ggplot_cm[0]
XCOL = "COMPASS_dist"
YCOL = "ont_dist"
marker_size = snakemake.params.marker_size
false_marker = snakemake.params.false_marker
ILLUMINA_DIST_THRESHOLDS = snakemake.params.illumina_thresholds
NANOPORE_DIST_THRESHOLDS = snakemake.params.nanopore_thresholds
FONT_SIZE = snakemake.params.font_size


class AsymmetrixMatrixError(Exception):
    pass


def load_matrix(fpath, delim: str = DELIM, name: str = "") -> pd.DataFrame:
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


# load the data
compass_df = load_matrix(snakemake.input.compass_matrix, name=XCOL)
ont_df = load_matrix(snakemake.input.bcftools_matrix, name=YCOL)
# merge the matrices
data = pd.concat([compass_df, ont_df], axis=1)
data = data.reset_index().rename(
    columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
)

threshold = snakemake.params.close_threshold
inset_data = data.query(f"{XCOL} <= @threshold")

fig, ax = plt.subplots()

ax.set_xlabel(snakemake.params.xaxis_label, fontsize=FONT_SIZE)
ax.set_ylabel(snakemake.params.yaxis_label, fontsize=FONT_SIZE)

xs = inset_data[XCOL]

# plot identity line
sns.lineplot(x=xs, y=xs, ax=ax, label="y=x", color=BLACK, linestyle="--", linewidth=1.0)

# set the inset axis limits
ax.set_xlim(auto=True)
ax.set_ylim(auto=True)

ax.tick_params(color=BLACK, labelcolor=BLACK)

# add threshold-based partitions
illumina_t = ILLUMINA_DIST_THRESHOLDS[1]
nanopore_t = NANOPORE_DIST_THRESHOLDS[1]
lower_illumina_t = ILLUMINA_DIST_THRESHOLDS[0]
lower_nanopore_t = NANOPORE_DIST_THRESHOLDS[0]

kwargs = dict(alpha=0.9, linewidth=0.1, s=marker_size)
# partition data into the quadrants
ont_fns = inset_data.query(f"{YCOL}>@nanopore_t and {XCOL}<=@illumina_t")
ont_fps = inset_data.query(f"{YCOL}<=@nanopore_t and {XCOL}>@illumina_t")
good_data = inset_data[~inset_data.index.isin(ont_fns.index)]
good_data = good_data[~good_data.index.isin(ont_fps.index)]
ax = sns.scatterplot(
    data=good_data,
    x=XCOL,
    y=YCOL,
    ax=ax,
    color=BLUE,
    alpha=0.5,
    linewidth=0.1,
    s=marker_size,
)
ax = sns.scatterplot(data=ont_fns, x=XCOL, y=YCOL, ax=ax, color=RED, **kwargs)
ax = sns.scatterplot(data=ont_fps, x=XCOL, y=YCOL, ax=ax, color=ggplot_cm[3], **kwargs)

# https://stackoverflow.com/questions/46961465/different-background-colour-areas-on-matplotlib-plot
ymin, ymax = ax.get_ylim()
xmin, xmax = ax.get_xlim()
ax.axhspan(
    ymin=nanopore_t,
    ymax=ymax,
    xmin=0,
    xmax=((illumina_t - xmin) / (xmax - xmin)),
    facecolor=RED,
    alpha=0.2,
)
ax.axvspan(
    xmin=illumina_t,
    xmax=xmax,
    ymin=0,
    ymax=((nanopore_t - ymin) / (ymax - ymin)),
    facecolor=ggplot_cm[3],
    alpha=0.2,
)
ax.set(ylim=(ymin, ymax), xlim=(xmin, xmax))

inset_ticklabels = sorted(
    {0, 2, lower_nanopore_t, lower_illumina_t, illumina_t, nanopore_t, threshold}
)
inset_fs = FONT_SIZE
ax.set_xticks(inset_ticklabels)
ax.set_xticklabels(inset_ticklabels, fontsize=inset_fs)
inset_ticklabels.append(inset_data.max()[YCOL])
ax.set_yticks(inset_ticklabels)
ax.set_yticklabels(inset_ticklabels, fontsize=inset_fs)

ax.axhspan(
    ymin=lower_nanopore_t,
    ymax=nanopore_t,
    xmin=0,
    xmax=((lower_illumina_t - xmin) / (xmax - xmin)),
    facecolor=RED,
    alpha=0.3,
    hatch="//",
    edgecolor="black",
)
ax.axvspan(
    xmin=lower_illumina_t,
    xmax=illumina_t,
    ymin=0,
    ymax=((lower_nanopore_t - ymin) / (ymax - ymin)),
    facecolor=ggplot_cm[3],
    alpha=0.3,
    hatch="//",
    edgecolor="black",
)
ont_fns = good_data.query(f"{YCOL}>@lower_nanopore_t and {XCOL}<=@lower_illumina_t")
ont_fps = good_data.query(f"{YCOL}<=@lower_nanopore_t and {XCOL}>@lower_illumina_t")
good_data = good_data[~good_data.index.isin(ont_fns.index)]
good_data = good_data[~good_data.index.isin(ont_fps.index)]
ax = sns.scatterplot(
    data=good_data,
    x=XCOL,
    y=YCOL,
    ax=ax,
    color=BLUE,
    alpha=0.5,
    linewidth=0.1,
    s=marker_size,
)
ax = sns.scatterplot(
    data=ont_fns, x=XCOL, y=YCOL, ax=ax, color=BLUE, marker=false_marker, **kwargs
)
ax = sns.scatterplot(
    data=ont_fps, x=XCOL, y=YCOL, ax=ax, color=BLUE, marker=false_marker, **kwargs
)

red_patch = mpatches.Patch(color=RED, label="FN")
black_patch = mpatches.Patch(color=ggplot_cm[3], label="FP")
ax.legend(handles=[red_patch, black_patch], loc="upper left", prop=dict(size=FONT_SIZE))

for fpath in snakemake.output:
    fig.savefig(fpath)
