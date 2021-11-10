import sys

sys.stderr = open(snakemake.log[0], "w")
from sklearn import linear_model
from typing import List, Tuple
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

plt.style.use("ggplot")
ggplot_cm = plt.rcParams["axes.prop_cycle"].by_key()["color"]
DELIM = ","
PAIR_IDX = ("sample1", "sample2")
BLACK = "#4c566a"
BLUE = ggplot_cm[1]
PURPLE = ggplot_cm[2]
GREY = ggplot_cm[3]
YELLOW = ggplot_cm[4]
RED = ggplot_cm[0]
XCOL = "xdist"
YCOL = "ydist"


def robust_regression(x: List[float], y: List[float]) -> List[float]:
    """Returns the prediction of the model. This prediction is used with X to plot
    line of best fit and get equation for that line.
    """
    X = [[v] for v in x]
    ransac = linear_model.RANSACRegressor()
    ransac.fit(X, y)
    pred = ransac.predict(X)
    return pred


def fit_model(
    x: List[float], y_pred: List[float]
) -> Tuple[float, float, float, float, float]:
    """Returns: gradient, intercept, r_value, p_value, std_err"""
    return stats.linregress(x, y_pred)


class AsymmetrixMatrixError(Exception):
    pass


if snakemake.params.mixed_dist:

    def load_matrix(fpath: str, delim: str = DELIM, name: str = "") -> pd.DataFrame:
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


else:

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
xdf = load_matrix(snakemake.input.x_matrix, name=XCOL)
ydf = load_matrix(snakemake.input.y_matrix, name=YCOL)
# merge the matrices
data = pd.concat([xdf, ydf], axis=1)
data = data.reset_index().rename(
    columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
)

fig, ax = plt.subplots(figsize=snakemake.params.figsize, dpi=snakemake.params.dpi)

# plot the full set of pairs
kwargs = snakemake.params.scatter_kws
kwargs["alpha"] /= 2
sns.scatterplot(data=data, x=XCOL, y=YCOL, color=BLUE, ax=ax, **kwargs)

xs = data[XCOL]
# plot identity line
sns.lineplot(x=xs, y=xs, ax=ax, label="y=x", color=BLACK, linestyle="--")
ax.set_xlabel(snakemake.params.xlabel)
ax.set_ylabel(snakemake.params.ylabel)
ax.legend(loc=snakemake.params.legend_loc)

if snakemake.params.plot_linear_model:
    # fit robust linear regression
    ys = data[YCOL]
    y_pred = robust_regression(xs, ys)
    slope, intercept, *_ = fit_model(xs, y_pred)
    equation = f"y={slope:.3f}x+{intercept:.3f}"
    # plot line of best fit as defined by the model
    sns.lineplot(x=xs, y=y_pred, ax=ax, label=equation, color=BLUE)


# set the lower left corner of the zoom inset (x0, y0, width, height)
axins = ax.inset_axes(bounds=snakemake.params.inset_bounds)
# sub region of the original image
threshold = snakemake.params.inset_threshold
inset_data = data.query(f"{XCOL} <= @threshold")
kwargs["alpha"] *= 2

# plot line of best fit as defined by the model
xs = inset_data[XCOL]
# plot identity line
sns.lineplot(x=xs, y=xs, ax=axins, label="y=x", color=BLACK, linestyle="--")

if snakemake.params.plot_linear_model:
    ys = inset_data[YCOL]
    y_pred = robust_regression(xs, ys)
    slope, intercept, *_ = fit_model(xs, y_pred)
    equation = f"y={slope:.3f}x+{intercept:.3f}"

    # plot line of best fit as defined by the model
    sns.lineplot(x=xs, y=y_pred, ax=axins, label=equation, color=BLUE)


# set the inset axis limits
axins.set_xlim(auto=True)
axins.set_ylim(auto=True)

# remove axis labels on inset
axins.set_xlabel("")
axins.set_ylabel("")

# styling of the inset zoom lines
rectangle_patch, connector_lines = ax.indicate_inset_zoom(
    axins, alpha=0.9, edgecolor=BLACK
)
# explicitly make the lines for the zoom window connector from the lower-left and -right corners
for i in range(len(connector_lines)):
    if i % 2 == 0:  # lower-left and lower-right
        connector_lines[i].set_visible(True)
    else:
        connector_lines[i].set_visible(False)

axins.tick_params(color=BLACK, labelcolor=BLACK)
# add a border to the inset
for spine in axins.spines.values():
    spine.set_edgecolor(BLACK)
    spine.set_alpha(0.5)

# add threshold-based partitions
xthreshold = snakemake.params.xthreshold
ythreshold = snakemake.params.ythreshold

kwargs = dict(alpha=0.9, linewidth=0.1)
# partition data into the quadrants
fns = inset_data.query(f"{YCOL}>@ythreshold and {XCOL}<=@xthreshold")
fps = inset_data.query(f"{YCOL}<=@ythreshold and {XCOL}>@xthreshold")
good_data = inset_data[~inset_data.index.isin(fns.index)]
good_data = good_data[~good_data.index.isin(fps.index)]
axins = sns.scatterplot(
    data=good_data, x=XCOL, y=YCOL, ax=axins, color=BLUE, alpha=0.5, linewidth=0.1
)
axins = sns.scatterplot(data=fns, x=XCOL, y=YCOL, ax=axins, color=RED, **kwargs)
axins = sns.scatterplot(data=fps, x=XCOL, y=YCOL, ax=axins, color=GREY, **kwargs)

# https://stackoverflow.com/questions/46961465/different-background-colour-areas-on-matplotlib-plot
ymin, ymax = axins.get_ylim()
xmin, xmax = axins.get_xlim()
axins.axhspan(
    ymin=ythreshold,
    ymax=ymax,
    xmin=0,
    xmax=((xthreshold - xmin) / (xmax - xmin)),
    facecolor=RED,
    alpha=0.2,
)
axins.axvspan(
    xmin=xthreshold,
    xmax=xmax,
    ymin=0,
    ymax=((ythreshold - ymin) / (ymax - ymin)),
    facecolor=GREY,
    alpha=0.2,
)
axins.set(ylim=(ymin, ymax), xlim=(xmin, xmax))

inset_ticklabels = list(snakemake.params.inset_ticklabels)
inset_ticklabels.append(threshold)
inset_ticklabels = list(sorted(set(inset_ticklabels)))
inset_ticklabel_fontsize = 7
axins.set_xticks(inset_ticklabels)
axins.set_xticklabels(inset_ticklabels, fontsize=inset_ticklabel_fontsize)
axins.set_yticks(inset_ticklabels)
axins.set_yticklabels(inset_ticklabels, fontsize=inset_ticklabel_fontsize)
# remove axis labels on inset
axins.set_xlabel("")
axins.set_ylabel("")

red_patch = mpatches.Patch(color=RED, label="FN")
black_patch = mpatches.Patch(color=GREY, label="FP")
axins.legend(handles=[red_patch, black_patch], loc=snakemake.params.inset_legend_loc)

fig.savefig(snakemake.output.plot)
