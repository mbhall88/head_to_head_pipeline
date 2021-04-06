import sys

sys.stderr = open(snakemake.log[0], "w")
from sklearn import linear_model
from typing import List, Tuple
from scipy import stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use("ggplot")
DELIM = ","
PAIR_IDX = ("sample1", "sample2")
BLACK = "#4c566a"
BLUE = "#5e81ac"
XCOL = "COMPASS_dist"
YCOL = "ont_dist"


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
ont_df = load_matrix(snakemake.input.ont_matrix, name=YCOL)
# merge the matrices
data = pd.concat([compass_df, ont_df], axis=1)
data = data.reset_index().rename(
    columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
)

fig, ax = plt.subplots(figsize=snakemake.params.figsize, dpi=snakemake.params.dpi)

# plot the full set of pairs
kwargs = snakemake.params.scatter_kws
kwargs["alpha"] /= 2
sns.scatterplot(data=data, x=XCOL, y=YCOL, ax=ax, **kwargs)

# fit robust linear regression
xs = data[XCOL]
ys = data[YCOL]
y_pred = robust_regression(xs, ys)
slope, intercept, *_ = fit_model(xs, y_pred)
equation = f"y={slope:.3f}x+{intercept:.3f}"

# plot line of best fit as defined by the model
sns.lineplot(x=xs, y=y_pred, ax=ax, label=equation, color=BLUE)
# plot identity line
sns.lineplot(x=xs, y=xs, ax=ax, label="y=x", color=BLACK, linestyle="--")
ax.set_xlabel(snakemake.params.xlabel)
ax.set_ylabel(snakemake.params.ylabel)
ax.legend(loc=snakemake.params.legend_loc)

# set the lower left corner of the zoom inset (x0, y0, width, height)
axins = ax.inset_axes(bounds=snakemake.params.inset_bounds)
# sub region of the original image
threshold = snakemake.params.inset_threshold
inset_data = data.query(f"{XCOL} <= @threshold")
kwargs["alpha"] *= 2
axins = sns.scatterplot(data=inset_data, x=XCOL, y=YCOL, ax=axins, **kwargs)

# plot line of best fit as defined by the model
xs = inset_data[XCOL]
ys = inset_data[YCOL]
y_pred = robust_regression(xs, ys)
slope, intercept, *_ = fit_model(xs, y_pred)
equation = f"y={slope:.3f}x+{intercept:.3f}"

# plot line of best fit as defined by the model
sns.lineplot(x=xs, y=y_pred, ax=axins, label=equation, color=BLUE)

# plot identity line
sns.lineplot(x=xs, y=xs, ax=axins, label="y=x", color=BLACK, linestyle="--")

# set the inset axis limits
axins.set_xlim(auto=True)
axins.set_ylim(auto=True)

# remove axis labels on inset
axins.set_xlabel("")
axins.set_ylabel("")

# styling of the inset zoom lines
ax.indicate_inset_zoom(axins, alpha=0.8, edgecolor=BLACK)
axins.tick_params(color=BLACK, labelcolor=BLACK)
# add a border to the inset
for spine in axins.spines.values():
    spine.set_edgecolor(BLACK)
    spine.set_alpha(0.5)

axins.legend(loc=snakemake.params.inset_legend_loc)

fig.savefig(snakemake.output.plot)
