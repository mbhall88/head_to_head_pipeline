import sys

sys.stderr = open(snakemake.log[0], "w")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi

ggplot_cm = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FONT_SIZE = snakemake.params.font_size
XCOL = "COMPASS_dist"
YCOL = "ont_dist"
BLACK = "#4c566a"
BLUE = ggplot_cm[1]
PAIR_IDX = ("sample1", "sample2")


class AsymmetrixMatrixError(Exception):
    pass


def load_matrix(fpath, delim: str = ",", name: str = "") -> pd.DataFrame:
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
ont_df = load_matrix(snakemake.input.bcftools.matrix, name=YCOL)
# merge the matrices
data = pd.concat([compass_df, ont_df], axis=1)
data = data.reset_index().rename(
    columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
)

fig, ax = plt.subplots()

# plot the full set of pairs
kwargs = dict(
    alpha=snakemake.params.alpha,
    linewidth=snakemake.params.linewidth,
    s=snakemake.params.marker_size,
)
sns.scatterplot(data=data, x=XCOL, y=YCOL, color=BLUE, ax=ax, **kwargs)

xs = data[XCOL]

# plot identity line
sns.lineplot(x=xs, y=xs, ax=ax, label="y=x", color=BLACK, linestyle="--")
ax.set_xlabel(snakemake.params.xaxis_label, fontsize=FONT_SIZE)
ax.set_ylabel(snakemake.params.yaxis_label, fontsize=FONT_SIZE)
ax.legend(loc="lower right", prop=dict(size=FONT_SIZE))

for fpath in snakemake.output:
    fig.savefig(fpath)
