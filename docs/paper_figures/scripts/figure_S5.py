import sys

sys.stderr = open(snakemake.log[0], "w")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.scale import scale_factory
import matplotlib.patches as mpatches


GGPLOT_CM = plt.rcParams["axes.prop_cycle"].by_key()["color"]
FONT_SIZE = snakemake.params.font_size
XCOL = "COMPASS_dist"
YCOL = "mixed_dist"
BLACK = "#4c566a"
BLUE = GGPLOT_CM[1]
RED = GGPLOT_CM[0]
PAIR_IDX = ("sample1", "sample2")
ALPHA = snakemake.params.alpha
LW = snakemake.params.line_width
ILLUMINA_DIST_THRESHOLDS = snakemake.params.illumina_thresholds
MIXED_DIST_THRESHOLDS = snakemake.params.mixed_thresholds


def load_mixed_matrix(fpath: str, delim: str = ",", name: str = "") -> pd.DataFrame:
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


def main():
    # load the data
    compass_df = load_mixed_matrix(snakemake.input.compass_matrix, name=XCOL)
    mixed_df = load_mixed_matrix(snakemake.input.mixed_matrix, name=YCOL)
    # merge the matrices
    data = pd.concat([compass_df, mixed_df], axis=1)
    data = data.reset_index().rename(
        columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
    )

    # set aesthetics
    plt.style.use(snakemake.params.style)
    plt.rcParams["figure.figsize"] = snakemake.params.figsize
    plt.rcParams["figure.dpi"] = snakemake.params.dpi

    fig, ax = plt.subplots()

    # plot the full set of pairs
    kwargs = dict(alpha=ALPHA, linewidth=LW)
    sns.scatterplot(data=data, x=XCOL, y=YCOL, color=BLUE, ax=ax, **kwargs)

    xs = data[XCOL]

    # plot identity line
    sns.lineplot(x=xs, y=xs, ax=ax, label="y=x", color=BLACK, linestyle="--")
    ax.set_xlabel(snakemake.params.xaxis_label, fontsize=FONT_SIZE)
    ax.set_ylabel(snakemake.params.yaxis_label, fontsize=FONT_SIZE)
    ax.legend(loc="lower right", prop=dict(size=FONT_SIZE))

    # set the lower left corner of the zoom inset (x0, y0, width, height)
    axins = ax.inset_axes(bounds=[0.05, 0.52, 0.5, 0.48])
    # sub region of the original image
    threshold = snakemake.params.close_threshold
    inset_data = data.query(f"{XCOL} <= @threshold")

    kwargs["alpha"] *= 2

    xs = inset_data[XCOL]

    # plot identity line
    sns.lineplot(x=xs, y=xs, ax=axins, label="y=x", color=BLACK, linestyle="--")

    # set the inset axis limits
    axins.set_xlim(auto=True)
    axins.set_ylim(auto=True)
    axins.set_yscale("symlog")

    symlog = scale_factory("symlog", axins).get_transform().transform

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
    illumina_t = ILLUMINA_DIST_THRESHOLDS[1]
    mixed_t = MIXED_DIST_THRESHOLDS[1]
    lower_illumina_t = ILLUMINA_DIST_THRESHOLDS[0]
    lower_mixed_t = MIXED_DIST_THRESHOLDS[0]

    kwargs = dict(alpha=0.9, linewidth=0.1)
    # partition data into the quadrants
    mixed_fns = inset_data.query(f"{YCOL}>@mixed_t and {XCOL}<=@illumina_t")
    mixed_fps = inset_data.query(f"{YCOL}<=@mixed_t and {XCOL}>@illumina_t")
    good_data = inset_data[~inset_data.index.isin(mixed_fns.index)]
    good_data = good_data[~good_data.index.isin(mixed_fps.index)]
    axins = sns.scatterplot(
        data=good_data, x=XCOL, y=YCOL, ax=axins, color=BLUE, alpha=0.5, linewidth=0.1
    )
    axins = sns.scatterplot(
        data=mixed_fns, x=XCOL, y=YCOL, ax=axins, color=RED, **kwargs
    )
    axins = sns.scatterplot(
        data=mixed_fps, x=XCOL, y=YCOL, ax=axins, color=GGPLOT_CM[3], **kwargs
    )

    # https://stackoverflow.com/questions/46961465/different-background-colour-areas-on-matplotlib-plot
    ymin, ymax = axins.get_ylim()
    xmin, xmax = axins.get_xlim()
    axins.axhspan(
        ymin=mixed_t,
        ymax=ymax,
        xmin=0,
        xmax=((illumina_t - xmin) / (xmax - xmin)),
        facecolor=RED,
        alpha=0.2,
    )
    axins.axvspan(
        xmin=illumina_t,
        xmax=xmax,
        ymin=0,
        ymax=((symlog(mixed_t) - symlog(ymin)) / (symlog(ymax) - symlog(ymin))),
        facecolor=GGPLOT_CM[3],
        alpha=0.2,
    )
    # fig.savefig(snakemake.output.plot)
    axins.set(ylim=(ymin, ymax), xlim=(xmin, xmax))

    inset_ticklabels = sorted(
        {0, 2, lower_mixed_t, lower_illumina_t, illumina_t, mixed_t, threshold}
    )
    inset_fs = FONT_SIZE - 5
    axins.set_xticks(inset_ticklabels)
    axins.set_xticklabels(inset_ticklabels, fontsize=inset_fs)
    inset_ticklabels.append(inset_data.max()[YCOL])
    axins.set_yticks(inset_ticklabels)
    axins.set_yticklabels(inset_ticklabels, fontsize=inset_fs)

    axins.axhspan(
        ymin=lower_mixed_t,
        ymax=mixed_t,
        xmin=0,
        xmax=((lower_illumina_t - xmin) / (xmax - xmin)),
        facecolor=RED,
        alpha=0.1,
        hatch="//",
        edgecolor="black",
    )
    axins.axvspan(
        xmin=lower_illumina_t,
        xmax=illumina_t,
        ymin=0,
        ymax=((symlog(lower_mixed_t) - symlog(ymin)) / (symlog(ymax) - symlog(ymin))),
        facecolor=GGPLOT_CM[3],
        alpha=0.1,
        hatch="//",
        edgecolor="black",
    )
    mixed_fns = good_data.query(f"{YCOL}>@lower_mixed_t and {XCOL}<=@lower_illumina_t")
    mixed_fps = good_data.query(f"{YCOL}<=@lower_mixed_t and {XCOL}>@lower_illumina_t")
    good_data = good_data[~good_data.index.isin(mixed_fns.index)]
    good_data = good_data[~good_data.index.isin(mixed_fps.index)]
    axins = sns.scatterplot(
        data=good_data, x=XCOL, y=YCOL, ax=axins, color=BLUE, alpha=0.5, linewidth=0.1
    )
    axins = sns.scatterplot(
        data=mixed_fns, x=XCOL, y=YCOL, ax=axins, color=BLUE, marker="s", **kwargs
    )
    axins = sns.scatterplot(
        data=mixed_fps, x=XCOL, y=YCOL, ax=axins, color=BLUE, marker="s", **kwargs
    )

    # remove axis labels on inset
    axins.set_xlabel("")
    axins.set_ylabel("")

    red_patch = mpatches.Patch(color=RED, label="FN")
    black_patch = mpatches.Patch(color=GGPLOT_CM[3], label="FP")
    axins.legend(
        handles=[red_patch, black_patch], loc="upper left", prop=dict(size=inset_fs)
    )

    for fpath in snakemake.output:
        fig.savefig(fpath)


main()
