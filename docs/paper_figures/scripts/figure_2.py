import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import pandas as pd
from typing import List
from pathlib import Path

# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi
GGPLOT_CM = plt.rcParams["axes.prop_cycle"].by_key()["color"]
RED = GGPLOT_CM[0]
BLUE = GGPLOT_CM[1]
FONT_SIZE = snakemake.params.font_size
FILTERS = snakemake.params.filters


def extract_results(result_files: List[Path]) -> pd.DataFrame:
    columns = ["sample", "technology", "filters", "precision", "recall", "FN", "FP"]
    data = []
    for p in result_files:
        sample = p.name.split(".")[0]
        filters = "compass" if "compass" in str(p) else p.parts[-3]

        if filters != "compass":
            tech = "Nanopore"
        else:
            tech = "Illumina"

        d = (
            pd.read_csv(p)
            .query("Type=='SNP' and Filter=='PASS'")
            .to_dict(orient="records")[0]
        )
        precision = float(d["METRIC.Precision"])
        recall = float(d["METRIC.Recall"])
        fns = int(d["TRUTH.FN"])
        fps = int(d["QUERY.FP"])
        data.append((sample, tech, filters, precision, recall, fns, fps))

    return pd.DataFrame(data, columns=columns)


def main():
    dfs = [
        extract_results(list(map(Path, snakemake.input.compass_results))),
        extract_results(list(map(Path, snakemake.input.bcftools_results))),
    ]
    df = pd.concat(dfs)

    x = "filters"
    fixed_labels = ["compass", "unfiltered"]
    order = [
        *fixed_labels,
        *[label for label in sorted(df[x].unique()) if label not in fixed_labels],
    ]
    colours = [RED]
    colours.extend([BLUE] * len(FILTERS))
    boxprops = dict(linewidth=0.75, fliersize=0, showcaps=False, palette=colours)
    stripprops = dict(dodge=True, color="black", alpha=0.6)
    xlabels = ["COMPASS", *[t[1] for t in FILTERS]]

    fig, axes = plt.subplots(ncols=2, sharey=True, sharex=True)
    for ax, y in zip(
        axes.flatten(), [snakemake.params.left_fig, snakemake.params.right_fig]
    ):
        sns.boxplot(x=x, y=y, data=df, ax=ax, order=order, **boxprops)
        sns.stripplot(x=x, y=y, data=df, ax=ax, order=order, **stripprops)
        ax.set(title=y.capitalize(), xlabel="", ylabel="")
        ax.set_xlabel("")
        ax.tick_params("x", labelrotation=snakemake.params.xrotation)
        ax.tick_params("both", labelsize=FONT_SIZE)
        ax.set_yscale("logit")
        ylabs = snakemake.params.yticks
        ax.set_yticks(ylabs)
        ax.set_yticklabels(ylabs)
        _ = ax.set_xticklabels(xlabels)

    axes.flatten()[0].yaxis.set_major_formatter(PercentFormatter(1.0, decimals=2))
    red_patch = mpatches.Patch(color=RED, label="Illumina")
    blue_patch = mpatches.Patch(color=BLUE, label="Nanopore")
    _ = ax.legend(
        handles=[red_patch, blue_patch], loc="best", prop=dict(size=FONT_SIZE)
    )
    fig.tight_layout()

    for fpath in snakemake.output:
        fig.savefig(fpath)


main()
