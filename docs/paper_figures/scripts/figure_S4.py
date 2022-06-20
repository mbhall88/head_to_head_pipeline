import sys

sys.stderr = open(snakemake.log[0], "w")
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PAIR_IDX = ("sample1", "sample2")


def main():
    # set aesthetics
    plt.style.use(snakemake.params.style)
    plt.rcParams["figure.figsize"] = snakemake.params.figsize
    plt.rcParams["figure.dpi"] = snakemake.params.dpi

    matrix = []
    with open(snakemake.input.mixed_matrix) as instream:
        header = next(instream).rstrip()
        names = np.array(header.split(",")[1:])
        idx = np.argsort(names)
        sorted_names = names[idx]
        for row in map(str.rstrip, instream):
            # sort row according to the name sorting
            sorted_row = np.array(row.split(",")[1:], dtype=int)[idx]
            matrix.append(sorted_row)

    m = np.array(matrix)[idx]
    self_dists = m.diagonal()

    print(pd.DataFrame(self_dists).describe())

    largest_outlier = sorted_names[self_dists.argmax()]
    print(
        f"Largest outlier is {largest_outlier} with a self-distance of {self_dists.max()}",
    )

    fig, ax = plt.subplots()
    num_bins = max(self_dists) + 1
    sns.histplot(self_dists, bins=num_bins, discrete=True)
    ax.set_xlabel(snakemake.params.xaxis_label)
    xticks = sorted(set(Counter(self_dists).keys()))
    xticks.extend([10, 15])
    xticks.sort()
    ax.set_xticks(xticks)
    ax.set_xlim([min(xticks) - 1, max(xticks) + 1])
    ax.set_yscale("log")
    yticks = sorted(set(Counter(self_dists).values()))
    ax.set_yticks(yticks)
    _ = ax.set_yticklabels(yticks)

    for fpath in snakemake.output:
        fig.savefig(fpath)

    for i in range(1, 10):
        numerator = sum(self_dists < i)
        print(
            f"{numerator / len(self_dists):.0%} ({numerator}/{len(self_dists)}) samples below {i} distance"
        )


main()
