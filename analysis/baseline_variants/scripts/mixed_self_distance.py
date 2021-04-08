import sys

sys.stderr = open(snakemake.log[0], "w")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use("ggplot")
PAIR_IDX = ("sample1", "sample2")
FIGSIZE = (13, 8)
DPI = 300
DELIM = ","

matrix = []
with open(snakemake.input.matrix) as instream:
    header = next(instream).rstrip()
    names = np.array(header.split(DELIM)[1:])
    idx = np.argsort(names)
    sorted_names = names[idx]
    for row in map(str.rstrip, instream):
        # sort row according to the name sorting
        sorted_row = np.array(row.split(DELIM)[1:], dtype=np.int)[idx]
        matrix.append(sorted_row)

m = np.array(matrix)[idx]
n_samples = len(sorted_names)

df = pd.DataFrame(m, columns=sorted_names, index=sorted_names)
df = df.stack().astype(int)
df = df.rename_axis(PAIR_IDX)
self_dists = m.diagonal()

print("Summary statistics of distribution", file=sys.stderr)
print(pd.DataFrame(self_dists).describe(), file=sys.stderr)

largest_outlier = sorted_names[self_dists.argmax()]
print(
    f"Largest outlier is {largest_outlier} with a self-distance of {self_dists.max()}",
    file=sys.stderr,
)

fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
num_bins = max(self_dists) + 1
sns.histplot(self_dists, bins=num_bins, discrete=True)
ax.set_xlabel("SNP distance")
xticks = np.append(ax.get_xticks(), range(10))
xticks.sort()
ax.set_xticks(xticks[1:])
ax.set_xlim([-1, max(xticks)])
fig.savefig(snakemake.output.plot)
