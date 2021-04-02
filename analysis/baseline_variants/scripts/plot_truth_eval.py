import json
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use("ggplot")

FIGSIZE = snakemake.params.figsize
DPI = snakemake.params.dpi

json_files = [Path(p) for p in snakemake.input.json_files]
recall_key = snakemake.params.recall_key
precision_key = snakemake.params.precision_key

data: List[Tuple[str, str, str, float, float]] = []
for p in json_files:
    d = json.loads(p.read_text())
    precision = float(d["Precision"][precision_key])
    recall = float(d["Recall"][recall_key])
    data.append((*p.parts[-4:-1], precision, recall))

columns = ["filters", "sample", "tool", "precision", "recall"]
df = pd.DataFrame(data, columns=columns)

x = "tool"
y = "recall"
ax_idx = 0
fig, axes = plt.subplots(ncols=2, figsize=FIGSIZE, dpi=DPI, sharex=True)
sns.boxenplot(x, y, data=df, ax=axes[ax_idx])
sns.stripplot(
    x,
    y,
    dodge=True,
    data=df,
    ax=axes[ax_idx],
    color="black",
)
axes[ax_idx].set(title=y.capitalize(), xlabel="variant caller")
axes[ax_idx].label_outer()

x = "tool"
y = "precision"
ax_idx = 1
sns.boxenplot(x, y, data=df, ax=axes[ax_idx])
sns.stripplot(
    x,
    y,
    dodge=True,
    data=df,
    ax=axes[ax_idx],
    color="black",
)
axes[ax_idx].set(title=y.capitalize(), xlabel="variant caller")
axes[ax_idx].label_outer()
fig.savefig(snakemake.output.plot)
