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
    precision = float(d["Precision"]["FILT"][precision_key])
    recall = float(d["Recall"]["FILT"][recall_key])
    prg_name = p.parts[-3]
    sample = p.parts[-2]
    tool = "pandora" if "pandora" in str(p) else "compass"
    data.append((prg_name, sample, tool, precision, recall))

columns = ["prg", "sample", "tool", "precision", "recall"]
df = pd.DataFrame(data, columns=columns)

x = "tool"
y = "recall"
hue = "prg"
fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
sns.boxenplot(x, y, data=df, ax=ax, palette="Set2", hue=hue)
sns.stripplot(x, y, dodge=True, data=df, ax=ax, color="black", hue=hue)
ax.set(title=f"{y.capitalize()} for different variant callers", xlabel="variant caller")
ax.set_axis_labels(y_var=recall_key.replace("_", " "))
fig.savefig(snakemake.output.recall_plot)

y = "precision"
fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
sns.boxenplot(x, y, data=df, ax=ax, palette="Set2", hue=hue)
sns.stripplot(x, y, dodge=True, data=df, ax=ax, color="black", hue=hue)
ax.set(title=f"{y.capitalize()} for different variant callers", xlabel="variant caller")
ax.set_axis_labels(y_var=precision_key.replace("_", " "))
fig.savefig(snakemake.output.precision_plot)
