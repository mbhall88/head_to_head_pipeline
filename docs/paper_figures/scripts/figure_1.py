import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import upsetplot


# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi

frames = []
for p in map(Path, snakemake.input.concordance):
    sample, tool = p.name.split(".")[0:2]
    tech = p.parts[-3]
    site = p.parts[-2]
    table = pd.read_csv(p)
    table["sample"] = sample
    table["tool"] = tool
    table["tech"] = tech
    table["site"] = site
    frames.append(table)

calls = pd.concat(frames)
calls.reset_index(drop=True, inplace=True)
valid_samples = set(calls["sample"])
DRUGS = sorted(set(calls["drug"]))

pheno = pd.read_csv(snakemake.input.phenosheet).melt(
    id_vars=["sample"], var_name="drug", value_name="phenotype"
)
arr = []
for r in pheno["phenotype"]:
    if pd.isna(r):
        arr.append(r)
    elif r.upper() in ("R", "S"):
        arr.append(r.upper())
    else:
        arr.append(None)
pheno["phenotype"] = arr
pheno = pheno.loc[pheno["sample"].isin(valid_samples)]
pheno.set_index(["sample", "drug"], drop=False, inplace=True)
pheno.sort_index(inplace=True)

d = {}
exclude = snakemake.params.exclude
samples_with_pheno = []
for drug in map(str.lower, DRUGS):
    if drug in exclude:
        continue
    samples = list(pheno.query("drug == @drug").dropna()["sample"])
    samples_with_pheno.extend(samples)
    if samples:
        d[drug.upper()] = samples

upset_data = upsetplot.from_contents(d)

fig, ax = plt.subplots()
p = upsetplot.plot(
    upset_data,
    fig=fig,
    sort_by="cardinality",
    orientation="horizontal",
    show_counts=True,
)
p["intersections"].set_ylabel("count")
ax.axis("off")
p["totals"].set_xticks(snakemake.params.xticks)

for fpath in snakemake.output:
    fig.savefig(fpath)
