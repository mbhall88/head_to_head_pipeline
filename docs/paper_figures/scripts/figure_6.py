import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict
from itertools import product
from math import sqrt
from scipy import stats
from dataclasses import dataclass
from typing import Tuple
from pathlib import Path
from enum import Enum

CONF = snakemake.params.conf_interval


@dataclass
class ConfusionMatrix:
    tp: int = 0
    tn: int = 0
    fp: int = 0
    fn: int = 0

    def ravel(self) -> Tuple[int, int, int, int]:
        """Return the matrix as a flattened tuple.
        The order of return is TN, FP, FN, TP
        """
        return self.tn, self.fp, self.fn, self.tp

    def as_matrix(self) -> np.ndarray:
        """Returns a 2x2 matrix [[TN, FP], [FN, TP]]"""
        return np.array([[self.tn, self.fp], [self.fn, self.tp]])

    def num_positive(self) -> int:
        """Number of TPs and FNs - i.e. actual condition positive"""
        return self.tp + self.fn

    def num_negative(self) -> int:
        """Number of TNs and FPs - i.e. actual condition negative"""
        return self.tn + self.fp

    def ppv(self) -> Tuple[float, float, float]:
        """Also known as precision"""
        try:
            ppv = self.tp / (self.tp + self.fp)
            lwr_bound, upr_bound = confidence_interval(n_s=self.tp, n_f=self.fp)
            return ppv, lwr_bound, upr_bound
        except ZeroDivisionError:
            return [None, None, None]

    def npv(self) -> Tuple[float, float, float]:
        """Negative predictive value"""
        try:
            npv = self.tn / (self.tn + self.fn)
            lwr_bound, upr_bound = confidence_interval(n_s=self.tn, n_f=self.fn)
            return npv, lwr_bound, upr_bound
        except ZeroDivisionError:
            return [None, None, None]

    def sensitivity(self) -> float:
        """Also known as recall and true positive rate (TPR)"""
        try:
            return self.tp / self.num_positive()
        except ZeroDivisionError:
            return None

    def specificity(self) -> float:
        """Also known as selectivity and true negative rate (TNR)"""
        try:
            return self.tn / self.num_negative()
        except ZeroDivisionError:
            return None

    def fnr(self) -> Tuple[float, float, float]:
        """False negative rate or VME (very major error rate)"""
        try:
            fnr = self.fn / self.num_positive()
            lwr_bound, upr_bound = confidence_interval(n_s=self.fn, n_f=self.tp)
            return fnr, lwr_bound, upr_bound
        except ZeroDivisionError:
            return [None, None, None]

    def fpr(self) -> Tuple[float, float, float]:
        """False positive rate or ME (major error rate)"""
        try:
            fpr = self.fp / self.num_negative()
            lwr_bound, upr_bound = confidence_interval(n_s=self.fp, n_f=self.tn)
            return fpr, lwr_bound, upr_bound
        except ZeroDivisionError:
            return [None, None, None]

    @staticmethod
    def from_series(s: pd.Series) -> "ConfusionMatrix":
        tp = s.get("TP", 0)
        fp = s.get("FP", 0)
        fn = s.get("FN", 0)
        tn = s.get("TN", 0)
        return ConfusionMatrix(tp=tp, fn=fn, fp=fp, tn=tn)


def confidence_interval(n_s: int, n_f: int, conf: float = CONF) -> Tuple[float, float]:
    """Calculate the Wilson score interval.
    Equation take from https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
    n_s: Number of successes or, in the case of confusion matrix statistics, the numerator
    n_f: Number of failures or, in the case of confusion matrix statistics, the denominator minus the numerator
    conf: the confidence level. i.e. 0.95 is 95% confidence
    """
    n = n_f + n_s
    z = stats.norm.ppf(1 - (1 - conf) / 2)  # two-sided
    z2 = z ** 2
    nz2 = n + z2
    A = (n_s + (0.5 * z2)) / nz2
    B = z / nz2
    C = sqrt(((n_s * n_f) / n) + (z2 / 4))
    CI = B * C
    return A - CI, A + CI


class Prediction(Enum):
    Resistant = "R"
    Susceptible = "S"
    MinorResistance = "r"
    Unknown = "U"
    Failed = "F"

    def __str__(self) -> str:
        return self.value


class Classification(Enum):
    TruePositive = "TP"
    FalsePositive = "FP"
    TrueNegative = "TN"
    FalseNegative = "FN"

    def __str__(self) -> str:
        return self.value


class Classifier:
    def __init__(
        self,
        minor_is_susceptible: bool = False,
        unknown_is_resistant: bool = False,
        failed_is_resistant: bool = False,
    ):
        self.minor_is_susceptible = minor_is_susceptible
        self.unknown_is_resistant = unknown_is_resistant
        self.failed_is_resistant = failed_is_resistant
        self.susceptible = {Prediction.Susceptible}
        self.resistant = {Prediction.Resistant}
        if self.minor_is_susceptible:
            self.susceptible.add(Prediction.MinorResistance)
        else:
            self.resistant.add(Prediction.MinorResistance)

        if self.unknown_is_resistant:
            self.resistant.add(Prediction.Unknown)
        else:
            self.susceptible.add(Prediction.Unknown)

        if self.failed_is_resistant:
            self.resistant.add(Prediction.Failed)
        else:
            self.susceptible.add(Prediction.Failed)

    def from_predictions(
        self, y_true: Prediction, y_pred: Prediction
    ) -> Classification:
        if y_true in self.susceptible:
            expected_susceptible = True
        elif y_true in self.resistant:
            expected_susceptible = False
        else:
            raise NotImplementedError(f"Don't know how to classify {y_true} calls yet")

        if y_pred in self.susceptible:
            called_susceptible = True
        elif y_pred in self.resistant:
            called_susceptible = False
        else:
            raise NotImplementedError(f"Don't know how to classify {y_pred} calls yet")

        if expected_susceptible and called_susceptible:
            return Classification.TrueNegative
        elif expected_susceptible and not called_susceptible:
            return Classification.FalsePositive
        elif not expected_susceptible and not called_susceptible:
            return Classification.TruePositive
        else:
            return Classification.FalseNegative


# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi

ggplot_cm = plt.rcParams["axes.prop_cycle"].by_key()["color"]
red = ggplot_cm[0]
blue = ggplot_cm[1]
purple = ggplot_cm[2]
black = ggplot_cm[3]
edgecol = "black"
cmap = {"illumina": blue, "nanopore": purple}
drug_abbrev = {
    "ethambutol": "E",
    "isoniazid": "H",
    "pyrazinamide": "Z",
    "rifampicin": "R",
    "streptomycin": "S",
    "kanamycin": "Km",
    "amikacin": "Am",
    "ofloxacin": "Ofx",
    "capreomycin": "Cm",
    "moxifloxacin": "Mfx",
    "ciprofloxacin": "Cfx",
}


def main():
    covdf = pd.read_csv(snakemake.input.coverage, index_col="sample")

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

    pheno_clf = []
    minor_is_susceptible = snakemake.params.minor_is_susceptible
    unknown_is_resistant = False
    ignore_drugs = snakemake.params.ignore_drugs
    classifier = Classifier(
        unknown_is_resistant=unknown_is_resistant,
        minor_is_susceptible=minor_is_susceptible,
    )

    for ix, row in calls.iterrows():
        drug = row["drug"].lower()
        if drug in ignore_drugs:
            continue

        sample = row["sample"]
        try:
            ph = pheno.loc[(sample, drug), "phenotype"]
            if pd.isna(ph).all():
                continue
            else:
                truth = Prediction(ph[0])
        except KeyError:
            continue

        tech = row["tech"]
        if tech == "illumina":
            covg = covdf.loc[sample]["illumina_covg"]
        else:
            covg = covdf.loc[sample]["nanopore_covg"]

        pred = Prediction(row["test_call"])
        clf = classifier.from_predictions(truth, pred)

        pheno_clf.append((sample, drug, str(clf), tech, covg, row["tool"], row["site"]))
    clf_df = pd.DataFrame(
        pheno_clf,
        columns=[
            "sample",
            "drug",
            "classification",
            "technology",
            "coverage",
            "tool",
            "site",
        ],
    )

    pheno_cms = defaultdict()
    TOOLS = ["mykrobe"]
    TECHS = ["nanopore", "illumina"]
    PHENO_DRUGS = set()

    for drug, tech, tool in product(set(clf_df["drug"]), TECHS, TOOLS):
        s = clf_df.query(
            "drug == @drug and technology == @tech and tool == @tool"
        ).value_counts(subset=["classification"])
        cm = ConfusionMatrix.from_series(s)
        pheno_cms[(drug, tech, tool)] = cm
        PHENO_DRUGS.add(drug)

    PHENO_DRUGS = sorted(PHENO_DRUGS)

    fig, axes = plt.subplots(ncols=2, sharey=True)
    axR = axes.flatten()[0]
    axS = axes.flatten()[1]

    # plot details
    bar_width = 0.25
    epsilon = 0.05
    line_width = 0.5

    all_positions = []
    for i, (tech, tool) in enumerate(product(TECHS, TOOLS)):
        tps = [pheno_cms[(d, tech, tool)].tp for d in PHENO_DRUGS]
        fps = [pheno_cms[(d, tech, tool)].fp for d in PHENO_DRUGS]
        tns = [pheno_cms[(d, tech, tool)].tn for d in PHENO_DRUGS]
        fns = [pheno_cms[(d, tech, tool)].fn for d in PHENO_DRUGS]

        positions = [p + ((bar_width + epsilon) * i) for p in np.arange(len(tps))]
        all_positions.append(positions)

        colour = cmap[tech]

        # resistance bar plots
        tps_bar = axR.bar(
            positions,
            tps,
            bar_width,
            color=colour,
            edgecolor=edgecol,
            linewidth=line_width,
        )
        fns_bar = axR.bar(
            positions,
            fns,
            bar_width,
            bottom=tps,
            color=red,
            edgecolor=edgecol,
            linewidth=line_width,
        )

        # susceptible bar plots
        tns_bar = axS.bar(
            positions,
            tns,
            bar_width,
            color=colour,
            edgecolor=edgecol,
            linewidth=line_width,
        )
        fps_bar = axS.bar(
            positions,
            fps,
            bar_width,
            bottom=tns,
            color=red,
            edgecolor=edgecol,
            linewidth=line_width,
        )

        tps_bar.set_label(tech.capitalize())

    fps_bar.set_label("FP")
    fns_bar.set_label("FN")
    labels = [drug_abbrev[d] for d in PHENO_DRUGS]
    label_pos = [np.mean(ps) for ps in zip(*all_positions)]
    plt.xticks(label_pos, labels, rotation=0, fontsize=12)
    axR.set_ylabel("samples")
    axR.set_xticks(label_pos)
    axR.set_xticklabels(axS.get_xticklabels(), rotation=0, fontsize=12)
    axR.tick_params("both", labelsize=12)
    axR.set_title("Resistant")
    axS.set_title("Susceptible")

    axS.legend(loc="best", prop={"size": 12})
    leghandles, leglabels = axR.get_legend_handles_labels()

    axR.legend(leghandles, leglabels, loc="best", prop={"size": 12})
    sns.despine()
    plt.tight_layout()

    for fpath in snakemake.output:
        fig.savefig(fpath)


main()
