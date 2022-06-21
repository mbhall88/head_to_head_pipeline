import sys

sys.stderr = open(snakemake.log[0], "w")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from enum import Enum
from typing import Tuple
from pathlib import Path
from scipy import stats
from math import sqrt
from dataclasses import dataclass

# set aesthetics
plt.style.use(snakemake.params.style)
plt.rcParams["figure.figsize"] = snakemake.params.figsize
plt.rcParams["figure.dpi"] = snakemake.params.dpi
TOOLS = ["mykrobe"]
CONF = snakemake.params.conf_interval
FONT_SIZE = snakemake.params.font_size
GGPLOT_CM = plt.rcParams["axes.prop_cycle"].by_key()["color"]
BLUE = GGPLOT_CM[1]
MARKERSCALE = 2.0  # legend marker size
LINESTYLES = ["solid", "dotted", "dashed", "dashdot"]


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

    def f_score(self, beta: float = 1.0) -> float:
        """Harmonic mean of precision and recall.
        When beta is set to 0, you get precision. When beta is set to 1, you get the
        unweighted F-score which is the harmonic mean of precision and recall. Setting
        beta to 2 weighs recall twice as much as precision. Setting beta to 0.5 weighs
        precision twice as much as recall.
        """
        ppv = self.precision()
        tpr = self.recall()
        if ppv is None or tpr is None:
            return None
        beta2 = beta ** 2

        return ((beta2 + 1) * ppv * tpr) / ((beta2 * ppv) + tpr)

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


def round_up_to_base(x, base=10):
    return int(x + (base - x) % base)


def round_down_to_base(x, base=10):
    return int(x - (x % base))


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


def plot_a(clf_df):
    data = []

    step = 10
    tech = "nanopore"
    start = round_down_to_base(
        clf_df.query("technology==@tech")["coverage"].min(), base=step
    )
    stop = round_up_to_base(
        clf_df.query("technology==@tech")["coverage"].max(), base=step
    )
    print(stop)
    coverages = np.arange(start=start, stop=stop, step=step)
    for i, covg in enumerate(coverages[:-1]):
        next_cov = coverages[i + 1]
        subdf = clf_df.query(
            "coverage > @covg and coverage <= @next_cov and technology == @tech"
        )
        n_samples = len(set(subdf["sample"]))
        for tool in TOOLS:
            s = subdf.query("tool == @tool").value_counts(subset="classification")
            cm = ConfusionMatrix.from_series(s)
            total = sum(cm.ravel())

            data.append((covg, "FP", cm.fp / total, tech, n_samples, tool))
            data.append((covg, "TP", cm.tp / total, tech, n_samples, tool))
            data.append((covg, "FN", cm.fn / total, tech, n_samples, tool))
            data.append((covg, "TN", cm.tn / total, tech, n_samples, tool))

    tech = "illumina"
    start = round_down_to_base(
        clf_df.query("technology==@tech")["coverage"].min(), base=step
    )
    stop = round_up_to_base(
        clf_df.query("technology==@tech")["coverage"].max(), base=step
    )
    coverages = np.arange(start=start, stop=stop, step=step)
    for i, covg in enumerate(coverages[:-1]):
        next_cov = coverages[i + 1]
        subdf = clf_df.query(
            "coverage > @covg and coverage <= @next_cov and technology == @tech"
        )
        n_samples = len(set(subdf["sample"]))
        for tool in TOOLS:
            s = subdf.query("tool == @tool").value_counts(subset="classification")
            cm = ConfusionMatrix.from_series(s)
            total = sum(cm.ravel())

            data.append((covg, "FP", cm.fp / total, tech, n_samples, tool))
            data.append((covg, "TP", cm.tp / total, tech, n_samples, tool))
            data.append((covg, "FN", cm.fn / total, tech, n_samples, tool))
            data.append((covg, "TN", cm.tn / total, tech, n_samples, tool))

    summary = pd.DataFrame(
        data,
        columns=[
            "coverage",
            "classification",
            "proportion",
            "technology",
            "total",
            "tool",
        ],
    )

    with sns.plotting_context(rc={"lines.linewidth": 1.1}):
        fig, ax = plt.subplots()
        for tool, ax in zip(TOOLS, [ax]):
            x = "coverage"
            y = "proportion"
            hue = "classification"
            tech = "nanopore"
            data = summary.query("technology == @tech and tool == @tool")
            sns.pointplot(
                data=data,
                x=x,
                y=y,
                hue=hue,
                ax=ax,
                linestyles=LINESTYLES,
                alpha=0.9,
            )
            ax.set_xlabel(snakemake.params.xaxis_label, fontsize=FONT_SIZE)

            ax2 = ax.twinx()
            sns.barplot(data=data, x=x, y="total", ax=ax2, color=BLUE, alpha=0.2)
            ax2.set_ylabel("number of samples", fontsize=FONT_SIZE)
            ax.set_ylabel(ax.get_ylabel(), fontsize=FONT_SIZE)

        xlabels = ax2.get_xticklabels()
        xlabels[-1].set_text(xlabels[-1].get_text() + "+")
        ax2.set_xticklabels(xlabels)
        ax.tick_params("both", labelsize=FONT_SIZE)
        ax.set_zorder(1)
        ax2.set_zorder(1)
        ax.legend(
            loc=6, bbox_to_anchor=(0, 0.45), fontsize=FONT_SIZE, markerscale=MARKERSCALE
        )

        for fpath in snakemake.output.a:
            fig.savefig(fpath)


def plot_b(clf_df):
    data = []

    step = 10
    tech = "nanopore"
    start = round_down_to_base(
        clf_df.query("technology==@tech")["coverage"].min(), base=step
    )
    stop = round_up_to_base(
        clf_df.query("technology==@tech")["coverage"].max(), base=step
    )
    print(stop)
    coverages = np.arange(start=start, stop=stop, step=step)
    for i, covg in enumerate(coverages[:-1]):
        subdf = clf_df.query("coverage >= @covg and technology == @tech")
        n_samples = len(set(subdf["sample"]))
        for tool in TOOLS:
            s = subdf.query("tool == @tool").value_counts(subset="classification")
            cm = ConfusionMatrix.from_series(s)
            total = sum(cm.ravel())

            data.append((covg, "FP", cm.fp / total, tech, n_samples, tool))
            data.append((covg, "TP", cm.tp / total, tech, n_samples, tool))
            data.append((covg, "FN", cm.fn / total, tech, n_samples, tool))
            data.append((covg, "TN", cm.tn / total, tech, n_samples, tool))

    tech = "illumina"
    start = round_down_to_base(
        clf_df.query("technology==@tech")["coverage"].min(), base=step
    )
    stop = round_up_to_base(
        clf_df.query("technology==@tech")["coverage"].max(), base=step
    )
    coverages = np.arange(start=start, stop=stop, step=step)
    for i, covg in enumerate(coverages[:-1]):
        subdf = clf_df.query("coverage >= @covg and technology == @tech")
        n_samples = len(set(subdf["sample"]))
        for tool in TOOLS:
            s = subdf.query("tool == @tool").value_counts(subset="classification")
            cm = ConfusionMatrix.from_series(s)
            total = sum(cm.ravel())

            data.append((covg, "FP", cm.fp / total, tech, n_samples, tool))
            data.append((covg, "TP", cm.tp / total, tech, n_samples, tool))
            data.append((covg, "FN", cm.fn / total, tech, n_samples, tool))
            data.append((covg, "TN", cm.tn / total, tech, n_samples, tool))

    decum_summary = pd.DataFrame(
        data,
        columns=[
            "coverage",
            "classification",
            "proportion",
            "technology",
            "total",
            "tool",
        ],
    )

    with sns.plotting_context(rc={"lines.linewidth": 1.1}):
        fig, ax = plt.subplots()
        for tool, ax in zip(TOOLS, [ax]):
            x = "coverage"
            y = "proportion"
            hue = "classification"
            tech = "nanopore"
            data = decum_summary.query("technology == @tech and tool == @tool")
            sns.pointplot(
                data=data,
                x=x,
                y=y,
                hue=hue,
                ax=ax,
                linestyles=LINESTYLES,
                alpha=0.9,
            )
            ax.set_xlabel("")

            ax2 = ax.twinx()
            sns.barplot(data=data, x=x, y="total", ax=ax2, color=BLUE, alpha=0.2)
            ax2.set_ylabel("number of samples", fontsize=FONT_SIZE)
            ax.set_ylabel(ax.get_ylabel(), fontsize=FONT_SIZE)

        xlabels = ax2.get_xticklabels()
        xlabels[-1].set_text(xlabels[-1].get_text() + "+")
        ax2.set_xticklabels(xlabels)
        ax.set_xlabel(snakemake.params.xaxis_label, fontsize=FONT_SIZE)
        ax.tick_params("both", labelsize=FONT_SIZE)
        ax.set_zorder(1)
        ax2.set_zorder(1)
        ax.legend(
            loc=6,
            bbox_to_anchor=(0.89, 0.59),
            fontsize=FONT_SIZE,
            markerscale=MARKERSCALE,
        )

        for fpath in snakemake.output.b:
            fig.savefig(fpath)


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
    classifier = Classifier(
        unknown_is_resistant=unknown_is_resistant,
        minor_is_susceptible=minor_is_susceptible,
    )

    for ix, row in calls.iterrows():
        drug = row["drug"].lower()
        if drug in snakemake.params.ignore_drugs:
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

    plot_a(clf_df)
    plot_b(clf_df)


main()
