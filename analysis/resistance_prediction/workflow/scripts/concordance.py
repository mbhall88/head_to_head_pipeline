import json
import logging
from enum import Enum
from typing import TextIO, Dict

import click


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


def load_susceptibility(stream: TextIO) -> dict:
    """Extract the susceptibility info from the JSON"""
    data = json.load(stream)
    try:
        return data[next(iter(data.keys()))]["susceptibility"]
    except (KeyError, TypeError):
        return data["susceptibility"]


def extract_calls(data: dict) -> Dict[str, Prediction]:
    """Takes a susceptibility dictionary and returns the calls for each drug"""
    calls = dict()

    for drug in data:
        call = data[drug]["predict"]
        calls[drug] = Prediction(call)

    return calls


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-a",
    "--true",
    "truth_istream",
    help="Resistance classification JSON that is considered truth",
    type=click.File(),
    required=True,
)
@click.option(
    "-b",
    "--test",
    "test_istream",
    help="Resistance classification JSON that is being tested for concordance",
    type=click.File(),
    required=True,
)
@click.option(
    "-o",
    "--output",
    help="File to write output to [default: stdout]",
    type=click.File(mode="w"),
)
@click.option(
    "-r",
    "--minor-is-susceptible",
    help="Minor resistance 'r' should be treated as susceptible 'S'",
    is_flag=True,
)
@click.option(
    "-U",
    "--unknown-is-resistant",
    help="Unknown calls are considered resistant. By default they're considered susceptible",
    is_flag=True,
)
@click.option(
    "-F",
    "--failed-is-resistant",
    help="Failed calls are considered resistant. By default they're considered susceptible",
    is_flag=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    truth_istream: TextIO,
    test_istream: TextIO,
    output: TextIO,
    verbose: bool,
    minor_is_susceptible: bool,
    unknown_is_resistant: bool,
):
    """Concordance of drug resistant predictions"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    truth_susceptibility = load_susceptibility(truth_istream)
    true_calls = extract_calls(truth_susceptibility)

    test_susceptibility = load_susceptibility(test_istream)
    test_calls = extract_calls(test_susceptibility)

    if len(true_calls.keys() - test_calls.keys()) != 0:
        raise KeyError("Drug names aren't consistent between the two datasets")

    logging.info("Classifying prediction concordance...")
    print("drug,classification,true_call,test_call", file=output)
    classifier = Classifier(
        unknown_is_resistant=unknown_is_resistant,
        minor_is_susceptible=minor_is_susceptible,
        failed_is_resistant=failed_is_resistant,
    )
    for drug, y_true in true_calls.items():
        y_pred = test_calls[drug]
        classification = classifier.from_predictions(y_true, y_pred)
        print(f"{drug},{classification},{y_true},{y_pred}", file=output)

    logging.info("Done!")


if __name__ == "__main__":
    main()
