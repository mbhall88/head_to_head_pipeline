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


class Classifier(Enum):
    TruePositive = "TP"
    FalsePositive = "FP"
    TrueNegative = "TN"
    FalseNegative = "FN"

    def __str__(self) -> str:
        return self.value

    @staticmethod
    def from_predictions(y_true: Prediction, y_pred: Prediction) -> "Classifier":
        if y_true is Prediction.Susceptible:
            return (
                Classifier.TrueNegative
                if y_pred is Prediction.Susceptible
                else Classifier.FalsePositive
            )
        elif y_true is Prediction.Resistant:
            return (
                Classifier.TruePositive
                if y_pred is Prediction.Resistant
                else Classifier.FalseNegative
            )
        else:
            raise NotImplementedError(
                "Don't know how to classify minor resistant calls yet"
            )


def load_susceptibility(stream: TextIO) -> dict:
    """Extract the susceptibility info from the JSON"""
    data = json.load(stream)
    try:
        return data[next(iter(data.keys()))]["susceptibility"]
    except (KeyError, TypeError):
        return data["susceptibility"]


def extract_calls(data: dict, convert_minor: bool = False) -> Dict[str, Prediction]:
    """Takes a susceptibility dictionary and returns the calls for each drug"""
    calls = dict()

    for drug in data:
        call = data[drug]["predict"]
        prediction = Prediction(call.upper()) if convert_minor else Prediction(call)
        calls[drug] = prediction

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
    "--minor-is-major",
    help="Minor resistance 'r' should be treated as major resistance 'R'",
    is_flag=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    truth_istream: TextIO,
    test_istream: TextIO,
    output: TextIO,
    verbose: bool,
    minor_is_major: bool,
):
    """Concordance of drug resistant predictions"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    truth_susceptibility = load_susceptibility(truth_istream)
    true_calls = extract_calls(truth_susceptibility, convert_minor=minor_is_major)

    test_susceptibility = load_susceptibility(test_istream)
    test_calls = extract_calls(test_susceptibility, convert_minor=minor_is_major)

    if len(true_calls.keys() - test_calls.keys()) != 0:
        raise KeyError("Drug names aren't consistent between the two datasets")

    logging.info("Classifying prediction concordance...")
    print("drug,classification,true_call,test_call", file=output)
    for drug, y_true in true_calls.items():
        y_pred = test_calls[drug]
        classification = Classifier.from_predictions(y_true, y_pred)
        print(f"{drug},{classification},{y_true},{y_pred}", file=output)

    logging.info("Done!")


if __name__ == "__main__":
    main()
