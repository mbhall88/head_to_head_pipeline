import subprocess
from collections import defaultdict
from pathlib import Path
from typing import TextIO, Tuple

import click
import jinja2
import pandas as pd

DEFAULT_CONTAM_WARNING = 5.0
DEFAULT_UNMAPPED_WARNING = 5.0
NORD0 = "#2e3440"
NORD11 = "#bf616a"
NORD12 = "#d08770"


class RipgrepError(Exception):
    pass


def ripgrep_search(file: Path) -> Tuple[int, int, int]:
    pattern = "\[INFO\]:\s(?P<num>\d+)\sread"
    extra_params = ["--replace", "$num", "--only-matching", "--no-line-number"]
    process = subprocess.Popen(
        ["rg", *extra_params, pattern, str(file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    exit_code = process.wait()
    if exit_code != 0:
        raise RipgrepError(
            f"Failed to execute the rg command on the file {file} due to the "
            f"following error:\n{process.stderr.read()}"
        )
    values = map(int, map(str.rstrip, process.stdout.readlines()))
    to_keep, contam, unmapped = list(values)
    return to_keep, contam, unmapped


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-a",
    "--assignment-dir",
    help="Directory containing lineage assignment CSV files.",
    type=click.Path(exists=True, file_okay=False),
    required=True,
)
@click.option(
    "-l",
    "--logs-dir",
    help="Directory containing the filter contamination script's log files.",
    type=click.Path(exists=True, file_okay=False),
    required=True,
)
@click.option(
    "-t",
    "--template",
    help="Jinja template to insert the resulting composition table into.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-o",
    "--outfile",
    type=click.File(mode="w", lazy=True),
    default="-",
    show_default=True,
    help="The filepath to write the output HTML file to.",
)
@click.option(
    "--contam-warning",
    help="Contamination above this percentage is highlighted in the HTML table.",
    default=DEFAULT_CONTAM_WARNING,
    show_default=True,
)
@click.option(
    "--unmapped-warning",
    help="Unmapped above this percentage is highlighted in the HTML table.",
    default=DEFAULT_UNMAPPED_WARNING,
    show_default=True,
)
def main(
    assignment_dir: str,
    logs_dir: str,
    template: str,
    outfile: TextIO,
    contam_warning: float,
    unmapped_warning: float,
):
    """This script generates a HTML file with a table containing information about the
    composition and lineage of each sample.
    """

    def highlight_high_contam(col: pd.Series):
        """Highlights cells in a column if their contamination level is over the
        threshold
        """
        return [
            f"background-color: {NORD12}" if val > contam_warning else "" for val in col
        ]

    def highlight_high_unmapped(col: pd.Series):
        """Highlights cells if their unmapped level is over the threshold
        """
        return [
            f"background-color: {NORD12}" if val > unmapped_warning else ""
            for val in col
        ]

    def highlight_abnormal_lineages(col: pd.Series):
        """Highlights cells if their lineage is not one of the numbered majors.
        """
        return [f"background-color: {NORD12}" if v.isalpha() else "" for v in col]

    data = defaultdict(dict)
    logfiles = Path(logs_dir).rglob("*.log")
    for file in logfiles:
        sample = file.name.split(".")[0]
        tech = file.name.split(".")[1]
        num_keep, num_contam, num_unmapped = ripgrep_search(file)
        total = sum([num_keep, num_contam, num_unmapped])
        data[sample].update(
            {
                f"{tech}_keep": num_keep,
                f"{tech}_keep%": num_keep / total,
                f"{tech}_contam": num_contam,
                f"{tech}_contam%": num_contam / total,
                f"{tech}_unmapped": num_unmapped,
                f"{tech}_unmapped%": num_unmapped / total,
                f"{tech}_total": total,
            }
        )

    assignment_files = Path(assignment_dir).rglob("*.csv")
    for file in assignment_files:
        fields = file.read_text().split("\n")[1].split(",")
        sample = file.name.split(".")[0]
        data[sample]["major_lineage"] = fields[1]
        data[sample]["full_lineage"] = fields[2]
        data[sample]["found_lineages"] = " ".join(fields[3].split(";"))

    df = pd.DataFrame(data).T
    df.index.name = "sample"

    percent_format_cols = [
        s for s in data[list(data.keys())[0]].keys() if s.endswith("%")
    ]
    df_styled = (
        df.style.apply(
            highlight_high_contam, subset=["illumina_contam%", "nanopore_contam%"]
        )
        .apply(highlight_abnormal_lineages, subset=["major_lineage"])
        .apply(
            highlight_high_unmapped, subset=["nanopore_unmapped%", "illumina_unmapped%"]
        )
        .format("{:.2%}", subset=percent_format_cols)
    )

    table_html = df_styled.render()
    template_content = Path(template).read_text()
    html = jinja2.Template(template_content).render(
        table=table_html,
        contam_warning=contam_warning,
        unmapped_warning=unmapped_warning,
    )
    outfile.write(html)
    outfile.close()


if __name__ == "__main__":
    main()
