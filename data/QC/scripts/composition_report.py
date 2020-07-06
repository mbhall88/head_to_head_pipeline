from collections import defaultdict
from pathlib import Path
from typing import TextIO, Dict

import click
import pandas as pd


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
    "-k",
    "--krona-dir",
    help="Directory containing krona input TSV files.",
    type=click.Path(exists=True, file_okay=False),
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
def main(
    assignment_dir: str, krona_dir: str, outfile: TextIO,
):
    """This script generates a HTML file with a table containing information about the
    composition and lineage of each sample.
    """
    krona_files = Path(krona_dir).rglob("*.tsv")
    for file in krona_files:
        sample = file.name.split(".")[0]
        with file.open() as istream:
            counts = count_krona(istream)
            #todo - figure out how to combine all these into a df



def count_krona(istream: TextIO) -> Dict[str, int]:
    counts = defaultdict(int)
    for line in map(str.rstrip, istream):
        count, classification = line.split("\t")[:2]
        counts[classification] += int(count)

    return counts

if __name__ == "__main__":
    main()
