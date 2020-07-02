from collections import Counter, defaultdict
from typing import TextIO, Dict

import click
import pandas as pd
import pysam

DELIM = "\t"


def get_ena_accession(accession: str) -> str:
    return accession.split("|")[1]


def get_organism_assignments(
    alnfile: str, metadata: pd.DataFrame, ignore_secondary: bool = True
) -> Dict[str, Counter]:
    organism_assignments = defaultdict(Counter)
    with pysam.AlignmentFile(alnfile) as bam:
        for record in bam:
            if record.is_secondary and ignore_secondary:
                continue
            if record.is_unmapped:
                organism_assignments["UNMAPPED"]["NA"] += 1
                continue
            ref = record.reference_name
            organism = metadata.loc[ref]["organism"]
            organism_assignments[organism][ref] += 1

    return organism_assignments


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--samfile",
    help="{B,CR,S}AM file of reads mapped to decontamination database.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-m",
    "--metadata",
    help=(
        "TSV file containing information about each reference in the database. Column "
        "1 is the category, column 2 is whether the reference is contamination, column "
        "3 is the accession ID for the reference."
    ),
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-o",
    "--outfile",
    type=click.File(mode="w", lazy=True),
    default="-",
    show_default=True,
    help=(
        "The filepath to write the output to. This file is formatted for input to "
        "krona."
    ),
)
@click.option(
    "--ignore-secondary/--include-secondary",
    help="Ignore organism assignments for secondary alignments?",
    default=True,
    show_default=True,
)
def main(
    samfile: str, metadata: str, outfile: TextIO, ignore_secondary: bool,
):
    """This script generates the text file input required to make a krona plot.
    """
    metadata_df = pd.read_table(
        metadata,
        header=None,
        names=["organism", "contamination", "accession"],
        index_col="accession",
    )
    organism_assignments = get_organism_assignments(
        samfile, metadata_df, ignore_secondary
    )

    for category, assignments in organism_assignments.items():
        for accession, count in assignments.items():
            if accession == "NA":
                content = DELIM.join([str(count), "Unmapped"])
            else:
                content = DELIM.join([str(count), category, accession])
            print(content, file=outfile)


if __name__ == "__main__":
    main()
