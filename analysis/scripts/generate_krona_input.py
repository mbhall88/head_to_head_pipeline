"""This script generates the text file input required to make a krona plot.
It takes a bam file of reads mapped to a reference/database. Additionally, it
requires a metadata file (tsv) that has information about each reference in the
database. Column 1 is the category, column 2 is whether the reference is
contamination, column 3 is the accession ID for the reference."""

from collections import Counter, defaultdict

import pandas as pd
import pysam


def get_ena_accession(accession_string):
    """Extract the accession ID from a header string"""
    return accession_string.split("|")[1]


def get_organism_assignments(bamfile: str, metadata):
    """Count the assignments for each read in a bam file.
    :param bamfile: Path to bam file.
    :param metadata: Pandas dataframe of the contamination database info.
    :return: A dictionary with the keys being organisms and the values being
    Counters with keys being accession IDs and values being counts.
    """
    organism_assignments = defaultdict(Counter)
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for record in bam:
            ref = record.reference_name
            if ref is not None:
                organism = metadata.loc[ref]["organism"]
                organism_assignments[organism][ref] += 1
            else:
                organism_assignments["UNMAPPED"]["NA"] += 1
    return organism_assignments


def run():
    metadata = pd.read_table(
        snakemake.input.metadata,
        header=None,
        names=["organism", "contamination", "accession"],
        index_col="accession",
    )
    organism_assignments = get_organism_assignments(snakemake.input.bam, metadata)

    with open(snakemake.output[0], "w") as f_out:
        for category, assignments in organism_assignments.items():
            for accession, count in assignments.items():
                if accession == "NA":
                    print("{count}\tUnmapped".format(count=count), file=f_out)
                else:
                    print("{count}\t{category}\t{accession}".format(
                        count=count, category=category, accession=accession), file=f_out)


# ========================================================================
# ========================================================================
run()
