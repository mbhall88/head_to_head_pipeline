from collections import Counter, defaultdict

import pandas as pd
import pysam


def get_ena_accession(accession_string):
    """Extract the accession ID from a header string"""
    return accession_string.split("|")[1]


def get_organism_assignments(bamfile: str, df):
    """Count the assignments for each read in a bam file.
    :param bamfile: Path to bam file.
    :param df: Pandas dataframe of the contamination database info.
    :return: A dictionary with the keys being organisms and the values being
    Counters with keys being accession IDs and values being counts.
    """
    organism_assignments = defaultdict(Counter)
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        ref = record.reference_name
        for record in bam:
            if ref is not None:
                organism = df.loc[ref]["organism"]
                organism_assignments[organism][ref] += 1
            else:
                organism_assignments["UNMAPPED"]["NA"] += 1
    return organism_assignments


def run():
    df = pd.read_table(
        snakemake.input.metadata,
        header=None,
        names=["organism", "contamination", "accession"],
        index_col="accession",
    )
    organism_assignments = get_organism_assignments(snakemake.input.bam, df)

    with open(snakemake.output[0], "w") as f_out:
        for category, d in organism_assignments.items():
            for accession, count in d.items():
                if accession == "NA":
                    print("{count}\tUnmapped".format(count=count), file=f_out)
                else:
                    print("{count}\t{category}\t{accession}".format(
                        count=count, category=category, accession=accession), file=f_out)


# ========================================================================
# ========================================================================
run()
