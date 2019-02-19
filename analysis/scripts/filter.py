import pandas as pd
import pysam
import gzip
from contextlib import ExitStack

# =======================================================
# Functions
# =======================================================

def to_fastq_string(record):
    """Extracts the information required for a fastq entry from a BAM entry
    :param record: BAM record
    :return: A string representation of the four fastq lines
    Tested here: https://github.com/mbhall88/bam2fastq/blob/master/test_bam2fastq.py
    """
    qual = record.get_forward_qualities()
    seq = record.get_forward_sequence()

    if seq is None or qual is None:
        return ""
    else:
        qual = "".join([chr(q + 33) for q in qual])

    header = "@{name}".format(name=record.query_name)

    assert len(seq) == len(qual), "Sequence is not the same length as quality string"

    return "{header}\n{seq}\n+\n{qual}".format(header=header, seq=seq, qual=qual)


def is_wanted(record, df):
    """Determines whether a record is to be kept or not.
    :param record: BAM record.
    :param df: Pandas dataframe holding information about what to filter on.
    Should have at least a contamination column with 0 or 1 representing whether
    the organism for that row is contamination or not.
    :return: A boolean indicating whether the read should be kept.
    """
    if record.reference_name is None or not is_primary(record):
        return False
    return not bool(df.loc[record.reference_name]["contamination"])


def is_primary(record):
    """Tests whether a record is a primary alignment."""
    return not any([record.is_unmapped, record.is_secondary, record.is_supplementary])


# =======================================================

# =======================================================
# Main
# =======================================================
df = pd.read_table(
    snakemake.input.metadata,
    header=None,
    names=["organism", "contamination", "accession"],
    index_col="accession",
)

with ExitStack() as stack:
    bam = stack.enter_context(pysam.AlignmentFile(snakemake.input.bam, "rb"))
    fastq_fh = stack.enter_context(gzip.open(snakemake.output.fastq, "wb"))

    for record in bam:
        keep_read = is_wanted(record, df)

        if keep_read:
            fastq_fh.write("{}\n".format(to_fastq_string(record)).encode("utf-8"))
