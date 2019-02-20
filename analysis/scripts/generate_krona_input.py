from collections import Counter
import json
import pysam


def get_ena_accession(accession_string):
    """Extract the accession ID from a header string"""
    return accession_string.split("|")[1]


def get_organism_assignments(bamfile: str):
    """Count the assignments for each read in a bam file.
    :param bamfile: Path to bam file.
    :return: A dictionary with the keys being organisms and the values being
    Counters with keys being accession IDs and values being counts.
    """
    organism_assignments = Counter()
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for record in bam:
            if record.reference_name is not None:
                organism_assignments[record.reference_name] += 1
            else:
                organism_assignments["UNMAPPED"] += 1
    return organism_assignments


def get_krona_entry(accession: str, count: int, lookup_table):
    """Generate a row entry for a text file required to plot in krona."""
    if accession == "UNMAPPED":
        return "{count}\tUnmapped".format(count=count)
    else:
        if accession.startswith("ENA"):
            accession = get_ena_accession(accession)
        taxonomy = "\t".join(lookup_table[accession])
        return "{count}\t{taxonomy}".format(count=count, taxonomy=taxonomy)

def load_json(fname):
    with open(fname, 'r') as f_in:
        data = json.loads(f_in.read(), fp)
    return data[0]

def run():
    organism_assignments = get_organism_assignments(snakemake.input.bam)
    lookup_table = load_json(snakemake.input.taxonomy)

    with open(snakemake.output[0], "w") as f_out:
        for accession, count in organism_assignments.items():
            krona_entry = get_krona_entry(accession, count, lookup_table)
            print(krona_entry, f_out)


# ========================================================================
# ========================================================================
run()
