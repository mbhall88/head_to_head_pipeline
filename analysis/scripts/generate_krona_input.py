import random
import string
import time
from collections import Counter
from urllib.error import HTTPError

import pysam
from Bio import Entrez, SeqIO


def get_taxonomy(accession: str):
    """Fetch taxonomy information for a given accession ID.
    :param accession: Accession ID
    :return: A list representing the taxonomy for the accession ID

    https://stackoverflow.com/a/28355078/5299417"""
    random_user = ''.join(random.choice(string.ascii_lowercase)
                          for _ in range(10))
    Entrez.email = "{}@gmail.com".format(random_user)
    while True:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gb", retmode="text")
            break
        except HTTPError:
            time.sleep(10)
            continue

    record = SeqIO.read(handle, "genbank")
    taxonomy = record.annotations["taxonomy"]
    taxonomy.extend([record.annotations["organism"]])
    return taxonomy


def get_ena_accession(accession_string):
    """Extract the accession ID from a header string"""
    return accession_string.split("|")[-1]


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


def get_krona_entry(accession: str, count: int):
    """Generate a row entry for a text file required to plot in krona."""
    if accession == "UNMAPPED":
        return "{count}\tUnmapped".format(count=count)
    else:
        if accession.startswith("ENA"):
            accession = get_ena_accession(accession)
        taxonomy = "\t".join(get_taxonomy(accession))
        return "{count}\t{taxonomy}".format(count=count, taxonomy=taxonomy)


def test_krona_entry():
    assert get_krona_entry("UNMAPPED", 3) == "3\tUnmapped"
    assert get_krona_entry("ENA|BBHB01000083|BBHB01000083.1",
                           5) == "5\tBacteria\tActinobacteria\tCorynebacteriales\tMycobacteriaceae\tMycobacterium\tMycobacterium marinum str. Europe"


def run():
    organism_assignments = get_organism_assignments(snakemake.input.bam)

    with open(snakemake.output[0], "w") as f_out:
        for accession, count in organism_assignments.items():
            krona_entry = get_krona_entry(accession, count)
            print(krona_entry, f_out)


# ========================================================================
# ========================================================================
test_krona_entry()
run()
