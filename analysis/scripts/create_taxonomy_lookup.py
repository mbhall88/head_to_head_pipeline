import json
from urllib.error import HTTPError

from Bio import Entrez, SeqIO

Entrez.email = snakemake.params.email

Entrez.api_key = snakemake.params.api_key


def get_accessions(fname: str):
    """Get the accession IDs from the given TSV file"""
    accessions = dict()

    with open(fname, 'r') as tsv:
        for line in tsv:
            category, contam, accession = line.strip().split()
            if accession.startswith("ENA"):
                accession = get_ena_accession(accession)
            if accession not in accessions:
                accessions[accession] = [category]

    return accessions


def get_taxonomy(accession: str):
    """Fetch taxonomy information for a given accession ID.
    :param accession: Accession ID
    :return: A list representing the taxonomy for the accession ID
    https://stackoverflow.com/a/28355078/5299417"""
    handle = Entrez.efetch(db="nucleotide", id=accession,
                           rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    taxonomy = record.annotations["taxonomy"]
    taxonomy.extend([record.annotations["organism"]])
    return taxonomy


def get_ena_accession(accession_string):
    """Extract the accession ID from a header string"""
    return accession_string.split("|")[1]


def run():
    accessions = get_accessions(snakemake.input.metadata)
    for acc, category in accessions.items():
        try:
            taxonomy = get_taxonomy(acc)
        except HTTPError as err:
            if err.code == 400:
                taxonomy = category
            else:
                raise(err)
        accessions[acc] = taxonomy

    with open(snakemake.output.taxonomy, 'w') as f_out:
        json.dump(accessions, f_out)


# ========================================================================
# ========================================================================
run()
