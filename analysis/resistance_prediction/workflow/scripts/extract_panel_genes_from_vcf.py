import sys

sys.stderr = open(snakemake.log[0], "w")
from typing import TextIO, Set, Dict

from loguru import logger
from intervaltree import IntervalTree
from cyvcf2 import VCF, Writer


def extract_genes_from_panel(stream: TextIO) -> Set[str]:
    genes = set()
    for line in map(str.rstrip, stream):
        if not line:
            continue
        fields = line.split("\t")
        if gene := fields[0]:
            genes.add(gene)
    return genes


def extract_intervals_for_genes_from_gff(
    genes: Set[str], gff_stream: TextIO, padding: int = 0
) -> IntervalTree:
    intervals = []
    for row in map(str.rstrip, gff_stream):
        if row.startswith("#") or not row:
            continue
        fields = row.split("\t")
        if fields[2].lower() != "gene":
            continue

        attributes = attributes_dict_from_str(fields[8])
        name = attributes.get("gene", attributes.get("Name", None))
        if name is None:
            logger.warning(f"No gene/Name attribute for ID {attributes['ID']}")
            continue

        if name not in genes:
            continue

        start = (int(fields[3]) - 1) - padding  # GFF start is 1-based inclusive
        end = (int(fields[4]) - 1) + padding  # GFF end is 1-based exclusive
        intervals.append((start, end, name))

    return IntervalTree.from_tuples(intervals)


def attributes_dict_from_str(s: str) -> Dict[str, str]:
    d = dict()
    for pairs in s.split(";"):
        k, v = pairs.split("=")
        if k in d:
            raise KeyError(f"Attribute key {k} appears twice")
        d[k] = v
    return d


##########################################################
# MAIN
##########################################################
padding: int = snakemake.params.padding

logger.info("Extracting gene names from panel...")
with open(snakemake.input.panel) as istream:
    genes = extract_genes_from_panel(istream)

logger.success(f"Extracted {len(genes)} genes from the panel")

logger.info("Extracting intervals for genes from GFF...")
with open(snakemake.input.annotation) as istream:
    ivtree = extract_intervals_for_genes_from_gff(genes, istream, padding)
logger.success(f"Intervals extracted for {len(ivtree)} genes")

logger.info(
    "Extracting those VCF records that fall within the gene intervals and altering "
    "their CHROM and POS accordingly..."
)
vcf_reader = VCF(snakemake.input.vcf)

logger.debug("Adding genes to header...")
for iv in ivtree:
    vcf_reader.add_to_header(f"##contig=<ID={iv.data},length={iv.end-iv.begin}>")
logger.debug("Genes added to header")

vcf_writer = Writer(snakemake.output.vcf, tmpl=vcf_reader)

for record in vcf_reader:
    iv = ivtree[record.start]
    if len(iv) > 1:
        raise IndexError(
            f"VCF record at POS {record.POS} overlaps with more than 1 gene: {iv}"
        )
    if len(iv) == 0:
        continue
    iv = list(iv)[0]
    chrom = iv.data
    norm_pos = record.start - iv.begin
    record.set_pos(norm_pos)
    record.CHROM = chrom
    vcf_writer.write_record(record)

vcf_writer.close()
vcf_reader.close()

logger.success("Done!")
