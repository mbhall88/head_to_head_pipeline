"""Snakemake script to extract genotype matrix from a VCF, paying attention to the
filter (FT) tag for each sample. If there is a FT tag present, the genotype is encoded
as -2"""
import sys

sys.stderr = open(snakemake.log[0], "w")

import numpy as np
from cyvcf2 import VCF


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


in_vcf = snakemake.input.vcf
outstream = open(snakemake.output.matrix, "w")
delim = snakemake.params.delim
FILTERED = -2

vcf_reader = VCF(in_vcf)
samples = vcf_reader.samples
N = len(samples)

# print header
print(delim.join(["CHROM", "POS", *samples]), file=outstream)

for variant in vcf_reader:
    chrom = variant.CHROM
    pos = variant.POS
    genotypes = np.full(N, FILTERED)
    if variant.FILTER is None:  # i.e. PASS or .
        mask = variant.format("FT") == "PASS"
        # replace gt for samples that pass filter with their correct gt
        genotypes[mask] = variant.genotype.array()[:, 0][mask]

    row = delim.join([chrom, pos, *map(str, genotypes)])
    print(row, file=outstream)


vcf_reader.close()
outstream.close()


eprint("Done!")
