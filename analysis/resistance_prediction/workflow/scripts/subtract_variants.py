import sys
from typing import NamedTuple, Optional, List

sys.stderr = open(snakemake.log[0], "w")
from cyvcf2 import VCF, Writer, Variant


class Genotype(NamedTuple):
    allele1: int
    allele2: int

    def is_null(self) -> bool:
        """Is the genotype null. i.e. ./."""
        return self.allele1 == -1 and self.allele2 == -1

    def is_hom(self) -> bool:
        """Is the genotype homozygous"""
        if self.is_null():
            return False
        if self.allele1 == -1 or self.allele2 == -1:
            return True
        return self.allele1 == self.allele2

    def is_het(self) -> bool:
        """Is the genotype heterozyhous"""
        return not self.is_null() and not self.is_hom()

    def is_hom_ref(self) -> bool:
        """Is genotype homozygous reference?"""
        return self.is_hom() and (self.allele1 == 0 or self.allele2 == 0)

    def is_hom_alt(self) -> bool:
        """Is genotype homozygous alternate?"""
        return self.is_hom() and (self.allele1 > 0 or self.allele2 > 0)

    def alt_index(self) -> Optional[int]:
        """If the genotype is homozygous alternate, returns the 0-based index of the
        alt allele in the alternate allele array.
        """
        if not self.is_hom_alt():
            return None
        return max(self.allele1, self.allele2) - 1

    def allele_index(self) -> Optional[int]:
        """The index of the called allele"""
        if self.is_hom_ref() or self.is_null():
            return 0
        elif self.is_hom_alt():
            return self.alt_index() + 1
        else:
            raise NotImplementedError(f"Het Genotype is unexpected: {self}")

    @staticmethod
    def from_arr(arr: List[int]) -> "Genotype":
        alleles = [a for a in arr if type(a) is int]
        if len(alleles) < 2:
            alleles.append(-1)
        return Genotype(*alleles)


def region_for_record(record: Variant) -> str:
    # region strings are 1-based inclusive [start, end], sp [1, 1] return pos 1 only
    return f"{record.CHROM}:{record.POS}-{record.POS+len(record.REF)-1}"


def overlaps_match(q: str, qpos: int, p: str, ppos: int) -> bool:
    qend = qpos + len(q)
    pend = ppos + len(p)
    qidx = slice(max(0, ppos - qpos), pend - qpos)
    pidx = slice(max(0, qpos - ppos), qend - ppos)
    return q[qidx] == p[pidx]


panel_vcf = VCF(snakemake.input.panel)
query_vcf = VCF(snakemake.input.query)
output_vcf = Writer(snakemake.output.vcf, tmpl=query_vcf)

for qrecord in query_vcf:
    gt = Genotype.from_arr(qrecord.genotypes[0])
    alt_idx = gt.alt_index()
    if alt_idx is None:
        continue
    qseq = qrecord.ALT[alt_idx]
    region = region_for_record(qrecord)
    has_match_in_panel = False

    for precord in panel_vcf(region):
        for i, pseq in enumerate(precord.ALT):
            if overlaps_match(qseq, qrecord.POS, pseq, precord.POS):
                print(
                    f"Query record at {qrecord.CHROM}:{qrecord.POS} has an overlapping "
                    f"match with panel record {precord.CHROM}:{precord.POS} ALT number "
                    f"{i}",
                    file=sys.stderr
                )
                has_match_in_panel = True
                break
        if has_match_in_panel:
            break

    if has_match_in_panel:
        continue
    else:
        output_vcf.write_record(qrecord)

output_vcf.close()


"""Tests for overlaps_match
q = "ATC"
qpos = 10
p = "TCA"
ppos = 11
assert overlaps_match(q, qpos, p, ppos)
q = "ATC"
qpos = 11
p = "CAT"
ppos = 10
assert overlaps_match(q, qpos, p, ppos)
q = "ATC"
qpos = 10
p = "ATC"
ppos = 10
assert overlaps_match(q, qpos, p, ppos)
q = "AT"
qpos = 10
p = "ATC"
ppos = 10
assert overlaps_match(q, qpos, p, ppos)
q = "ATC"
qpos = 10
p = "AT"
ppos = 10
assert overlaps_match(q, qpos, p, ppos)
q = "ATCGA"
qpos = 10
p = "TCG"
ppos = 11
assert overlaps_match(q, qpos, p, ppos)
q = "TCG"
qpos = 11
p = "ATCGA"
ppos = 10
assert overlaps_match(q, qpos, p, ppos)
"""
