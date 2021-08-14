import sys
from collections import defaultdict
from typing import NamedTuple, Optional, List

sys.stderr = open(snakemake.log[0], "w")
from cyvcf2 import VCF, Writer, Variant

TP = "TP"
FP = "FP"
FN = "FN"


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


TAG = "CLF"
truth_rdr = VCF(snakemake.input.truth)
query_rdr = VCF(snakemake.input.query)
query_rdr.add_info_to_header(
    {
        "ID": TAG,
        "Description": "Classification of record",
        "Type": "String",
        "Number": ".",
    }
)
query_wtr = Writer(snakemake.output.annotated_query_vcf, tmpl=query_rdr)

classifications = []
classified_qrecords = set()
classified_trecords = set()

for query_record in query_rdr:
    if query_record.FILTER is not None:
        continue
    record_clfs = []
    query_gt = Genotype.from_arr(query_record.genotypes[0])
    qalt_idx = query_gt.alt_index()
    if qalt_idx is None:
        continue
    qseq = query_record.ALT[qalt_idx]
    region = region_for_record(query_record)

    overlapping_regions = truth_rdr(region)
    i = 0
    for truth_record in overlapping_regions:
        i += 1
        if truth_record in classified_trecords:
            continue
        truth_gt = Genotype.from_arr(truth_record.genotypes[0])
        talt_idx = truth_gt.alt_index()
        if talt_idx is None:
            continue
        tseq = truth_record.ALT[talt_idx]
        if overlaps_match(qseq, query_record.POS, tseq, truth_record.POS):
            clf = TP
            print(
                f"Query record at {query_record.CHROM}:{query_record.POS} has an overlapping "
                f"match with truth record {truth_record.CHROM}:{truth_record.POS} ALT number "
                f"{i}: {clf}",
                file=sys.stderr,
            )
        else:
            clf = FP
            print(
                f"Query record at {query_record.CHROM}:{query_record.POS} has no overlapping "
                f"match with truth record {truth_record.CHROM}:{truth_record.POS} ALT number "
                f"{i}: {clf}",
                file=sys.stderr,
            )

        classifications.append(clf)
        record_clfs.append(clf)
        classified_trecords.add((truth_record.CHROM, truth_record.POS))

    if i == 0:
        classifications.append(FP)
        record_clfs.append(FP)
        print(
            f"Query record at {query_record.CHROM}:{query_record.POS} has no overlapping "
            f"truth record: {FP}",
            file=sys.stderr,
        )

    classified_qrecords.add(query_record.ID)
    query_record.INFO[TAG] = ",".join(record_clfs)
    query_wtr.write_record(query_record)

query_wtr.close()

truth_rdr = VCF(snakemake.input.truth)
query_rdr = VCF(snakemake.input.query)
truth_rdr.add_info_to_header(
    {
        "ID": TAG,
        "Description": "Classification of record",
        "Type": "String",
        "Number": ".",
    }
)
truth_wtr = Writer(snakemake.output.annotated_truth_vcf, tmpl=truth_rdr)

for truth_record in truth_rdr:
    if truth_record.FILTER is not None:
        continue
    record_clfs = []
    truth_gt = Genotype.from_arr(truth_record.genotypes[0])
    talt_idx = truth_gt.alt_index()
    if talt_idx is None:
        continue
    tseq = truth_record.ALT[talt_idx]
    region = region_for_record(truth_record)

    overlapping_regions = query_rdr(region)
    i = 0
    for query_record in overlapping_regions:
        if query_record.FILTER is not None:
            continue
        i += 1

        query_gt = Genotype.from_arr(query_record.genotypes[0])
        qalt_idx = query_gt.alt_index()
        if qalt_idx is None:
            continue
        qseq = query_record.ALT[qalt_idx]
        if overlaps_match(qseq, query_record.POS, tseq, truth_record.POS):
            clf = TP
            print(
                f"Query record at {query_record.CHROM}:{query_record.POS} has an overlapping "
                f"match with truth record {truth_record.CHROM}:{truth_record.POS} ALT number "
                f"{i}: {clf}",
                file=sys.stderr,
            )
        else:
            clf = FN
            print(
                f"Query record at {query_record.CHROM}:{query_record.POS} has no overlapping "
                f"match with truth record {truth_record.CHROM}:{truth_record.POS} ALT number "
                f"{i}: {clf}",
                file=sys.stderr,
            )

        record_clfs.append(clf)

    if i == 0:
        classifications.append(FN)
        record_clfs.append(FN)
        print(
            f"Truth record at {truth_record.CHROM}:{truth_record.POS} has no overlapping "
            f"query record: {FN}",
            file=sys.stderr,
        )

    truth_record.INFO[TAG] = ",".join(record_clfs)
    truth_wtr.write_record(truth_record)


truth_wtr.close()
