import json
import sys
from typing import NamedTuple, Optional, List

sys.stderr = open(snakemake.log[0], "w")


class Genotype(NamedTuple):
    allele1: int
    allele2: int

    def __str__(self) -> str:
        if self.is_null():
            return "NULL"
        elif self.is_het():
            return "HET"
        elif self.is_hom_ref():
            return "REF"
        elif self.is_hom_alt():
            return "ALT"
        else:
            raise ValueError(
                f"Unknown genotype category: ({self.allele1}, {self.allele2})"
            )

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


def extract_call_for_mutation(d: dict) -> str:
    filters = d["info"]["filter"]
    if filters:
        return "FILT"
    else:
        return str(Genotype.from_arr(d["genotype"]))


ont_json = snakemake.input.nanopore_predictions
ill_json = snakemake.input.illumina_predictions
delim = snakemake.params.get("delim", ",")

outstream = open(snakemake.output.table, "w")

sample = snakemake.wildcards.sample
with open(ill_json) as f:
    ill_data = json.load(f)[sample]["variant_calls"]

with open(ont_json) as f:
    ont_data = json.load(f)[sample]["variant_calls"]

print(f"mutation{delim}illumina_call{delim}nanopore_call", file=outstream)
for var, ill_d in ill_data.items():
    ont_d = ont_data[var]
    print(
        delim.join(
            (var, extract_call_for_mutation(ill_d), extract_call_for_mutation(ont_d))
        ),
        file=outstream,
    )

outstream.close()
