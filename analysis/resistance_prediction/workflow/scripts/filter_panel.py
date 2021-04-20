import sys

sys.stderr = open(snakemake.log[0], "w")

import re
from typing import Set
from dataclasses import dataclass

frameshift_lengths: Set[int] = set(snakemake.params.frameshift_lengths)
exclude_frameshifts_in: Set[str] = set(snakemake.params.exclude)


@dataclass
class Variant:
    before: str
    pos: int
    after: str

    @staticmethod
    def from_str(var: str) -> "Variant":
        before, pos, after = re.split(r"(-?\d+)", var, maxsplit=1)
        return Variant(before=before, after=after, pos=int(pos))

    def is_frameshift(self) -> bool:
        len_diff = len(self.before) - len(self.after)
        return abs(len_diff) in frameshift_lengths


num_filtered = 0

with open(snakemake.input.panel) as istream, open(
    snakemake.output.panel, "w"
) as ostream:
    for row in map(str.rstrip, istream):
        gene, var, alpha = row.split("\t")[:3]
        is_dna = alpha == "DNA"
        if (
            gene in exclude_frameshifts_in
            and is_dna
            and Variant.from_str(var).is_frameshift()
        ):
            num_filtered += 1
            continue
        print(row, file=ostream)

print(f"Filtered out {num_filtered} rows", file=sys.stderr)
