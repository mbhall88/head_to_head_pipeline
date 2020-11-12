import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
from typing import TextIO, Union
from intervaltree import IntervalTree, Interval
import logging
from contextlib import ExitStack

PathLike = Union[Path, str]

log_level = logging.DEBUG if snakemake.params.verbose else logging.INFO
logging.basicConfig(format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level)


def load_loci_info(instream: TextIO) -> IntervalTree:
    intervals = []
    _ = next(instream)  # skip header
    for row in map(str.rstrip, instream):
        fields = row.split(",", maxsplit=6)
        start, end, name = fields[2:5]

        intervals.append((int(start), int(end), name))

    return IntervalTree.from_tuples(intervals)


def load_mask(instream: TextIO) -> IntervalTree:
    intervals = []
    for row in map(str.rstrip, instream):
        fields = row.split("\t")
        start = int(fields[1])  # BED is zero-based
        end = int(fields[2])  # BED end is exclusive
        chrom = fields[0]

        intervals.append((start, end, chrom))
    return IntervalTree.from_tuples(intervals)


with open(snakemake.input.loci_info) as instream:
    ivtree = load_loci_info(instream)

logging.info(f"Loaded {len(ivtree)} loci")

with open(snakemake.input.mask) as instream:
    mask = load_mask(instream)

logging.info(f"Loaded {len(mask)} mask regions")

masked_tree = IntervalTree(ivtree)

for iv in mask:
    masked_tree.remove_overlap(iv.begin, iv.end)

full_len = 0
for iv in ivtree:
    full_len += iv.length()

masked_len = 0
for iv in masked_tree:
    masked_len += iv.length()

logging.info(
    f"{len(ivtree)-len(masked_tree)} ({1-(len(masked_tree)/len(ivtree)):.2%}) loci removed"
)
logging.info(
    f"{full_len-masked_len}bp ({1-(masked_len/full_len):.2%}) of the genome removed"
)

logging.info("Creating masked loci info file...")

with ExitStack() as stack:
    outstream = stack.enter_context(open(snakemake.output.masked_loci_info, "w"))
    instream = stack.enter_context(open(snakemake.input.loci_info))
    outstream.write(next(instream))  # skip header
    for row in map(str.rstrip, instream):
        fields = row.split(",", maxsplit=6)
        start, end, name = fields[2:5]
        iv = Interval(int(start), int(end), data=name)
        if iv in masked_tree:
            print(row, file=outstream)  # we stripped the newline
        else:
            logging.debug(f"Discarding masked locus {name}")

logging.info("Done!")
