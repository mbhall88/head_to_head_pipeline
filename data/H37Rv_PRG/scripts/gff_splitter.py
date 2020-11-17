#!/usr/bin/env python3
import logging
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from itertools import repeat
from pathlib import Path
from typing import TextIO, Tuple, Dict, List

import click
import numpy as np
from intervaltree import IntervalTree, Interval

Contig = str
Seq = str
Index = Dict[Contig, Seq]


class Strand(Enum):
    Forward = "+"
    Reverse = "-"
    NotRelevant = "."
    Unknown = "?"

    def __str__(self) -> str:
        return str(self.value)


@dataclass
class GffFeature:
    seqid: Contig
    source: str
    method: str  # correct term is type, but that is a python reserved variable name
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    score: float
    strand: Strand
    phase: int
    attributes: Dict[str, str]

    @staticmethod
    def from_str(s: str) -> "GffFeature":
        fields = s.split("\t")
        score = 0 if fields[5] == "." else float(fields[5])
        phase = -1 if fields[7] == "." else int(fields[7])
        attr_fields = fields[-1].split(";")
        attributes = {k: v for k, v in map(str.split, attr_fields, repeat("="))}
        return GffFeature(
            seqid=fields[0],
            source=fields[1],
            method=fields[2],
            start=int(fields[3]),
            end=int(fields[4]),
            score=score,
            strand=Strand(fields[6]),
            phase=phase,
            attributes=attributes,
        )

    def slice(self, zero_based: bool = True) -> Tuple[int, int]:
        """Get a tuple for slicing a python object.
        The reason this method is required is that GFF uses 1-based INCLUSIVE
        coordinates. Meaning the end position is also included in the slice.
        """
        if zero_based:
            return self.start - 1, self.end
        return self.start, self.end + 1


class DuplicateContigsError(Exception):
    pass


class OverlapError(Exception):
    pass


class NoDataError(Exception):
    pass


def is_header(s: str) -> bool:
    if not s:
        return False
    return s[0] == ">"


def has_adjoining_neighbour_to_left(iv: Interval, tree: IntervalTree) -> bool:
    """Check if there is *no* gap between the given interval and the neighbour to the
    left.
    """
    return tree.overlaps(iv.begin - 1)


def has_adjoining_neighbour_to_right(iv: Interval, tree: IntervalTree) -> bool:
    """Check if there is *no* gap between the given interval and the neighbour to the
    right.
    Note: we don't add 1 to end as it is normally exclusive, but is inclusive in overlap
    queries.
    """
    return tree.overlaps(iv.end)


# test
# tree = IntervalTree.from_tuples([(1, 5), (10, 20)])
# iv = Interval(6, 8)
# assert not has_adjoining_neighbour_to_left(iv, tree)
# assert not has_adjoining_neighbour_to_right(iv, tree)

# iv = Interval(5, 9)
# assert has_adjoining_neighbour_to_left(iv, tree)
# assert not has_adjoining_neighbour_to_right(iv, tree)

# iv = Interval(7, 10)
# assert not has_adjoining_neighbour_to_left(iv, tree)
# assert has_adjoining_neighbour_to_right(iv, tree)


def load_mask(instream: TextIO) -> Dict[Contig, IntervalTree]:
    """Loads a bedfile mask into an interval tree"""
    intervals = defaultdict(IntervalTree)
    for row in map(str.rstrip, instream):
        fields = row.split("\t")
        start = int(fields[1])  # BED is zero-based
        end = int(fields[2])  # BED end is exclusive
        chrom = fields[0]

        intervals[chrom].add(Interval(start, end))
    return intervals


def index_fasta(stream: TextIO) -> Index:
    fasta_index: Index = dict()
    sequence: List[Seq] = []
    name: Contig = ""
    for line in map(str.rstrip, stream):
        if not line:
            continue
        if is_header(line):
            if sequence and name:
                fasta_index[name] = "".join(sequence)
                sequence = []
            name = line.split()[0][1:]
            if name in fasta_index:
                raise DuplicateContigsError(
                    f"Contig {name} occurs multiple times in the fasta file."
                )
            continue
        else:
            sequence.append(line)
    if name and sequence:
        fasta_index[name] = "".join(sequence)

    return fasta_index


def data_reducer(current_reduced_data: str, new_data: str) -> str:
    """This function is used when merging overlaps in the features index tree. By
    default, when merging overlaps, the data is removed. However, in this script we
    need the interval data as it hold the name of the interval. Therefore, this function
    tells intervalltree how to merge the data field in intervals.
    """
    return f"{current_reduced_data}+{new_data}"


def slice_seq(seq: Seq, interval: Interval) -> Seq:
    i = interval.begin
    j = interval.end
    return seq[i:j]


def construct_feature_trees(
    gff: TextIO, types: Tuple[str]
) -> Dict[Contig, IntervalTree]:
    feature_trees: Dict[Contig, IntervalTree] = defaultdict(IntervalTree)

    for line in map(str.rstrip, gff):
        if not line or line.startswith("#"):
            continue

        feature = GffFeature.from_str(line)
        if feature.method not in types:
            continue

        start, end = feature.slice(zero_based=True)

        if "Name" in feature.attributes:
            name = feature.attributes["Name"]
        elif "ID" in feature.attributes:
            name = feature.attributes["ID"]
        else:
            name = f"{feature.method};{start}-{end}"
            logging.warning(
                f"Can't find a Name or ID for feature {feature}. Using {name}"
            )

        iv = Interval(start, end, data=name)
        feature_trees[feature.seqid].add(iv)
    return feature_trees


def merge_short_intervals(tree: IntervalTree, min_len: int = 1) -> IntervalTree:
    """Merge short intervals with their neighbour."""
    merged_tree = IntervalTree()
    for i, iv in enumerate(sorted(tree)):
        if iv.length() >= min_len:
            merged_tree.add(iv)
            continue
        if has_adjoining_neighbour_to_right(iv, tree):
            # have to add 2 to end to force merge as is no inclusive
            merged_tree.add(Interval(iv.begin, iv.end + 2, data=iv.data))
        elif has_adjoining_neighbour_to_left(iv, tree):
            merged_tree.add(Interval(iv.begin - 1, iv.end, data=iv.data))

    merged_tree.merge_overlaps(data_reducer=data_reducer)

    if any(iv.length() < min_len for iv in merged_tree):
        return merge_short_intervals(merged_tree, min_len=min_len)
    else:
        return merged_tree


# test
# min_len = 3
# tree = IntervalTree.from_tuples([(1, 5), (6, 8), (10, 20)])
# actual = merge_short_intervals(tree, min_len)
# expected = IntervalTree.from_tuples([(1, 5), (10, 20)])
# assert actual == expected, actual

# tree = IntervalTree.from_tuples([(1, 5), (8, 10), (10, 20)])
# actual = merge_short_intervals(tree, min_len)
# expected = tree = IntervalTree.from_tuples([(1, 5), (8, 20, "None+None")])
# assert actual == expected, actual

# min_len = 50
# tree = IntervalTree.from_tuples([(1, 5), (7, 10), (10, 20)])
# actual = merge_short_intervals(tree, min_len)
# expected = IntervalTree()
# assert actual == expected, actual

# min_len = 5
# tree = IntervalTree.from_tuples([(1, 5), (6, 10), (10, 20)])
# actual = merge_short_intervals(tree, min_len)
# expected = IntervalTree.from_tuples([(6, 20, "None+None")])
# assert actual == expected, actual

# min_len = 5
# tree = IntervalTree.from_tuples([(1, 5), (5, 10), (10, 20)])
# actual = merge_short_intervals(tree, min_len)
# expected = IntervalTree.from_tuples([(1, 10, "None+None"), (10, 20)])
# assert actual == expected, actual

# min_len = 11
# tree = IntervalTree.from_tuples([(1, 5), (5, 10), (10, 20)])
# actual = merge_short_intervals(tree, min_len)
# expected = IntervalTree.from_tuples([(1, 20, "None+None+None")])
# assert actual == expected, actual


def infer_interval_type(interval: Interval) -> str:
    if not interval.data:
        raise NoDataError(f"Expected data in interval, but gone none: {interval}")
    names = str(interval.data).split("+")
    num_intervals = len(names)
    if num_intervals > 1:
        if all("IGR" in name for name in names):
            return "merged_igrs"
        elif any("IGR" in name for name in names):
            return "merged_feature_and_igr"
        else:
            return "merged_features"
    return "igr" if "IGR" in names[0] else "feature"


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-f",
    "--fasta",
    help="FASTA file to split.",
    type=click.File(mode="r"),
    default="-",
    show_default=True,
    required=True,
)
@click.option(
    "-g",
    "--gff",
    help="GFF3 file to base split coordinates on.",
    type=click.File(mode="r"),
    required=True,
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(file_okay=False, resolve_path=True, writable=True),
    default=".",
    show_default=True,
    help="The directory to write the output files to.",
)
@click.option(
    "-m",
    "--mask",
    type=click.File(mode="r"),
    help="BED file of positions to mask from the GFF.",
)
@click.option(
    "-F",
    "--min-overlap",
    help=(
        "Minimum overlap (fraction) between the mask and a locus/feature that causes "
        "the locus to be removed. Setting to 0 means even a 1bp overlap causes the "
        "locus to be removed"
    ),
    default=0.3,
    show_default=True,
)
@click.option(
    "--types",
    help=(
        "The feature types to split on. Separate types by a space or pass option "
        "mutiple times."
    ),
    multiple=True,
    default=["gene"],
    show_default=True,
)
@click.option(
    "--min-len",
    help=(
        "The minimum length of the chunks to output. If a chunk is shorter than this "
        "value, it is joined to a neighbouring chunk."
    ),
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "--max-len",
    help=(
        "The maximum length of the chunks to output. Use this option with caution. If "
        "a chunk is greater than this value it is discarded."
    ),
    metavar="INTEGER",
    type=float,
    default=float("inf"),
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    fasta: TextIO,
    gff: TextIO,
    mask: TextIO,
    outdir: str,
    min_overlap: float,
    types: Tuple[str],
    min_len: int,
    max_len: float,
    verbose: bool,
):
    """Splits a FASTA file into chunks based on a GFF3 file.
    The splits produced are based on the --types given and everything inbetween. For
    example, the default --types is 'gene'. In this case, the coordinates for each gene
    are cut out of the FASTA file, as well as the bits inbetween - intergenic regions
    (IGRs).
    In the output, start is 0-based inclusive and end is 0-based non-inclusive.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    logging.info("Indexing fasta file...")
    index: Index = index_fasta(fasta)
    logging.info(f"{len(index)} contig(s) indexed in the input file.")

    index_trees: Dict[Contig, IntervalTree] = {
        contig: IntervalTree([Interval(1, len(seq))]) for contig, seq in index.items()
    }

    masked_intervals: Dict[Contig, np.ndarray] = dict()
    if mask is not None:
        logging.info("Loading the mask...")
        ivtrees = load_mask(mask)
        for contig in index:
            gsize = len(index[contig])
            arr = np.zeros(gsize, dtype=np.bool)
            if contig in ivtrees:
                for i in range(len(arr)):
                    if ivtrees[contig][i]:
                        arr[i] = True
            masked_intervals[contig] = arr

    logging.info("Constructing interval tree for features...")
    feature_trees: Dict[Contig, IntervalTree] = construct_feature_trees(gff, types)

    for contig in index_trees:
        logging.info(f"Inferring intergenic region interval(s) for {contig}...")
        for iv in feature_trees[contig]:
            index_trees[contig].chop(iv.begin, iv.end)

        intervals_with_names = set()
        for iv in index_trees[contig]:
            name = f"IGR:{iv.begin}-{iv.end}"
            intervals_with_names.add(Interval(iv.begin, iv.end, data=name))

        index_trees[contig] = IntervalTree(intervals_with_names)

        logging.info(
            f"Found {len(index_trees[contig])} intergenic region interval(s) for "
            f"{contig}"
        )

    logging.debug("Joining IGR and feature trees...")
    trees = {
        contig: index_trees[contig].union(feature_trees[contig])
        for contig in feature_trees
    }

    # go through each interval, and if it overlaps by min_overlap or more, remove it
    ivs_to_remove: Dict[Contig, List[Interval]] = defaultdict(list)
    for contig, ivtree in trees.items():
        if contig not in masked_intervals:
            continue
        logging.info(f"Masking {contig}...")
        for iv in ivtree:
            frac_masked = masked_intervals[contig][iv.begin : iv.end].mean()
            if frac_masked > min_overlap:
                logging.debug(f"Removing {iv.data}")
                ivs_to_remove[contig].append(iv)

    for contig, ivs in ivs_to_remove.items():
        logging.info(f"{len(trees[contig])} feature(s) before masking {contig}...")
        for iv in ivs:
            trees[contig].remove(iv)
        logging.info(f"{len(trees[contig])} feature(s) after masking {contig}")

    for contig in feature_trees:
        logging.info(f"Merging overlapping features for {contig}...")
        feature_trees[contig].merge_overlaps(data_reducer=data_reducer)

    for contig in trees:
        logging.info(f"Merging overlapping features for {contig}...")
        logging.info(f"{len(trees[contig])} features before merging...")
        trees[contig].merge_overlaps(data_reducer=data_reducer)
        logging.info(f"{len(trees[contig])} features after merging")

    for contig, tree in trees.items():
        logging.info(f"Merging short intervals for {contig}...")
        trees[contig] = merge_short_intervals(tree, min_len=min_len)
        logging.info(f"{len(trees[contig])} interval(s) after merging short ones.")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    file_mapping_path = outdir / "loci-info.csv"
    mapping_stream = file_mapping_path.open("w")
    print(
        ",".join(["filename", "type", "start", "end", "name", "contig"]),
        file=mapping_stream,
    )

    for contig, tree in trees.items():
        logging.info(f"Writing output file(s) for {contig}...")
        for interval in tree:
            if interval.length() > max_len:
                logging.info(
                    f"Interval {interval.data} ({interval.length()}bp) is longer than "
                    f"the maximum allowed length. Skipping..."
                )
                continue

            contig_dir = outdir / contig
            contig_dir.mkdir(parents=True, exist_ok=True)
            filepath = contig_dir / f"{interval.data}.fa"
            if filepath.exists():
                raise FileExistsError(
                    f"A file already exists for {interval} at {filepath}"
                )
            header = (
                f">{interval.data} contig={contig}|start={interval.begin}|"
                f"end={interval.end}"
            )
            seq = slice_seq(index[contig], interval)
            filepath.write_text(f"{header}\n{seq}")

            interval_type = infer_interval_type(interval)
            print(
                ",".join(
                    map(
                        str,
                        [
                            "/".join(filepath.parts[-2:]),
                            interval_type,
                            interval.begin,
                            interval.end,
                            interval.data,
                            contig,
                        ],
                    )
                ),
                file=mapping_stream,
            )

            logging.debug(f"{interval} written to {filepath}")

    mapping_stream.close()
    logging.info(f"File mapping written to {file_mapping_path}")
    logging.info("All done!")


if __name__ == "__main__":
    main()
