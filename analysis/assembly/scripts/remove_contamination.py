import logging
import os
import sys
from pathlib import Path
from typing import Union, Set, NamedTuple

import click
import pysam

PathLike = Union[Path, str, os.PathLike]


class Classification(NamedTuple):
    read_id: str
    seq_id: str
    taxid: int
    score: float
    second_best_score: float
    hit_len: int
    query_len: int
    num_matches: int

    @staticmethod
    def from_line(line: str) -> "Classification":
        fields = line.strip().split("\t")
        read_id, seq_id = fields[:2]
        taxid = int(fields[2])
        score, second_best_score = map(float, fields[3:5])
        hit_len, query_len, num_matches = map(int, fields[5:])

        return Classification(
            read_id,
            seq_id,
            taxid,
            score,
            second_best_score,
            hit_len,
            query_len,
            num_matches,
        )


def extract_taxids_from_taxtree(taxtree: str) -> Set[int]:
    taxids = set()
    lines = map(str.strip, taxtree.split("\n"))
    for line in lines:
        fields = line.split()
        if fields:
            taxids.add(int(fields[0]))

    return taxids


@click.command()
@click.option(
    "-t",
    "--taxtree",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,),
    help=(
        "Taxon IDs NOT considered contamination. This file can be either a taxon ID on "
        "each line, or a taxonomic tree output by `taxonkit`."
    ),
)
@click.option(
    "-S",
    "--classification",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,),
    help="A centrifuge classification file for the input.",
)
@click.option(
    "-i",
    "--input",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True,),
    help="The fast{a,q} file to remove contamination from.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, writable=True, allow_dash=True
    ),
    help="The path to write the decontaminated sequences to.",
    default="-",
    show_default=True,
)
@click.option(
    "-v",
    "--invert",
    is_flag=True,
    help=(
        "Remove any record with a classification outside the taxtree. By default, any "
        "record with a classification within the taxtree is kept; regardless of "
        "whether it has off-taxtree classifications."
    ),
)
@click.help_option("--help", "-h")
@click.option("--verbose", help="Turns on debug-level logging.", is_flag=True)
@click.pass_context
def main(
    ctx: click.Context,
    classification: PathLike,
    taxtree: PathLike,
    input: PathLike,
    output: PathLike,
    verbose: bool,
    invert: bool,
):
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(format="[%(levelname)s] %(message)s", level=log_level)

    taxids_to_keep = extract_taxids_from_taxtree(Path(taxtree).read_text())

    records_to_keep = set()
    records_with_contamination = set()
    with Path(classification).open() as classification_fh:
        for line in classification_fh:
            if line.startswith("readID\tseqID"):  # skip header
                continue

            record_classification = Classification.from_line(line.strip())
            if record_classification.taxid in taxids_to_keep:
                records_to_keep.add(record_classification.read_id)
            else:
                records_with_contamination.add(record_classification.read_id)

    for read_id in records_to_keep.intersection(records_with_contamination):
        if not invert:
            logging.warning(
                (
                    f"{read_id} has classification(s) on and off the taxtree provided. "
                    "It will be kept - but you have been warned..."
                )
            )
        else:
            records_to_keep.remove(read_id)
            logging.info(
                (
                    f"Removing: {read_id} as it has both on- and off-taxtree "
                    "classifications."
                )
            )

    output_fh = sys.stdout if output == "-" else Path(output).open("w")

    logging.info("Removing contamination...")

    with pysam.FastxFile(input) as fastx_fh:
        for record in fastx_fh:
            if record.name in records_to_keep:
                logging.info(f"Keeping: {record.name}")
                print(str(record), file=output_fh)
            else:
                logging.info(f"Removing: {record.name}")

    logging.info("Finished!")


if __name__ == "__main__":
    main()
