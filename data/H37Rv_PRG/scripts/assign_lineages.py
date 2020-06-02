import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import TextIO, Dict, Union, Optional, NamedTuple

import click
from cyvcf2 import VCF

LINEAGE_REGEX = re.compile(r"^(lineage)?(?P<major>[\w]+)\.?(?P<minor>.*)?$")


class RowError(Exception):
    pass


class InvalidLineageString(Exception):
    pass


@dataclass
class Lineage:
    """A class representing a lineage and any of its sublineage information."""

    major: str
    minor: Optional[str] = None

    @staticmethod
    def from_str(s: str) -> "Lineage":
        """Create a Lineage object from a str."""
        match = LINEAGE_REGEX.search(s)
        if not match:
            raise InvalidLineageString(
                f"Lineage string {s} is not in the expected format."
            )
        major = match.group("major")
        minor = match.group("minor") or None
        return Lineage(major=major, minor=minor)


class Variant(NamedTuple):
    """A class that represents a lineage-defining variant."""

    lineage: Lineage
    position: int
    ref: str
    alt: str
    locus_id: Optional[str] = None
    gene_name: Optional[str] = None
    gene_coord: Optional[int] = None

    @staticmethod
    def from_row(row: str, delim: str = ",") -> "Variant":
        """A Variant factory constructor to build a Variant object from a row in a
        delimited file.
        """
        fields = row.split(delim)
        if not fields or not any(fields):
            raise RowError("Row is empty.")

        lineage = Lineage.from_str(fields[0])
        position = int(fields[1])
        ref, alt = fields[3].split("/")
        try:
            locus_id = fields[7] or None
        except IndexError:
            locus_id = None

        try:
            gene_name = fields[8] or None
        except IndexError:
            gene_name = None

        try:
            gene_coord = int(fields[2])
        except (IndexError, ValueError):
            gene_coord = None

        return Variant(lineage, position, ref, alt, locus_id, gene_name, gene_coord)


def load_panel(
    stream: TextIO, no_header: bool = True, delim: str = ","
) -> Dict[int, Variant]:
    """Creates a dictionary mapping a position to a Variant for each entry in the
    file that panel points to.
    """
    if not no_header:
        _ = next(stream)  # skip header

    index: Dict[int, Variant] = dict()
    for row in map(str.rstrip, stream):
        if not row:
            continue
        variant = Variant.from_row(row, delim=delim)
        if variant.position in index:
            raise IndexError(f"Duplicate position {variant.position} in panel.")
        index[variant.position] = variant

    return index


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--input",
    help="VCF file to call lineages for.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-p",
    "--panel",
    help="Panel containing variants that define the lineages.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.File(mode="w", lazy=True),
    default="-",
    show_default=True,
    help="The filepath to write the output to.",
)
@click.option(
    "-d",
    "--delim",
    help="Delimiting character used in the panel file.",
    default=",",
    show_default=True,
)
@click.option(
    "--no-header",
    help="Indicates there is no header line in the panel file.",
    is_flag=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    input: str, panel: str, output: TextIO, delim: str, no_header: bool, verbose: bool,
):
    """Call Mycobacterium tuberculosis lineages for samples in a VCF file."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    logging.info("Loading the panel...")
    with open(panel) as stream:
        panel_index = load_panel(stream, no_header=no_header, delim=delim)
    logging.info(f"Loaded panel with {len(panel_index)} variants successfully.")

    classification: Dict[str, Counter] = defaultdict(Counter)
    with VCF(input) as vcf:
        for variant in vcf:
            if variant.POS not in panel_index:
                continue

            panel_variant = panel_index[variant.POS]
            if panel_variant.ref != variant.REF:
                logging.warning(
                    f"Reference allele {variant.REF} at position {variant.POS} does "
                    f"not match panel variant {panel_variant.ref} at that position. "
                    f"Skipping this position..."
                )
                continue

            try:
                alt_idx = variant.ALT.index(panel_variant.alt)
            except ValueError:
                logging.debug(
                    f"No samples contain the panel ALT variant at position "
                    f"{variant.POS}"
                )
                continue

            # todo: get genotypes and see if any are the alt we want


if __name__ == "__main__":
    main()
