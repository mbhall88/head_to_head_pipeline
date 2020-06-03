import logging
import re
from itertools import starmap

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import TextIO, Dict, Optional, NamedTuple, List

import click
from cyvcf2 import VCF, Variant

Indices = List[int]
LINEAGE_REGEX = re.compile(r"^(lineage)?(?P<major>[\w]+)\.?(?P<minor>.*)?$")


class RowError(Exception):
    pass


class InvalidLineageString(Exception):
    pass


class Genotype(NamedTuple):
    allele1: int
    allele2: int
    phased: bool

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


class PanelVariant(NamedTuple):
    """A class that represents a lineage-defining variant."""

    lineage: Lineage
    position: int
    ref: str
    alt: str
    locus_id: Optional[str] = None
    gene_name: Optional[str] = None
    gene_coord: Optional[int] = None

    @staticmethod
    def from_row(row: str, delim: str = ",") -> "PanelVariant":
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

        return PanelVariant(
            lineage, position, ref, alt, locus_id, gene_name, gene_coord
        )


class Classifier:
    """A class to handle validating and classifying VCF entries with respect to the
    panel."""

    def __init__(
        self, index: Optional[Dict[int, PanelVariant]] = None, include_het: bool = False
    ):
        if index is None:
            index = dict()
        self.index = index
        self.include_het = include_het

    def _is_position_valid(self, pos: int) -> bool:
        return pos in self.index

    def is_variant_valid(self, variant: Variant) -> bool:
        if not self._is_position_valid(variant.POS):
            return False

        panel_variant = self.index[variant.POS]
        if panel_variant.ref != variant.REF:
            logging.warning(
                f"Reference allele {variant.REF} at position {variant.POS} does "
                f"not match panel variant {panel_variant.ref} at that position. "
                f"Skipping this position..."
            )
            return False

        if panel_variant.alt not in variant.ALT:
            logging.debug(
                f"No samples contain the panel ALT variant at position "
                f"{variant.POS}. Skipping..."
            )
            return False

        return True

    def samples_with_lineage_variant(self, variant: Variant) -> Indices:
        """Returns the indices of samples in the variant that have the same alternate
        allele as the panel variant for that position."""
        indices = []
        panel_variant = self.index.get(variant.POS, None)
        if panel_variant is None:
            return indices

        for sample_idx, gt in enumerate(starmap(Genotype, variant.genotypes)):
            if gt.is_hom_alt():
                alt_idx = gt.alt_index()
                alt_base = variant.ALT[alt_idx]
                if alt_base == panel_variant.alt:
                    indices.append(sample_idx)
                    continue
            if self.include_het and gt.is_het():
                alt_idxs = [a - 1 for a in [gt.allele1, gt.allele2] if a > 0]
                alt_bases = {variant.ALT[i] for i in alt_idxs}
                if panel_variant.alt in alt_bases:
                    indices.append(sample_idx)

        return indices


def load_panel(
    stream: TextIO, no_header: bool = True, delim: str = ","
) -> Dict[int, PanelVariant]:
    """Creates a dictionary mapping a position to a Variant for each entry in the
    file that panel points to.
    """
    if not no_header:
        _ = next(stream)  # skip header

    index: Dict[int, PanelVariant] = dict()
    for row in map(str.rstrip, stream):
        if not row:
            continue
        variant = PanelVariant.from_row(row, delim=delim)
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
    "--default-lineage",
    help="Lineage to use when no panel variants are found for a sample.",
    default="unknown",
    show_default=True,
)
@click.option(
    "--output-delim",
    help="Delimiting character used in the output file.",
    default=",",
    show_default=True,
)
@click.option(
    "--no-header",
    help="Indicates there is no header line in the panel file.",
    is_flag=True,
)
@click.option(
    "--include-het",
    help="Consider heterozygous (with relevant ALT) as lineage-defining.",
    is_flag=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    input: str,
    panel: str,
    output: TextIO,
    delim: str,
    default_lineage: str,
    output_delim: str,
    no_header: bool,
    verbose: bool,
    include_het: bool,
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

    logging.info("Searching VCF for lineage-defining variants...")
    classification: Dict[str, List[Lineage]] = defaultdict(list)
    classifier = Classifier(panel_index, include_het=include_het)
    vcf = VCF(input)
    for variant in vcf:
        if not classifier.is_variant_valid(variant):
            continue

        samples_with_alt: Indices = classifier.samples_with_lineage_variant(variant)
        panel_variant = panel_index[variant.POS]
        for idx in samples_with_alt:
            sample = vcf.samples[idx]
            classification[sample].append(panel_variant.lineage)
            logging.debug(f"{sample} has variant for {panel_variant.lineage}")

    logging.info("Calling lineages for samples based on found lineage variants...")
    output_header = output_delim.join(["sample", "called_lineage", "found_lineages"])
    print(output_header, file=output)

    for sample in vcf.samples:
        if sample not in classification:
            # todo: what is a good default?
            logging.warning(
                f"No panel variants were found for {sample}. Defaulting to "
                f"{default_lineage}"
            )
            print(output_delim.join([sample, default_lineage, ""]), file=output)
            continue

        found_lineages = classification[sample]
        called_lineage = classifier.call_lineage(found_lineages)  # todo: implement method

    logging.info("All done!")
    vcf.close()


if __name__ == "__main__":
    main()
