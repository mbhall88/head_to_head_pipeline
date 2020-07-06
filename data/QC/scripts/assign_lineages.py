import logging
import re
from collections import Counter, defaultdict
from enum import Enum
from itertools import starmap, repeat, takewhile
from typing import TextIO, Dict, Optional, NamedTuple, List, Union

import click
from cyvcf2 import VCF, Variant

Indices = List[int]
LINEAGE_DELIM = "."
LINEAGE_REGEX = re.compile(
    rf"^(lineage)?(?P<major>[\w]+)\{LINEAGE_DELIM}?(?P<minor>.*)?$"
)


class RowError(Exception):
    pass


class InvalidLineageString(Exception):
    pass


class Filter(Enum):
    Pass = "PASS"
    Fail = "FAIL"

    @staticmethod
    def from_str(s: str) -> "Filter":
        return Filter.Pass if s in {"PASS", "."} else Filter.Fail


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


class Lineage:
    """A class representing a lineage and any of its sublineage information."""

    def __init__(
        self,
        major: str,
        minor: Optional[Union[str, List[str]]] = None,
        minor_delim: str = ".",
    ):
        self.major = major
        if minor is None:
            minor = tuple()
        if isinstance(minor, str):
            minor = tuple(minor.split(minor_delim))
        self.minor = tuple(minor)
        self._minor_delim = minor_delim

    def __str__(self):
        if not self.minor:
            return self.major
        return self.major + self._minor_delim + self._minor_delim.join(self.minor)

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

    def __eq__(self, other: "Lineage") -> bool:
        return self.major == other.major and self.minor == other.minor

    def __lt__(self, other: "Lineage") -> bool:
        if (not self.minor and not other.minor) or not self.minor:
            return False
        if not other.minor:
            return True
        # we now know that both have minors
        return len(self.minor) > len(other.minor)

    def mrca(self, other: "Lineage") -> Optional["Lineage"]:
        """Determine the most recent common ancestor between two Lineages.
        Returns None if there is no MRCA.
        """
        if self.major != other.major:
            return None
        if not self.minor or not other.minor:
            return Lineage(major=self.major)

        try:
            common_minors, _ = zip(
                *takewhile(lambda xy: xy[0] == xy[1], zip(self.minor, other.minor))
            )
            minor_str = LINEAGE_DELIM.join(common_minors)
        except ValueError:
            minor_str = None

        return Lineage(self.major, minor_str)

    @staticmethod
    def call(lineages: List["Lineage"]) -> Optional["Lineage"]:
        """Returns the Lineage with the most specific minor.
        if the majors are different, returns None. If minors are the same, then takes
        MRCA of the lineages with the same minor length.
        """
        if not lineages:
            return None
        if len(lineages) == 1:
            return lineages[0]
        lineages.sort()
        minors_of_same_len = filter(
            lambda l: len(l.minor) == len(lineages[0].minor), lineages
        )
        lineage = lineages[0]
        for lin in minors_of_same_len:
            lineage = lin.mrca(lineage)

        return lineage


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
        self,
        index: Optional[Dict[int, PanelVariant]] = None,
        max_het: int = 0,
        max_alt_lineages: int = 0,
        ref_lineage_position: int = 0,
    ):
        if index is None:
            index = dict()
        self.index = index
        self.max_het = max_het
        self.het_counts: Dict[int, int] = defaultdict(int)
        self.max_alt_lineages = max_alt_lineages
        self.ref_lineage_position = ref_lineage_position

    def _is_position_valid(self, pos: int) -> bool:
        return pos in self.index

    def is_variant_valid(self, variant: Variant) -> bool:
        if not self._is_position_valid(variant.POS):
            return False

        failed_filter = variant.FILTER is not None
        if failed_filter:
            logging.debug(f"Position {variant.POS} failed FILTER with {variant.FILTER}")
            return False

        panel_variant = self.index[variant.POS]
        if variant.POS == self.ref_lineage_position:
            if panel_variant.alt != variant.REF:
                logging.warning(
                    f"Panel alternate base does not match REF at reference lineage "
                    f"position {variant.POS}."
                )
                return False
            return True

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

        try:
            filters = map(Filter.from_str, variant.format("FT"))
        except KeyError:
            filters = repeat(Filter.Pass)

        for sample_idx, gt in enumerate(starmap(Genotype, variant.genotypes)):
            failed_filter = next(filters) is Filter.Fail
            if failed_filter:
                continue

            if variant.POS == self.ref_lineage_position and gt.is_hom_ref():
                indices.append(sample_idx)
                continue

            if gt.is_hom_alt():
                alt_idx = gt.alt_index()
                alt_base = variant.ALT[alt_idx]
                if alt_base == panel_variant.alt:
                    indices.append(sample_idx)
                    continue
            if gt.is_het():
                alt_idxs = [a - 1 for a in [gt.allele1, gt.allele2] if a > 0]
                alt_bases = {variant.ALT[i] for i in alt_idxs}
                if panel_variant.alt in alt_bases:
                    indices.append(sample_idx)
                    self.het_counts[sample_idx] += 1

        return indices

    def call_sample_lineage(
        self, lineages: List[Lineage], sample_idx: int, default: str = ""
    ) -> str:
        if not lineages:
            return default

        num_hets = self.het_counts[sample_idx]
        if num_hets > self.max_het:
            return "too_many_hets"

        if len(lineages) == 1:
            return str(lineages[0])

        majors = Counter([lin.major for lin in lineages])
        two_most_common = majors.most_common(n=2)
        most_common_count = two_most_common[0][-1]
        most_common_major = two_most_common[0][0]

        even_count_on_two_most_common_majors = (
            len(two_most_common) > 1 and most_common_count == two_most_common[-1][-1]
        )
        if even_count_on_two_most_common_majors:
            return "mixed"

        num_alt_lineages = len(lineages) - most_common_count
        if num_alt_lineages > self.max_alt_lineages:
            return "mixed"

        non_alt_lineages = [
            lin for lin in lineages if lin.major == most_common_major[0]
        ]
        return str(Lineage.call(non_alt_lineages))


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
    "--ref-lineage-position",
    help=(
        "Variant position that defines the lineage of the VCF reference. When "
        "classifying lineages a REF call at this position will be considered a call "
        "for the reference lineage as opposed to other variants which require an ALT "
        "call. Set to 0 to disable this option."
    ),
    default=0,
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
    help="Delimiting character used in the output file. Do not use `;`.",
    default=",",
    show_default=True,
)
@click.option(
    "--no-header",
    help="Indicates there is no header line in the panel file.",
    is_flag=True,
)
@click.option(
    "--max-het",
    help=(
        "Maximum allowed heterozygous lineage-defining variants before abandoning "
        "lineage assignment for a sample."
    ),
    default=0,
    show_default=True,
)
@click.option(
    "--max-alt-lineages",
    help=(
        "Maximum allowed number of variants from different major lineages. For example "
        "if a sample has 2 L4 variants and 1 L3 variant, this sample would be called "
        "L4 if this parameter is set to 1 or 'mixed' if set to 0"
    ),
    default=0,
    show_default=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    input: str,
    panel: str,
    output: TextIO,
    delim: str,
    default_lineage: str,
    ref_lineage_position: int,
    output_delim: str,
    no_header: bool,
    verbose: bool,
    max_het: int,
    max_alt_lineages: int,
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
    classifier = Classifier(
        panel_index,
        max_het=max_het,
        max_alt_lineages=max_alt_lineages,
        ref_lineage_position=ref_lineage_position,
    )
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
    output_header = output_delim.join(
        ["sample", "major_lineage", "full_lineage", "found_lineages"]
    )
    print(output_header, file=output)

    for idx, sample in enumerate(vcf.samples):
        if sample not in classification:
            # todo: what is a good default?
            logging.warning(
                f"No panel variants were found for {sample}. Defaulting to "
                f"{default_lineage}"
            )
            print(
                output_delim.join([sample, default_lineage, default_lineage, ""]),
                file=output,
            )
            continue

        found_lineages = classification[sample]
        called_lineage = Lineage.from_str(
            classifier.call_sample_lineage(
                found_lineages, sample_idx=idx, default=default_lineage
            )
        )
        print(
            output_delim.join(
                [
                    sample,
                    called_lineage.major,
                    str(called_lineage),
                    ";".join(map(str, found_lineages)),
                ]
            ),
            file=output,
        )

    logging.info("All done!")
    vcf.close()


if __name__ == "__main__":
    main()
