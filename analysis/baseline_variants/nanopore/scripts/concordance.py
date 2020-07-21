import logging
from enum import Enum
from typing import TextIO, Set, Optional, NamedTuple, List, Tuple

import click
from cyvcf2 import Variant, VCF


class DummyVariant(NamedTuple):
    POS: float


DUMMY = DummyVariant(float("inf"))


class Bed:
    """Positions in this class are 1-based - in contrast to BED positions being 0-based.
    This is to align it with VCF positions which are 1-based."""

    def __init__(self, file: Optional[str] = None, zero_based: bool = True):
        self.zero_based = zero_based
        self.positions: Set[int] = set()
        if file is not None:
            with open(file) as instream:
                for line in map(str.rstrip, instream):
                    # start is 0-based inclusive; end is 0-based non-inclusive
                    start, end = [int(i) for i in line.split()[1:3]]
                    if not self.zero_based:
                        start += 1
                        end += 1
                    self.positions.update(range(start, end))

    def __contains__(self, item: int) -> bool:
        return item in self.positions


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

    @staticmethod
    def from_arr(arr: List[int]) -> "Genotype":
        alleles = [a for a in arr if type(a) is int]
        if len(alleles) < 2:
            alleles.append(-1)
        return Genotype(*alleles)


class Classification(Enum):
    Ref = "REF"
    Alt = "ALT"
    Null = "NULL"
    Het = "HET"
    Missing = "MISSING"

    def __str__(self):
        return self.value

    @staticmethod
    def from_variant(variant: Variant) -> "Classification":
        gt = Genotype.from_arr(variant.genotypes[0])
        if gt.is_null():
            classification = Classification.Null
        elif gt.is_het():
            classification = Classification.Het
        elif gt.is_hom_ref():
            classification = Classification.Ref
        elif gt.is_hom_alt():
            classification = Classification.Alt
        else:
            raise NotImplementedError(
                f"Can't determine genotype for variant \n{str(variant)}"
            )
        return classification


class Outcome(Enum):
    Null = "NULL"
    FalseNull = "FALSE_NULL"
    TrueRef = "TRUE_REF"
    FalseRef = "FALSE_REF"
    TrueAlt = "TRUE_ALT"
    FalseAlt = "FALSE_ALT"
    DiffAlt = "DIFF_ALT"
    Masked = "MASKED"
    BothFailFilter = "BOTH_FAIL_FILTER"
    AFailFilter = "A_FAIL_FILTER"
    BFailFilter = "B_FAIL_FILTER"
    Het = "HET"
    MissingPos = "MISSING_POS"

    def __str__(self):
        return self.value

    @staticmethod
    def from_variants(a_variant: Variant, b_variant: Variant) -> "Outcome":
        class_a = Classification.from_variant(a_variant)
        class_b = Classification.from_variant(b_variant)

        if class_a is Classification.Null or class_b is Classification.Null:
            return Outcome.Null if class_a is Classification.Null else Outcome.FalseNull

        if class_a is Classification.Het or class_b is Classification.Het:
            return Outcome.Het

        if class_b is Classification.Ref:
            return (
                Outcome.TrueRef if class_a is Classification.Ref else Outcome.FalseRef
            )

        if class_b is Classification.Alt:
            if class_a is not Classification.Alt:
                return Outcome.FalseAlt

            a_idx = Genotype.from_arr(a_variant.genotypes[0]).alt_index()
            b_idx = Genotype.from_arr(b_variant.genotypes[0]).alt_index()
            a_base = a_variant.ALT[a_idx]
            b_base = b_variant.ALT[b_idx]

            return Outcome.TrueAlt if a_base == b_base else Outcome.DiffAlt

        # can't think of a way we could get here...
        raise NotImplementedError(
            f"Could not classify variants at position {a_variant.POS}"
        )


class Classifier:
    def __init__(
        self, mask: Optional[Bed] = None, apply_filter: bool = False,
    ):
        self.apply_filter = apply_filter

        if mask is None:
            mask = Bed()
        self.mask = mask

    def classify(
        self, a_variant: Variant, b_variant: Variant
    ) -> Tuple[Classification, Classification, Outcome]:
        if a_variant.POS != b_variant.POS:
            raise IndexError(
                f"Expected positions of variant to match but got a: {a_variant.POS} "
                f"and b: {b_variant.POS}"
            )
        pos = a_variant.POS
        a_class = Classification.from_variant(a_variant)
        b_class = Classification.from_variant(b_variant)
        outcome = Outcome.from_variants(a_variant, b_variant)

        if pos in self.mask:
            outcome = Outcome.Masked
        elif self.apply_filter:
            if a_variant.FILTER and b_variant.FILTER:
                outcome = Outcome.BothFailFilter
            elif a_variant.FILTER:
                outcome = Outcome.AFailFilter
            elif b_variant.FILTER:
                outcome = Outcome.BFailFilter

        return a_class, b_class, outcome


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-a",
    "--truth-vcf",
    help="VCF file that is considered 'truth'.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-b",
    "--query-vcf",
    help="VCF file to compare to 'truth'.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-m",
    "--mask",
    "bedfile",
    help="BED file containing positions to ignore.",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "-c",
    "--csv",
    type=click.File(mode="w", lazy=True),
    default="-",
    show_default=True,
    help="Write the classifications for each position to a CSV file.",
)
@click.option(
    "--apply-filter/--no-apply-filter",
    "-f/-F",
    default=True,
    show_default=True,
    help="Adhere to the VCF FILTER column.",
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    truth_vcf: str,
    query_vcf: str,
    bedfile: str,
    csv: TextIO,
    apply_filter: bool,
    verbose: bool,
):
    """Calculate concordance between two VCF files. An assumption is made that both VCFs
    have the exact same positions present."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    if bedfile:
        logging.info("Loading mask...")
        # make 1-based as this is the same as VCF
        mask = Bed(bedfile, zero_based=False)
        logging.info(f"Loaded {len(mask.positions)} positions to mask.")
    else:
        logging.info("No mask given. Classifying all positions...")
        mask = Bed()

    classifier = Classifier(mask=mask, apply_filter=apply_filter)

    logging.info("Classifying variants...")
    print("pos,a,b,classification", file=csv)
    a_vcf = VCF(truth_vcf)
    a_variant = next(a_vcf, DUMMY)
    a_exhausted = a_variant.POS == float("inf")
    b_vcf = VCF(query_vcf)
    b_variant = next(b_vcf, DUMMY)
    b_exhausted = b_variant.POS == float("inf")
    data = []

    while (not a_exhausted) and (not b_exhausted):
        if a_variant.POS == b_variant.POS:
            class_a, class_b, outcome = classifier.classify(a_variant, b_variant)
            row = [a_variant.POS, class_a, class_b, outcome]
            data.append(row)
            print(",".join(map(str, row)), file=csv)
            a_variant = next(a_vcf, DUMMY)
            a_exhausted = a_variant.POS == float("inf")
            b_variant = next(b_vcf, DUMMY)
            b_exhausted = b_variant.POS == float("inf")
            continue

        pos = min(a_variant.POS, b_variant.POS)
        class_a = (
            Classification.from_variant(a_variant)
            if pos == a_variant.POS
            else Classification.Missing
        )
        class_b = (
            Classification.from_variant(b_variant)
            if pos == b_variant.POS
            else Classification.Missing
        )
        classification = Outcome.Masked if pos in mask else Outcome.MissingPos
        row = [a_variant.POS, class_a, class_b, classification]
        data.append(row)
        print(",".join(map(str, row)), file=csv)

        if pos == a_variant.POS:
            a_variant = next(a_vcf, DUMMY)
            a_exhausted = a_variant.POS == float("inf")
        elif pos == b_variant.POS:
            b_variant = next(b_vcf, DUMMY)
            b_exhausted = b_variant.POS == float("inf")
        else:
            raise NotImplementedError(
                f"Failed to compare the following two variants:\n{str(a_variant)}\n"
                f"{str(b_variant)}"
            )


if __name__ == "__main__":
    main()
