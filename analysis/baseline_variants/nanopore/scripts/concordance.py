import logging
from enum import Enum
from typing import TextIO, Set, Optional, NamedTuple, List

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


class Outcome(Enum):
    TrueNull = "TRUE_NULL"
    FalseNull = "FALSE_NULL"
    TrueRef = "TRUE_REF"
    FalseRef = "FALSE_REF"
    TrueAlt = "TRUE_ALT"
    FalseAlt = "FALSE_ALT"
    DiffAlt = "DIFF_ALT"
    Masked = "MASKED"
    SkippedNull = "SKIPPED_NULL"
    BothFailFilter = "BOTH_FAIL_FILTER"
    AFailFilter = "A_FAIL_FILTER"
    BFailFilter = "B_FAIL_FILTER"
    BothHet = "BOTH_HET"
    AHet = "A_HET"
    BHet = "B_HET"
    AMissingPos = "A_MISSING_POS"
    BMissingPos = "B_MISSING_POS"


class Classifier:
    def __init__(
        self,
        mask: Optional[Bed] = None,
        skip_null: bool = False,
        apply_filter: bool = False,
    ):
        self.skip_null = skip_null
        self.apply_filter = apply_filter

        if mask is None:
            mask = Bed()
        self.mask = mask

    def classify(self, a_variant: Variant, b_variant: Variant) -> Outcome:
        if a_variant.POS != b_variant.POS:
            raise IndexError(
                f"Expected positions of variant to match but got a: {a_variant.POS} "
                f"and b: {b_variant.POS}"
            )
        pos = a_variant.POS

        if pos in self.mask:
            return Outcome.Masked

        if self.apply_filter:
            if a_variant.FILTER and b_variant.FILTER:
                return Outcome.BothFailFilter
            if a_variant.FILTER:
                return Outcome.AFailFilter
            if b_variant.FILTER:
                return Outcome.BFailFilter

        a_gt = Genotype.from_arr(a_variant.genotypes[0])
        b_gt = Genotype.from_arr(b_variant.genotypes[0])
        if self.skip_null and (a_gt.is_null() or b_gt.is_null()):
            return Outcome.SkippedNull

        if b_gt.is_null():
            return Outcome.TrueNull if a_gt.is_null() else Outcome.FalseNull

        if a_gt.is_het() and b_gt.is_het():
            return Outcome.BothHet
        elif a_gt.is_het():
            return Outcome.AHet
        elif b_gt.is_het():
            return Outcome.BHet

        if b_gt.is_hom_ref():
            return Outcome.TrueRef if a_gt.is_hom_ref() else Outcome.FalseRef

        if b_gt.is_hom_alt():
            if not a_gt.is_hom_alt():
                return Outcome.FalseAlt

            a_base = a_variant.ALT[a_gt.alt_index()]
            b_base = b_variant.ALT[b_gt.alt_index()]

            return Outcome.TrueAlt if a_base == b_base else Outcome.DiffAlt

        # can't think of a way we could get here...
        raise NotImplementedError(f"Could not classify variants at position {pos}")


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
    "-N",
    "--skip-null",
    help="If either VCF has a NULL call at a position, skip classification.",
    is_flag=True,
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
    skip_null: bool,
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

    classifier = Classifier(mask=mask, skip_null=skip_null, apply_filter=apply_filter)

    logging.info("Classifying variants...")
    print("pos,classification", file=csv)
    a_vcf = VCF(truth_vcf)
    a_variant = next(a_vcf, DUMMY)
    a_exhausted = a_variant.POS == float("inf")
    b_vcf = VCF(query_vcf)
    b_variant = next(b_vcf, DUMMY)
    b_exhausted = b_variant.POS == float("inf")

    while (not a_exhausted) and (not b_exhausted):
        if a_variant.POS == b_variant.POS:
            classification = classifier.classify(a_variant, b_variant)
            print(f"{a_variant.POS},{classification.value}", file=csv)
            a_variant = next(a_vcf, DUMMY)
            a_exhausted = a_variant.POS == float("inf")
            b_variant = next(b_vcf, DUMMY)
            b_exhausted = b_variant.POS == float("inf")
            continue

        if a_variant.POS > b_variant.POS:
            classification = (
                Outcome.Masked if b_variant.POS in mask else Outcome.AMissingPos
            )
            print(f"{b_variant.POS},{classification.value}", file=csv)
            b_variant = next(b_vcf, DUMMY)
            b_exhausted = b_variant.POS == float("inf")
            continue

        if a_variant.POS < b_variant.POS:
            classification = (
                Outcome.Masked if a_variant.POS in mask else Outcome.BMissingPos
            )
            print(f"{a_variant.POS},{classification.value}", file=csv)
            a_variant = next(a_vcf, DUMMY)
            a_exhausted = a_variant.POS == float("inf")
            continue


if __name__ == "__main__":
    main()
