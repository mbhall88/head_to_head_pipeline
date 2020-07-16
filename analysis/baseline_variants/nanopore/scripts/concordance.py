import logging
from enum import Enum
from typing import TextIO, Set, Optional
from cyvcf2 import VCF, Variant

import click

class Bed:
    """Positions in this class are 1-based - in contrast to BED positions being 0-based.
    This is to align it with VCF positions which are 1-based."""
    def __init__(self, file: Optional[str]=None, zero_based: bool=True):
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


class Classifier:
    def __init__(self, mask: Optional[Bed] = None, skip_null: bool=False):
        self.skip_null = skip_null

        if mask is None:
            mask = Bed()
        self.mask = mask

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

@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-a", "--truth-vcf",
    help="VCF file that is considered 'truth'.",
    type=click.Path(exists=True, dir_okay=False),
    required=True
)
@click.option(
    "-b", "--query-vcf",
    help="VCF file to compare to 'truth'.",
    type=click.Path(exists=True, dir_okay=False),
    required=True
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
    "-N", "--skip-null",
    help="If either VCF has a NULL call at a position, skip classification.",
    is_flag=True,
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    truth_vcf: str,
    query_vcf: str,
    bedfile: str,
    csv: TextIO,
    skip_null: bool,
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

    classifier = Classifier(mask=mask, skip_null=skip_null)



if __name__ == "__main__":
    main()
