import logging
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Optional, NamedTuple, List, Tuple

import click
import numpy as np
from cyvcf2 import Variant, VCF, Writer

HIST_BINS = 40


class Tags(Enum):
    Depth = "DP"
    LowDepth = "ld"
    HighDepth = "hd"
    LowQual = "lq"
    StrandBias = "sb"
    StrandDepth = "DP4"
    Pass = "PASS"

    def __str__(self) -> str:
        return str(self.value)


@dataclass
class FilterStatus:
    low_depth: bool = False
    high_depth: bool = False
    low_qual: bool = False
    strand_bias: bool = False
    delim: str = ";"

    def __str__(self) -> str:
        status = []
        if self.low_depth:
            status.append(str(Tags.LowDepth))
        if self.high_depth:
            status.append(str(Tags.HighDepth))
        if self.low_qual:
            status.append(str(Tags.LowQual))
        if self.strand_bias:
            status.append(str(Tags.StrandBias))

        return self.delim.join(status) if status else str(Tags.Pass)


@dataclass
class StrandDepths:
    ref_forward: int = 0
    ref_reverse: int = 0
    alt_forward: int = 0
    alt_reverse: int = 0

    @staticmethod
    def _ratio(depths: Tuple[int, int]) -> float:
        try:
            return min(depths) / sum(depths)
        except ZeroDivisionError:
            return 1.0

    @property
    def ref_depths(self) -> Tuple[int, int]:
        return self.ref_forward, self.ref_reverse

    @property
    def alt_depths(self) -> Tuple[int, int]:
        return self.alt_forward, self.alt_reverse

    @property
    def ref_ratio(self) -> float:
        return self._ratio(self.ref_depths)

    @property
    def alt_ratio(self) -> float:
        return self._ratio(self.alt_depths)

    def to_tuple(self) -> Tuple[int, int, int, int]:
        return self.ref_forward, self.ref_reverse, self.alt_forward, self.alt_reverse


class DepthTagError(Exception):
    pass


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


class Filter:
    def __init__(
        self,
        expected_depth: int = 0,
        min_depth: float = 0,
        max_depth: float = 0,
        min_strand_bias: int = 0,
        min_qual: float = 0,
    ):
        self.expected_depth = expected_depth
        self.min_depth_frac = min_depth
        self.min_depth = self.expected_depth * self.min_depth_frac
        self.max_depth_frac = max_depth
        self.max_depth = self.expected_depth * self.max_depth_frac
        self.min_strand_bias = min_strand_bias / 100
        self.min_qual = min_qual

        if self.min_depth and self.max_depth and self.min_depth > self.max_depth:
            raise ValueError(
                f"Minimum depth is more than maximum depth: "
                f"{self.min_depth:.1f} > {self.max_depth:.1f}"
            )

    def _is_low_depth(self, variant: Variant) -> bool:
        variant_depth = get_depth(variant)
        return variant_depth < self.min_depth if self.min_depth else False

    def _is_high_depth(self, variant: Variant) -> bool:
        variant_depth = get_depth(variant)
        return variant_depth > self.max_depth if self.max_depth else False

    def _is_low_qual(self, variant: Variant) -> bool:
        return variant.QUAL < self.min_qual

    def filter_status(self, variant: Variant) -> str:
        status = FilterStatus()
        if self.min_depth or self.max_depth:
            status.low_depth = self._is_low_depth(variant)
            status.high_depth = self._is_high_depth(variant)

        if self.min_qual:
            status.low_qual = self._is_low_qual(variant)

        if self.min_strand_bias:
            strand_depths = get_strand_depths(variant)
            assert strand_depths is not None, (
                f"Strand bias filter should be turned off if no {Tags.StrandDepth} "
                f"tag is present."
            )

            gt = Genotype.from_arr(variant.genotypes[0])
            if gt.is_hom_alt():
                ratio = strand_depths.alt_ratio
            elif gt.is_hom_ref():
                ratio = strand_depths.ref_ratio
            elif gt.is_het():
                ratio = min(strand_depths.ref_ratio, strand_depths.alt_ratio)
            elif gt.is_null():
                ratio = float("inf")
            else:
                raise NotImplementedError(
                    f"Don't know how to interpret genotype {gt} for variant at "
                    f"POS {variant.POS}"
                )
            status.strand_bias = ratio < self.min_strand_bias

        return str(status)

    def add_filters_to_header(self, vcf: VCF):
        if self.min_depth > 0:
            header = {
                "ID": str(Tags.LowDepth),
                "Description": (
                    f"Depth ({Tags.Depth}) less than {self.min_depth_frac:.1%} the "
                    f"expected depth of {self.expected_depth:.1f}. "
                    f"{Tags.Depth}<{self.min_depth:.1f}"
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. depth: {header}")

        if self.max_depth > 0:
            header = {
                "ID": str(Tags.HighDepth),
                "Description": (
                    f"Depth ({Tags.Depth}) more than {self.max_depth_frac:.1%} the "
                    f"expected depth of {self.expected_depth:.1f}. "
                    f"{Tags.Depth}>{self.max_depth:.1f}"
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for max. depth: {header}")

        if self.min_qual > 0:
            header = {
                "ID": str(Tags.LowQual),
                "Description": f"QUAL less than {self.min_qual}",
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for min. QUAL: {header}")

        if self.min_strand_bias > 0:
            header = {
                "ID": str(Tags.StrandBias),
                "Description": (
                    f"A strand on the called allele has less than  "
                    f"{self.min_strand_bias:.2%} of the high-quality depth for that "
                    f"allele. This is judged on the {Tags.StrandDepth} tag."
                ),
            }
            vcf.add_filter_to_header(header)
            logging.debug(f"Header for strand bias: {header}")


def get_depth(variant: Variant, default: int = 0) -> int:
    return variant.INFO.get(str(Tags.Depth), default)


def get_strand_depths(
    variant: Variant, default: Optional[StrandDepths] = None
) -> Optional[StrandDepths]:
    strand_depths = variant.INFO.get(str(Tags.StrandDepth), None)
    return StrandDepths(*strand_depths) if strand_depths is not None else default


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--in-vcf",
    help="VCF file to apply filters to.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-o",
    "--out-vcf",
    help="New VCF with filters.",
    type=click.Path(exists=False, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-d",
    "--min-depth",
    help=(
        "Minimum depth as a percentage of the expected (median) depth. This filter "
        f"has ID: {Tags.LowDepth}. Set to 0 to disable"
    ),
    default=0.2,
    show_default=True,
)
@click.option(
    "-D",
    "--max-depth",
    help=(
        "Maximum depth as a fraction of the expected (median) depth. This filter "
        f"has ID: {Tags.HighDepth}. Set to 0 to disable"
    ),
    default=2.0,
    show_default=True,
)
@click.option(
    "-s",
    "--min-strand-bias",
    help=(
        "Filter a variant if either strand has less than INT% of the (high-quality) "
        f"depth on the called allele ({Tags.StrandDepth}). This filter has ID: "
        f"{Tags.StrandBias}. Set to 0 to disable"
    ),
    type=click.IntRange(0, 50),
    metavar="INT",
    default=25,
    show_default=True,
)
@click.option(
    "-q",
    "--min-qual",
    help=(
        f"Filter a variant if QUAL is less than INT. This filter has ID: "
        f"{Tags.LowQual}. Set to 0 to disable"
    ),
    default=0.0,
    show_default=True,
)
@click.option(
    "--overwrite/--no-overwrite",
    "-f/-F",
    default=True,
    show_default=True,
    help="Overwrite existing information in FILTER field.",
)
@click.option(
    "-p/-P",
    "--hist/--no-hist",
    default=True,
    show_default=True,
    help="Print histograms of depth and QUAL. Histograms will be printed to stderr.",
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    in_vcf: str,
    out_vcf: str,
    overwrite: bool,
    verbose: bool,
    min_qual: float,
    min_depth: float,
    max_depth: float,
    min_strand_bias: int,
    hist: bool,
):
    """Apply the following filters to a VCF:\n
      - Minimum proportion of the expected (median) depth\n
      - Maximum proportion of the expected (median) depth\n
      - Minimum QUAL threshold\n
      - Minimum Strand bias percentage
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    vcf_reader = VCF(in_vcf)
    if not vcf_reader.contains(Tags.Depth.value):
        raise DepthTagError(f"Depth tag {Tags.Depth} not found in header")

    if (not vcf_reader.contains(str(Tags.StrandDepth))) and min_strand_bias:
        logging.warning(
            f"Strand depth tag {Tags.StrandDepth} not found in header. "
            f"Turning off strand bias filter..."
        )
        min_strand_bias = 0

    logging.info("Calculating expected (median) depth...")
    depths = []
    quals = []
    for v in vcf_reader:
        depths.append(get_depth(v))
        quals.append(v.QUAL)

    median_depth = np.median(depths)
    logging.info(f"Expected depth: {median_depth}")

    median_qual = np.median(quals)
    logging.info(f"Median QUAL: {median_qual}")

    if hist:
        import histoprint

        tick_format = "% .1f"
        logging.info("Depth histogram:")
        histoprint.print_hist(
            np.histogram(depths, bins=HIST_BINS),
            title="Depth histogram",
            summary=True,
            tick_format=tick_format,
            file=click.get_text_stream("stderr"),
        )

        logging.info("QUAL histogram")
        histoprint.print_hist(
            np.histogram(quals, bins=HIST_BINS),
            title="QUAL histogram",
            summary=True,
            tick_format=tick_format,
            file=click.get_text_stream("stderr"),
        )

    assessor = Filter(
        expected_depth=int(median_depth),
        min_qual=min_qual,
        min_depth=min_depth,
        max_depth=max_depth,
        min_strand_bias=min_strand_bias,
    )

    vcf_reader = VCF(in_vcf)
    assessor.add_filters_to_header(vcf_reader)
    vcf_writer = Writer(out_vcf, tmpl=vcf_reader)

    stats = Counter()
    logging.info("Filtering variants...")
    for variant in vcf_reader:
        filter_status = assessor.filter_status(variant)

        if (
            (not overwrite)
            and variant.FILTER is not None
            and filter_status != str(Tags.Pass)
        ):
            current_filter = variant.FILTER.rstrip(";")
            variant.FILTER = f"{current_filter};{filter_status}"
        else:
            variant.FILTER = filter_status

        vcf_writer.write_record(variant)

        stats.update(filter_status.split(";"))

    vcf_reader.close()
    vcf_writer.close()

    logging.info("FILTER STATISTICS")
    logging.info("=================")
    for filt, count in stats.items():
        logging.info(f"Filter: {filt}\tCount: {count}")

    logging.info("Done!")


if __name__ == "__main__":
    main()
