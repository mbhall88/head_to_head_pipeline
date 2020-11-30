import logging
from typing import TextIO, Dict, List, Set

import click
import pandas as pd
from cyvcf2 import VCF

LOCI_ID = "LOCI"
LOCI_POS_ID = "LPOS"


def add_info_headers(vcf: VCF):
    vcf.add_info_to_header(
        {
            "ID": LOCI_ID,
            "Number": "1",
            "Type": "String",
            "Description": "CHROM (loci) in the original VCF",
        }
    )
    vcf.add_info_to_header(
        {
            "ID": LOCI_POS_ID,
            "Number": "1",
            "Type": "Integer",
            "Description": f"{LOCI_ID} (see other INFO header) POS in the original VCF",
        }
    )


def reset_contigs_in_header(header: str, contigs: Set[str]) -> str:
    header_without_contigs = [
        row for row in header.rstrip().split("\n") if not row.startswith("##contig")
    ]
    new_contig_rows = "\n".join([f"##contig=<ID={s}>" for s in contigs])
    column_header = header_without_contigs[-1]
    header = "\n".join(header_without_contigs[:-1])
    header += f"\n{new_contig_rows}"
    header += f"\n{column_header}"
    return header


def load_reference(path: str) -> Dict[str, List[str]]:
    index = dict()
    with open(path) as instream:
        seqid = None
        for line in map(str.rstrip, instream):
            if line.startswith(">"):
                seqid = line.split()[0][1:]
                if seqid in index:
                    raise ReferenceError(
                        f"Found duplicate sequence identifier in the reference fasta: "
                        f"{seqid}"
                    )
                index[seqid] = []
                continue
            if seqid is None:
                raise ReferenceError(
                    "Did not find header before sequence in reference fasta"
                )
            index[seqid].extend(list(line))

    return index


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--in-vcf",
    help="VCF file to normalise.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-o",
    "--out-vcf",
    help=(
        "Path to write normalised VCF to. Note: can only output VCF currently - not BCF"
    ),
    type=click.File(mode="w"),
    default="-",
    show_default=True,
)
@click.option(
    "-l",
    "--loci-info",
    help=(
        "CSV file linking loci to the original genome. The start and end columns are "
        "BOTH assumed to be 0-based inclusive"
    ),
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-r",
    "--ref",
    help=(
        "Original reference genome that positions are being normalised to. If given, "
        "the reference sequence in the input VCF will be checked to ensure it matches "
        "the position and chromosome it is being normalised to."
    ),
    type=click.Path(exists=True, dir_okay=False),
)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    in_vcf: str, out_vcf: TextIO, loci_info: str, verbose: bool, ref: str,
):
    """Normalise the positions in a pandora VCF. As pandora VCF positions are with
    respect to the local PRG, they cannot be easily compared to VCFs from other tools
    that call variants with respect to the reference genome. Provided the --vcf-ref
    option was used in pandora, and the vcf reference sequence used is the same as the
    reference genome, this script will output a VCF where the POS column relates to the
    original reference genome instead of the local PRG. Information about the local PRG
    is moved to the INFO field - so no information is lost.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    ref_index = dict()
    if ref:
        logging.info("Loading reference fasta...")
        ref_index = load_reference(ref)
        logging.info(f"Loaded {len(ref_index)} sequence(s) from the reference fasta")

    logging.info("Loading loci information...")
    loci_df = pd.read_csv(loci_info, index_col="name")

    logging.info(f"Loaded information for {len(loci_df)} loci")

    vcf_reader = VCF(in_vcf)
    add_info_headers(vcf_reader)
    vcf_header = vcf_reader.raw_header
    vcf_header = reset_contigs_in_header(vcf_header, set(loci_df.contig))
    print(vcf_header.rstrip(), file=out_vcf)

    logging.info("Normalising positions...")
    for variant in vcf_reader:
        if variant.CHROM not in loci_df.index:
            raise KeyError(f"{variant.CHROM} is not in the loci information file")

        loci = loci_df.loc[variant.CHROM]
        loci_len = loci.end - loci.start
        if variant.POS > loci_len:
            raise IndexError(
                f"Position {variant.POS} for loci {variant.CHROM} is outside the "
                f"expected length of the loci"
            )

        if ref:
            if loci.contig not in ref_index:
                raise IndexError(f"Contig {loci.contig} not in the reference genome")

            ref_start = loci.start + (variant.POS - 1)  # 0-based inclusive
            ref_stop = ref_start + len(variant.REF)  # 0-based non-inclusive
            ref_seq = "".join(ref_index[loci.contig][ref_start:ref_stop])
            if ref_seq != variant.REF:
                raise ReferenceError(
                    f"VCF REF {variant.REF} does not match the expected reference "
                    f"sequence {ref_seq} at loci {variant.CHROM} position {variant.POS}"
                )

        variant.INFO[LOCI_ID] = variant.CHROM
        variant.INFO[LOCI_POS_ID] = variant.POS
        normalised_pos = loci.start + variant.POS
        variant.set_pos(normalised_pos - 1)  # this function takes 0-based
        var_fields = str(variant).rstrip().split("\t")
        var_fields[0] = loci.contig
        out_variant = "\t".join(var_fields)
        print(out_variant, file=out_vcf)

    logging.info("Positions normalised!")
    vcf_reader.close()


if __name__ == "__main__":
    main()
