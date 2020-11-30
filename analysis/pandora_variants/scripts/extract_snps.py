import logging

import click
from cyvcf2 import VCF, Writer


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--in-vcf",
    help="Pandora VCF file to extract SNPs from.",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
    default="-",
    show_default=True,
)
@click.option(
    "-o",
    "--out-vcf",
    help="Path to write output VCF to. Note: can only output VCF currently",
    type=click.Path(dir_okay=False, allow_dash=True, writable=True),
    default="-",
    show_default=True,
)
@click.option("-M", "--keep-mnps", help="Also keep MNPs", is_flag=True)
@click.option("-v", "--verbose", help="Turns on debug-level logging.", is_flag=True)
def main(
    in_vcf: str, out_vcf: str, keep_mnps: bool, verbose: bool,
):
    """Extract SNPs from a pandora VCF. It keeps records regardless of the called
    allele, provided they record is a SNP. If the ALT allele is a '.' and the REF is a
    single character, it is considered a SNP.
    It is assumed that each record has no more than 1 alternate allele.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )
    vcf_reader = VCF(in_vcf)
    vcf_writer = Writer(out_vcf, tmpl=vcf_reader)

    logging.info("Checking for records to keep...")

    for record in vcf_reader:
        keep_this_record = False
        ref_len = len(record.REF)
        ref_is_single_base = ref_len == 1
        empty_alt = not bool(record.ALT)
        if empty_alt and (keep_mnps or ref_is_single_base):
            keep_this_record = True
        elif ref_len == len(record.ALT[0]) and (ref_is_single_base or keep_mnps):
            keep_this_record = True

        if keep_this_record:
            vcf_writer.write_record(record)
        else:
            logging.debug(
                f"Discarding record CHROM: {record.CHROM} at POS: {record.POS}"
            )

    vcf_writer.close()
    vcf_reader.close()
    logging.info("Done!")


if __name__ == "__main__":
    main()
