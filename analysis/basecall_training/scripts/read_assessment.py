from typing import TextIO, Iterator

import click
from pafpy import PafFile, PafRecord


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--infile",
    help="PAF file to assess.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.File(mode="w"),
    help="The path to write the output file to. Use '-' for stdout",
)
@click.option(
    "--delim",
    help="The column delimiter to use in the output file.",
    default=",",
    show_default=True,
)
@click.option("--primary-only", help="Only assess primary alignments.", is_flag=True)
@click.option(
    "--min-cov",
    help=(
            "Minimum read coverage required for a record to be assessed. Read coverage "
            "is defined as the proportion of the query sequence involved in the "
            "alignment. It is the aligned length minus the read length."
    ),
    type=float,
    default=0,
    show_default=True,
)
def main(
        infile: str, output: TextIO, primary_only: bool, min_cov: float, delim: str,
):
    """A script to produce data relevant to assessing the per-read accuracy for a PAF
    file. The output is a file with columns containing read identifier, read
    length, BLAST identity, and relative length of read compared to the reference."""

    def is_valid(record: PafRecord) -> bool:
        aln_type_is_valid = True
        if primary_only:
            aln_type_is_valid = record.is_primary()

        return aln_type_is_valid and record.query_coverage >= min_cov

    header = delim.join(["read_id", "read_len", "identity", "relative_len"])
    print(header, file=output)

    with PafFile(infile) as paf:
        valid_records: Iterator[PafRecord] = filter(is_valid, paf)
        for record in valid_records:
            fields = [
                record.qname,  # read_id
                record.qlen,  # read length
                record.blast_identity(),
                record.relative_length,
            ]
            line = delim.join(map(str, fields))
            print(line, file=output)


if __name__ == "__main__":
    main()
