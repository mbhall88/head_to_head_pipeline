from typing import TextIO

import click
import pysam


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--infile",
    help="Assembly to chop up.",
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
    "--chunk-size",
    help="The length of the pieces to chop the assembly into.",
    type=int,
    default=10_000,
    show_default=True,
)
@click.option(
    "--min-tail-size",
    help=(
        "When splitting the assembly into chunks, only keep the last chunk if it is "
        "longer than (or equal to) this value."
    ),
    default=500,
    show_default=True,
)
def main(infile: str, output: TextIO, chunk_size: int, min_tail_size: int):
    """Takes an assembly as input and produces an output of 'reads': the
    assembly chopped into pieces (chunks).
    """
    with pysam.FastxFile(infile) as contigs:
        for contig in contigs:
            chunk_num = 0
            seq = contig.sequence
            for i in range(0, len(seq), chunk_size):
                chunk = seq[i : i + chunk_size]
                if len(chunk) >= min_tail_size:
                    chunk_num += 1
                    chunk_name = f"{contig.name}_chunk={chunk_num}"
                    print(f">{chunk_name}", file=output)
                    print(chunk, file=output)


if __name__ == "__main__":
    main()
