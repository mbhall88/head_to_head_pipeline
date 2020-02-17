import click
import logging
from pathlib import Path
import pysam


def validate_file(ctx, param, value):
    if value.exists():
        return value
    else:
        raise click.BadParameter(f"{value} does not exist!")


def count_bases(path: Path) -> int:
    with pysam.FastxFile(path) as fastx:
        return sum(len(record.sequence) for record in fastx)


def get_region_length(line: str) -> int:
    fields = line.split()
    start = int(fields[1])
    end = int(fields[2])
    return end - start


def test():
    line = "2       0       128"

    actual = get_region_length(line)
    expected = 128

    assert actual == expected

    line = "2       3655    3656"
    actual = get_region_length(line)
    expected = 1

    assert actual == expected


@click.command()
@click.option(
    "--bed", help="Bed file.", type=Path, required=True, callback=validate_file
)
@click.option(
    "--assembly",
    help="Fasta of assembly bed file is related to.",
    type=Path,
    required=True,
    callback=validate_file,
)
def main(bed: Path, assembly: Path):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

    num_bases = count_bases(assembly)
    logging.info(f"{num_bases} bases found in assembly.")

    with bed.open() as fh:
        num_masked_positions = sum(get_region_length(line) for line in fh)

    logging.info(f"There are {num_masked_positions} positions in the BED file.")

    perc_masked = round(num_masked_positions / num_bases * 100, 5)

    print(f"{perc_masked}% of the assembly is masked")


if __name__ == "__main__":
    test()
    main()
