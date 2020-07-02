import click
import pandas as pd


class RequiredIf(click.Option):
    def __init__(self, *args, **kwargs):
        self.required_if = kwargs.pop("required_if")
        assert self.required_if, "'required_if' parameter required"
        kwargs["help"] = (
            kwargs.get("help", "")
            + " NOTE: This argument is mutually inclusive with %s" % self.required_if
        ).strip()
        super(RequiredIf, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        we_are_present = self.name in opts
        other_present = self.required_if in opts

        if not other_present:
            if we_are_present:
                raise click.UsageError(
                    "Illegal usage: `%s` is mutually inclusive with `%s`"
                    % (self.name, self.required_if)
                )

        return super(RequiredIf, self).handle_parse_result(ctx, opts, args)


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-i",
    "--samfile",
    help="{B,CR,S}AM file of reads mapped to decontamination database.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-m",
    "--metadata",
    help=(
        "TSV file containing information about each reference in the database. Column "
        "1 is the category, column 2 is whether the reference is contamination, column "
        "3 is the accession ID for the reference."
    ),
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-1",
    "--read1",
    help="The original fastq file the reads come from.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-2",
    "--read2",
    help="The mate for --read1. Only relevant if using paired reads.",
    type=click.Path(exists=True, dir_okay=False),
    cls=RequiredIf,
    required_if="outfile2",
)
@click.option(
    "-o1",
    "--outfile1",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Path for the output fastq.",
)
@click.option(
    "-o2",
    "--outfile2",
    type=click.Path(dir_okay=False, writable=True),
    help="Path for the mate of --outfile1 fastq. Only relevant if using paired reads.",
    cls=RequiredIf,
    required_if="read2",
)
@click.option(
    "--ignore-secondary/--include-secondary",
    help="Ignore organism assignments for secondary alignments?",
    default=True,
    show_default=True,
)
def main(
    samfile: str,
    metadata: str,
    read1: str,
    read2: str,
    outfile1: str,
    outfile2: str,
    ignore_secondary: bool,
):
    """This script generates the text file input required to make a krona plot.
    """
    metadata_df = pd.read_table(
        metadata,
        header=None,
        names=["organism", "contamination", "accession"],
        index_col="accession",
    )
    #todo - https://github.com/mbhall88/head_to_head_pipeline/blob/576809b8d67f663e38aa72928bea9d93dd22569e/analysis/scripts/filter.py
    #todo - https://github.com/iqbal-lab-org/clockwork/blob/master/python/clockwork/contam_remover.py

if __name__ == "__main__":
    main()
