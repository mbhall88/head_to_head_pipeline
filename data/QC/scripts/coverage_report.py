import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List

import click
import pandas as pd
from bokeh.models import ColumnDataSource
from bokeh.models import Legend
from bokeh.palettes import Set2
from bokeh.plotting import figure, save, output_file
from bokeh.transform import factor_cmap, factor_mark

MARKERS = ["circle", "triangle", "square"]
TOOLS = "pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,undo,redo,save,hover"
PIXEL_INCHES = 96
HEIGHT = 8 * PIXEL_INCHES
WIDTH = 13 * PIXEL_INCHES


class RipgrepError(Exception):
    pass


def ripgrep_extract_covg(file: Path) -> float:
    pattern = r"Actual cov.*\s(?P<covg>\d+?\.?\d+)x"
    extra_params = ["--replace", "$covg", "--only-matching", "--no-line-number"]
    process = subprocess.Popen(
        ["rg", *extra_params, pattern, str(file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    exit_code = process.wait()
    if exit_code != 0:
        raise RipgrepError(
            f"Failed to execute the rg command on the file {file} due to the "
            f"following error:\n{process.stderr.read()}"
        )
    return float(process.stdout.read().strip())


@click.command()
@click.argument(
    "log-dirs", type=click.Path(exists=True, file_okay=False), nargs=-1,
)
@click.help_option("--help", "-h")
@click.option(
    "-a",
    "--assignment-dir",
    help="Directory containing lineage assignment CSV files.",
    type=click.Path(exists=True, file_okay=False),
    required=True,
)
@click.option(
    "-o",
    "--outfile",
    type=click.Path(exists=False, dir_okay=False),
    default="-",
    show_default=True,
    help="The filepath to write the output HTML file to.",
)
def main(
    assignment_dir: str, log_dirs: List[str], outfile: str,
):
    """This script generates a HTML file containing a plot of the coverage for each
    sample after they have been through the QC pipeline.\n
    LOG_DIRS: Director(y/ies) containing the subsampling log files. (Coverage is
    extracted from these)
    ),
    """
    logfiles = []
    for d in log_dirs:
        logfiles.extend(Path(d).rglob("*.log"))

    data = defaultdict(dict)
    for file in logfiles:
        tech = "illumina" if "illumina" in file.parts[2] else "nanopore"
        sample = file.with_suffix("").name
        site = file.parts[-2]
        covg = ripgrep_extract_covg(file)
        data[sample].update({f"{tech}_covg": covg, "site": site})

    assignment_files = Path(assignment_dir).rglob("*.csv")
    for file in assignment_files:
        fields = file.read_text().split("\n")[1].split(",")
        sample = file.name.split(".")[0]
        data[sample]["lineage"] = fields[1]

    df = pd.DataFrame(data).T
    df.index.name = "sample"
    cds = ColumnDataSource(df)

    sites = list(set(df["site"]))
    lineages = list(set(df["lineage"]))
    tooltips = [
        ("index", "@sample"),
        ("Illumina", "@illumina_covg"),
        ("Nanopore", "@nanopore_covg"),
        ("site", "@site"),
        ("lineage", "@lineage"),
    ]
    title = "Sample coverage for different technologies after quality control"
    palette = Set2[len(lineages)]
    legend_var = "lineage"

    # inline effectively allows the plot to work offline
    output_file(outfile, title=title, mode="inline")

    p = figure(
        tools=TOOLS,
        height=HEIGHT,
        width=WIDTH,
        tooltips=tooltips,
        active_drag="box_zoom",
        active_scroll="wheel_zoom",
        title=title,
    )
    p.yaxis.axis_label = "Nanopore Coverage"
    p.xaxis.axis_label = "Illumina Coverage"
    legend = Legend(
        click_policy="hide",
        location="top_left",
        title=legend_var.capitalize(),
        background_fill_alpha=0.1,
    )
    p.add_layout(legend)

    p.scatter(
        source=cds,
        y="nanopore_covg",
        x="illumina_covg",
        size=10,
        fill_alpha=0.4,
        marker=factor_mark("site", MARKERS, sites),
        color=factor_cmap("lineage", palette, lineages),
        legend_field=legend_var,
    )

    save(p)


if __name__ == "__main__":
    main()
