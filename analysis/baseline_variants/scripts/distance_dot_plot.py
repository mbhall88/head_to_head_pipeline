from pathlib import Path

import click
import numpy as np
from scipy import stats
import pandas as pd
from bokeh.io import output_file, save
from bokeh.plotting import figure

TOOLS = "pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,undo,redo,save,hover"
PIXEL_INCHES = 96
HEIGHT = 8 * PIXEL_INCHES
WIDTH = 13 * PIXEL_INCHES
TITLE = "Distance matrix dot plot"
PAIR_IDX = ("sample1", "sample2")


def load_matrix(fpath: str, delim: str) -> pd.DataFrame:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        names = header.split(delim)[1:]
        for row in map(str.rstrip, instream):
            matrix.append(map(int, row.split(delim)[1:]))
    df = pd.DataFrame(matrix, index=names, columns=names)
    # remove the lower triangle of the matrix and the middle diagonal
    return df.where(np.triu(np.ones(df.shape), k=1).astype(np.bool))


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-x",
    "--x-matrix",
    help="Distance matrix to plot on the x-axis.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-y",
    "--y-matrix",
    help="Distance matrix to plot on the y-axis.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
)
@click.option(
    "-X",
    "--xname",
    help=(
        "Name for the x-axis matrix. If not given, the stem of the filename will be "
        "used."
    ),
)
@click.option(
    "-Y",
    "--yname",
    help=(
        "Name for the y-axis matrix. If not given, the stem of the filename will be "
        "used."
    ),
)
@click.option(
    "-o",
    "--output",
    help="Path to save HTML plot to.",
    type=click.Path(dir_okay=False, writable=True),
    default="dotplot.html",
    show_default=True,
)
@click.option(
    "-d",
    "--delim",
    help="Delimiter used in the matrix. [ default: '\\t' ]",
    default="\t",
)
@click.option(
    "-t", "--title", help="Title for the plot.", default=TITLE, show_default=True
)
@click.option("--width", help="Plot width in pixels", default=WIDTH, show_default=True)
@click.option(
    "--height", help="Plot height in pixels", default=HEIGHT, show_default=True
)
@click.option(
    "-T",
    "--threshold",
    help="Only plot pairs where SNP distance is no greater than this threshold.",
    type=int,
)
def main(
    x_matrix: str,
    y_matrix: str,
    xname: str,
    yname: str,
    output: str,
    delim: str,
    title: str,
    width: int,
    height: int,
    threshold: int,
):
    if not xname:
        xname = Path(x_matrix).name.split(".")[0]
    if not yname:
        yname = Path(y_matrix).name.split(".")[0]

    xmat = load_matrix(x_matrix, delim)
    xdf = xmat.stack().rename(xname)
    xdf = xdf.rename_axis(PAIR_IDX)

    ymat = load_matrix(y_matrix, delim)
    ydf = ymat.stack().rename(yname)
    ydf = ydf.rename_axis(PAIR_IDX)

    # merge the matrices
    df = pd.concat([xdf, ydf], axis=1)
    df = df.reset_index().rename(
        columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
    )
    if threshold:
        df = df.query(f"{xname} <= @threshold")

    dot_alpha = 0.25
    line_alpha = 0.75
    line_width = 2
    tooltips = [
        ("pair", f"@{PAIR_IDX[0]} x @{PAIR_IDX[1]}"),
        (f"{xname} dist.", f"@{xname}"),
        (f"{yname} dist.", f"@{yname}"),
    ]

    output_file(output, title=title, mode="inline")

    fig = figure(
        tools=TOOLS,
        height=height,
        width=width,
        tooltips=tooltips,
        active_drag="box_zoom",
        active_scroll="wheel_zoom",
        title=title,
    )
    fig.circle(x=xname, y=yname, alpha=dot_alpha, source=df)

    # determine best fit line
    x = df[xname]
    y = df[yname]
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    y_predicted = [gradient * i + intercept for i in x]

    # add line of best fit to the data
    fitted_xcoord = [min(x), max(x)]
    fitted_ycoord = [min(y_predicted), max(y_predicted)]
    fitted_equation = f"y={gradient:.2f}x + {intercept:.2f} (r={r_value:.3f})"
    fig.line(
        x=fitted_xcoord,
        y=fitted_ycoord,
        legend_label=fitted_equation,
        line_color="#4c566a",
        line_dash="dashed",
        line_width=line_width,
    )

    # add "expected" line if both caller have same distance for each pair
    expected_xcoord = [0, max(x)]
    expected_ycoord = expected_xcoord
    expected_equation = "y=x (r=1.0)"
    fig.line(
        expected_xcoord,
        expected_ycoord,
        line_color="red",
        line_dash="dashed",
        line_alpha=line_alpha,
        legend_label=expected_equation,
        line_width=line_width,
    )

    fig.legend.location = "top_left"
    fig.legend.label_text_font_size = "12pt"
    fig.xaxis.axis_label = f"{xname} SNP distance"
    fig.yaxis.axis_label = f"{yname} SNP distance"
    save(fig)


if __name__ == "__main__":
    main()
