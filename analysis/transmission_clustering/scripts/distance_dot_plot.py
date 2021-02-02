from pathlib import Path
from typing import List, Dict, Tuple

import click
import numpy as np
import pandas as pd
from bokeh.io import output_file, save
from bokeh.models import ColumnDataSource, Legend
from bokeh.palettes import Set2
from bokeh.plotting import figure
from scipy import stats

TOOLS = "pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,undo,redo,save,hover"
PIXEL_INCHES = 96
HEIGHT = 8 * PIXEL_INCHES
WIDTH = 13 * PIXEL_INCHES
TITLE = "Distance matrix dot plot"
PAIR_IDX = ("sample1", "sample2")


def load_matrix(fpath: str, delim: str, name: str) -> pd.DataFrame:
    matrix = []
    with open(fpath) as instream:
        header = next(instream).rstrip()
        names = header.split(delim)[1:]
        for row in map(str.rstrip, instream):
            matrix.append(map(int, row.split(delim)[1:]))
    df = pd.DataFrame(matrix, index=names, columns=names)
    # remove the lower triangle of the matrix and the middle diagonal
    df = df.where(np.triu(np.ones(df.shape), k=1).astype(np.bool)).sort_index()
    df = df.stack().rename(name).astype(int)
    df = df.rename_axis(PAIR_IDX)
    sorted_idx = [sorted(ix) for ix in df.index]
    df.index = pd.MultiIndex.from_tuples(sorted_idx, names=df.index.names)
    return df


class PlotFactory:
    def __init__(
        self,
        xname: str,
        xcol: str,
        yval_name: str,
        colour_by: str,
        palette: Dict[int, List[str]],
        data: pd.DataFrame,
        tools: str = TOOLS,
        height: int = HEIGHT,
        width: int = WIDTH,
        float_fmt: str = "0.0[000000]",
        toolbar_location: str = "below",
        legend_location: str = "right",
        point_size: int = 10,
        point_alpha: float = 0.4,
        minor_grid_colour: str = "black",
        minor_grid_alpha: float = 0.1,
        x_axis_type: str = "auto",
        y_axis_type: str = "auto",
        line_width: float = 2,
        line_alpha: float = 0.8,
    ):
        # float precision for tooltips can be found at https://docs.bokeh.org/en/latest/docs/reference/models/formatters.html#bokeh.models.formatters.NumeralTickFormatter.format
        self.float_fmt = float_fmt
        self.xname = xname
        self.xcol = xcol
        self.yval_name = yval_name
        self.colour_by = colour_by
        self.data = data
        self.categories = set(self.data[self.colour_by])
        self.palette = palette[max(3, len(self.categories))]
        self.tools = tools
        self.height = height
        self.width = width
        self.toolbar_location = toolbar_location
        self.legend_location = legend_location
        self.point_size = point_size
        self.point_alpha = point_alpha
        self.line_alpha = line_alpha
        self.line_width = line_width
        self.minor_grid_colour = minor_grid_colour
        self.minor_grid_alpha = minor_grid_alpha
        self.x_axis_type = x_axis_type
        self.y_axis_type = y_axis_type

    def _build_tooltips(self) -> List[Tuple[str, str]]:
        return [
            ("pair", f"@{PAIR_IDX[0]} x @{PAIR_IDX[1]}"),
            (f"{self.xname} dist.", f"@{self.xcol}"),
            ("Nanopore dist.", f"@{self.yval_name}"),
            ("Nanopore caller", f"@{self.colour_by}"),
        ]

    @property
    def legend_var(self) -> str:
        return self.colour_by

    def _create_figure(self, title: str) -> figure:
        tooltips = self._build_tooltips()
        return figure(
            tools=self.tools,
            height=self.height,
            width=self.width,
            tooltips=tooltips,
            active_drag="box_zoom",
            active_scroll="wheel_zoom",
            title=title,
            toolbar_location=self.toolbar_location,
            x_axis_type=self.x_axis_type,
            y_axis_type=self.y_axis_type,
        )

    def _create_legend(self) -> Legend:
        return Legend(
            click_policy="hide",
            title=" ".join(self.colour_by.capitalize().split("_")),
            background_fill_alpha=0.1,
            title_text_font_style="bold",
        )

    def generate_plot(
        self, outfile: str, x_var: str, y_var: str, title: str, xlabel: str, ylabel: str
    ):
        fig = self._create_figure(title)
        legend = self._create_legend()
        fig.xaxis.axis_label = xlabel
        fig.yaxis.axis_label = ylabel
        fig.add_layout(legend, self.legend_location)

        # inline effectively allows the plot to work offline
        output_file(outfile, title=title, mode="inline")

        for i, cat in enumerate(self.categories):
            cat_data = self.data[self.data[self.colour_by] == cat]
            source = ColumnDataSource(cat_data)
            fig.circle(
                x=x_var,
                y=y_var,
                source=source,
                color=self.palette[i],
                legend_label=cat,
                size=self.point_size,
                alpha=self.point_alpha,
            )

            # determine best fit line
            x = cat_data[x_var]
            y = cat_data[y_var]
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
                line_color=self.palette[i],
                line_dash="dashed",
                line_width=self.line_width,
            )

        # add "expected" line if both caller have same distance for each pair
        expected_xcoord = [0, max(x)]
        expected_ycoord = expected_xcoord
        expected_equation = "y=x (r=1.0)"
        fig.line(
            expected_xcoord,
            expected_ycoord,
            line_color="black",
            line_dash="dashed",
            line_alpha=self.line_alpha,
            legend_label=expected_equation,
            line_width=self.line_width,
        )

        fig.ygrid.minor_grid_line_color = self.minor_grid_colour
        fig.ygrid.minor_grid_line_alpha = self.minor_grid_alpha
        fig.xgrid.minor_grid_line_color = self.minor_grid_colour
        fig.xgrid.minor_grid_line_alpha = self.minor_grid_alpha

        save(fig)


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
    "--y-matrices",
    help="Distance matrix to plot on the y-axis.",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    multiple=True,
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
    "--ynames",
    help=(
        "Name for the y-axis matrix. If not given, the stem of the filename(s) will be "
        "used. For multiple names, use a comma-separated list in the same order as -Y"
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
    y_matrices: List[str],
    xname: str,
    ynames: str,
    output: str,
    delim: str,
    title: str,
    width: int,
    height: int,
    threshold: int,
):
    if not xname:
        xname = Path(x_matrix).name.split(".")[0]
    if not ynames:
        ynames = [p.name.split(".")[0] for p in map(Path, y_matrices)]
    else:
        ynames = ynames.split(",")

    xcol = f"{xname}_dist"
    xdf = load_matrix(x_matrix, delim, xcol)

    ydfs = []
    for yname, y_matrix in zip(ynames, y_matrices):
        ydf = load_matrix(y_matrix, delim, yname)
        ydfs.append(ydf)

    # merge the matrices
    yvar_name = "nanopore_caller"
    yval_name = "nanopore_dist"
    df = pd.concat([xdf, *ydfs], axis=1)
    df = df.reset_index().rename(
        columns={"level_0": PAIR_IDX[0], "level_1": PAIR_IDX[1]}
    )
    df = df.melt(
        id_vars=[*PAIR_IDX, xcol],
        value_vars=ynames,
        var_name=yvar_name,
        value_name=yval_name,
    )

    if threshold:
        df = df.query(f"{xcol} <= @threshold")

    plotter = PlotFactory(
        xname=xname,
        xcol=xcol,
        yval_name=yval_name,
        colour_by=yvar_name,
        palette=Set2,
        data=df,
        point_alpha=0.4,
        width=width,
        height=height,
    )

    plotter.generate_plot(
        outfile=output,
        x_var=xcol,
        y_var=yval_name,
        title=title,
        xlabel="COMPASS (Illumina) distance (bp)",
        ylabel="Nanopore distance (bp)",
    )


if __name__ == "__main__":
    main()
