import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, output_file, save

TOOLS = "pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,undo,redo,save,hover"
PIXEL_INCHES = 96
HEIGHT = 8 * PIXEL_INCHES
WIDTH = 13 * PIXEL_INCHES
json_dir = Path(sys.argv[1])
JSON_FILES: List[Path] = list(json_dir.rglob("*.json"))
outdir = Path(sys.argv[2])
outdir.mkdir(exist_ok=True)
INDEX: str = "sample"


class PlotFactory:
    def __init__(
        self,
        index: str,
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
    ):
        # float precision for tooltips can be found at https://docs.bokeh.org/en/latest/docs/reference/models/formatters.html#bokeh.models.formatters.NumeralTickFormatter.format
        self.float_fmt = float_fmt
        self.index = index
        self.data = data
        self.tools = tools
        self.height = height
        self.width = width
        self.toolbar_location = toolbar_location
        self.legend_location = legend_location
        self.point_size = point_size
        self.point_alpha = point_alpha
        self.minor_grid_colour = minor_grid_colour
        self.minor_grid_alpha = minor_grid_alpha
        self.x_axis_type = x_axis_type
        self.y_axis_type = y_axis_type

    def _build_tooltips(
        self, x_var: str, y_var: str, xlabel: str, ylabel: str
    ) -> List[Tuple[str, str]]:
        return [
            (self.index, f"@{self.index}"),
            (xlabel, f"@{x_var}{{({self.float_fmt})}}"),
            (ylabel, f"@{y_var}{{({self.float_fmt})}}"),
        ]

    def _create_figure(
        self, x_var: str, y_var, title: str, xlabel: str, ylabel: str
    ) -> figure:
        tooltips = self._build_tooltips(x_var, y_var, xlabel, ylabel)
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

    def generate_plot(
        self, outfile: str, x_var: str, y_var: str, title: str, xlabel: str, ylabel: str
    ):
        fig = self._create_figure(x_var, y_var, title, xlabel, ylabel)
        fig.xaxis.axis_label = xlabel
        fig.yaxis.axis_label = ylabel

        # inline effectively allows the plot to work offline
        output_file(outfile, title=title, mode="inline")

        source = ColumnDataSource(self.data)
        fig.circle(
            x=x_var,
            y=y_var,
            source=source,
            size=self.point_size,
            alpha=self.point_alpha,
        )

        fig.ygrid.minor_grid_line_color = self.minor_grid_colour
        fig.ygrid.minor_grid_line_alpha = self.minor_grid_alpha
        fig.xgrid.minor_grid_line_color = self.minor_grid_colour
        fig.xgrid.minor_grid_line_alpha = self.minor_grid_alpha

        save(fig)


def load_concordance_data(json_files: List[Path]) -> pd.DataFrame:
    data = defaultdict(dict)
    for path in json_files:
        sample = path.name.split(".")[0]
        with path.open() as fp:
            data[sample] = json.load(fp)

    df = pd.DataFrame(data).T
    df.reset_index(inplace=True)
    return df


concordance_df = load_concordance_data(JSON_FILES)
concordance_df.rename(columns={"index": INDEX}, inplace=True)
concordance_df.set_index(INDEX, drop=True, inplace=True)


plotter = PlotFactory(index=INDEX, data=concordance_df)

plotter.generate_plot(
    outfile=outdir / "alt_concordance.html",
    x_var="concordance",
    y_var="call_rate",
    title="Call rate vs. Concordance at ALT positions",
    xlabel="Concordance",
    ylabel="Call rate",
)

plotter.generate_plot(
    outfile=outdir / "gw_concordance.html",
    x_var="gw_concordance",
    y_var="gw_call_rate",
    title="Genome-wide Call rate vs. Concordance",
    xlabel="Concordance",
    ylabel="Call rate",
)
