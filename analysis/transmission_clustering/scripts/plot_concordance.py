import json
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Tuple

import pandas as pd
from bokeh.models import ColumnDataSource, Legend
from bokeh.palettes import Set2
from bokeh.plotting import figure, output_file, save

TOOLS = "pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,undo,redo,save,hover"
PIXEL_INCHES = 96
HEIGHT = 8 * PIXEL_INCHES
WIDTH = 13 * PIXEL_INCHES
JSON_FILES: List[Path] = [Path(p) for p in snakemake.input.jsons]
COLOUR_BY: str = snakemake.params.colour_by
INDEX: str = snakemake.params.index
LOG_SCALE: bool = snakemake.params.log_scale


class PlotFactory:
    def __init__(
        self,
        index: str,
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
    ):
        # float precision for tooltips can be found at https://docs.bokeh.org/en/latest/docs/reference/models/formatters.html#bokeh.models.formatters.NumeralTickFormatter.format
        self.float_fmt = float_fmt
        self.index = index
        self.colour_by = colour_by
        self.data = data
        _, cats = zip(*self.data.index)
        self.categories = set(cats)
        self.palette = palette[max(len(self.categories), 3)]
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
            (self.colour_by, f"@{self.colour_by}"),
            (xlabel, f"@{x_var}{{({self.float_fmt})}}"),
            (ylabel, f"@{y_var}{{({self.float_fmt})}}"),
        ]

    @property
    def legend_var(self) -> str:
        return self.colour_by

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

    def _create_legend(self) -> Legend:
        return Legend(
            click_policy="hide",
            title=self.colour_by.capitalize(),
            background_fill_alpha=0.1,
            title_text_font_style="bold",
        )

    def generate_plot(
        self, outfile: str, x_var: str, y_var: str, title: str, xlabel: str, ylabel: str
    ):
        fig = self._create_figure(x_var, y_var, title, xlabel, ylabel)
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

        fig.ygrid.minor_grid_line_color = self.minor_grid_colour
        fig.ygrid.minor_grid_line_alpha = self.minor_grid_alpha
        fig.xgrid.minor_grid_line_color = self.minor_grid_colour
        fig.xgrid.minor_grid_line_alpha = self.minor_grid_alpha

        save(fig)


def load_concordance_data(json_files: List[Path]) -> pd.DataFrame:
    data = defaultdict(dict)
    for path in json_files:
        site = path.parts[-2]
        sample = path.name.split(".")[0]
        if "pandora" in str(path):
            caller = "pandora"
        elif "baseline" in str(path):
            caller = "bcftools"
        else:
            raise NotImplementedError(f"{path} is not from a known caller")
        with path.open() as fp:
            data[site, sample, caller] = json.load(fp)

    df = pd.DataFrame(data).T
    df.reset_index(inplace=True)
    return df


concordance_df = load_concordance_data(JSON_FILES)
concordance_df.rename(
    columns={"level_0": "site", "level_1": INDEX, "level_2": COLOUR_BY}, inplace=True
)
concordance_df.set_index([INDEX, COLOUR_BY], drop=True, inplace=True)
covg_df = pd.read_csv(snakemake.input.coveragesheet, index_col=INDEX)
for sample in covg_df.index:
    if sample not in concordance_df.index:
        continue
    depth = covg_df.at[sample, "nanopore_covg"]
    concordance_df.at[sample, "depth"] = depth

xscale = "log" if LOG_SCALE else "auto"
yscale = "log" if LOG_SCALE else "auto"
plotter = PlotFactory(
    index=INDEX,
    colour_by=COLOUR_BY,
    palette=Set2,
    data=concordance_df,
    x_axis_type=xscale,
    y_axis_type=yscale,
)

plotter.generate_plot(
    outfile=snakemake.output.alt_plot,
    x_var="concordance",
    y_var="call_rate",
    title="Call rate vs. Concordance at ALT positions",
    xlabel="Concordance",
    ylabel="Call rate",
)

plotter.generate_plot(
    outfile=snakemake.output.gw_plot,
    x_var="gw_concordance",
    y_var="gw_call_rate",
    title="Genome-wide Call rate vs. Concordance",
    xlabel="Concordance",
    ylabel="Call rate",
)

if LOG_SCALE:
    # don't want depth to be log-scaled
    plotter.x_axis_type = "auto"

plotter.generate_plot(
    outfile=snakemake.output.depth_call_rate_plot,
    x_var="depth",
    y_var="call_rate",
    title="Effect of depth on call rate",
    xlabel="Median Depth",
    ylabel="Call rate",
)

plotter.generate_plot(
    outfile=snakemake.output.depth_gw_call_rate_plot,
    x_var="depth",
    y_var="gw_call_rate",
    title="Effect of depth on genome-wide call rate",
    xlabel="Median Depth",
    ylabel="Call rate",
)

plotter.generate_plot(
    outfile=snakemake.output.depth_concordance_plot,
    x_var="depth",
    y_var="concordance",
    title="Effect of depth on concordance",
    xlabel="Median Depth",
    ylabel="Concordance",
)

plotter.generate_plot(
    outfile=snakemake.output.depth_gw_concordance_plot,
    x_var="depth",
    y_var="gw_concordance",
    title="Effect of depth on genome-wide concordance",
    xlabel="Median Depth",
    ylabel="Concordance",
)
