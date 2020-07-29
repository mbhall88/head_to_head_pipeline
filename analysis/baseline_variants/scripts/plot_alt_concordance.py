import json
from collections import defaultdict
from pathlib import Path
from typing import List

import pandas as pd
from bokeh.models import ColumnDataSource
from bokeh.models import Legend
from bokeh.palettes import Set2
from bokeh.plotting import figure, output_file
from bokeh.plotting import save

TOOLS = "pan,wheel_zoom,box_zoom,reset,box_select,lasso_select,undo,redo,save,hover"
PIXEL_INCHES = 96
HEIGHT = 8 * PIXEL_INCHES
WIDTH = 13 * PIXEL_INCHES

JSON_FILES: List[Path] = [Path(p) for p in snakemake.input.jsons]
OUTFILE = snakemake.output.plot
TITLE = snakemake.params.title
Y_VAR = snakemake.params.y_var
X_VAR = snakemake.params.x_var
COLOUR_BY = "site"
INDEX = "sample"

data = defaultdict(dict)
for path in JSON_FILES:
    site = path.parts[-2]
    sample = path.name.split(".")[0]
    with path.open() as fp:
        data[site, sample] = json.load(fp)

df = pd.DataFrame(data).T
df.reset_index(inplace=True)
df.rename(columns={"level_0": COLOUR_BY, "level_1": INDEX}, inplace=True)
df.set_index(INDEX, drop=True, inplace=True)

# float precision for tooltips can be found at https://docs.bokeh.org/en/latest/docs/reference/models/formatters.html#bokeh.models.formatters.NumeralTickFormatter.format
float_fmt = "0.0[000000]"
tooltips = [
    (INDEX, f"@{INDEX}"),
    (COLOUR_BY, f"@{COLOUR_BY}"),
    (X_VAR.capitalize(), f"@{X_VAR}{{({float_fmt})}}"),
    (Y_VAR.capitalize().replace("_", " "), f"@{Y_VAR}{{({float_fmt})}}"),
]
legend_var = COLOUR_BY
categories = set(df[COLOUR_BY])
palette = Set2[len(categories)]

# inline effectively allows the plot to work offline
output_file(OUTFILE, title=TITLE, mode="inline")

plot = figure(
    tools=TOOLS,
    height=HEIGHT,
    width=WIDTH,
    tooltips=tooltips,
    active_drag="box_zoom",
    active_scroll="wheel_zoom",
    title=TITLE,
    toolbar_location="below",
)
plot.xaxis.axis_label = X_VAR.capitalize()
plot.yaxis.axis_label = Y_VAR.replace("_", " ").capitalize()
legend = Legend(
    click_policy="hide",
    title=legend_var.capitalize(),
    background_fill_alpha=0.1,
    title_text_font_style="bold",
)
plot.add_layout(legend, "right")

for i, cat in enumerate(categories):
    cat_data = df[df[COLOUR_BY] == cat]
    source = ColumnDataSource(cat_data)
    plot.circle(
        x=X_VAR,
        y=Y_VAR,
        source=source,
        color=palette[i],
        legend_label=cat,
        size=10,
        alpha=0.4,
    )

plot.ygrid.minor_grid_line_color = "black"
plot.ygrid.minor_grid_line_alpha = 0.05
plot.xgrid.minor_grid_line_color = "black"
plot.xgrid.minor_grid_line_alpha = 0.05

save(plot)
