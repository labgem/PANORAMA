#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from pathlib import Path
from typing import List

# installed libraries
import pandas as pd
from bokeh.io import output_file, save, export_png
from bokeh.plotting import figure
from bokeh.transform import transform, linear_cmap
from bokeh.palettes import Colorblind, Reds256, linear_palette
from bokeh.models import ColumnDataSource, LinearColorMapper, BasicTicker

# local libraries
from panorama.pangenomes import Pangenome


# def associate_systems_to_modules():

def get_association_df(pangenome: Pangenome, association: List[str]) -> pd.DataFrame:
    columns = ['system number', 'system name', 'families'] + ['modules'] if 'modules' in association else []
    association_list = {}
    for system in pangenome.systems:
        if 'modules' in association:
            association_list[system.ID] = [system.name, ",".join([fam.name for fam in system.families]),
                                           ",".join(map(str, system.modules))]

    association_df = pd.DataFrame.from_dict(association_list, orient='index', columns=columns[1:])
    association_df.index.name = columns[0]
    return association_df


def write_correlation_matrix(association_df: pd.DataFrame, output: Path, name: str, out_format: List[str] = None):
    out_format = out_format if out_format is not None else ['html']

    association_split = association_df.drop(columns=['families']).join(
        association_df['modules'].str.get_dummies(sep=',')).drop(columns=['modules'])

    correlation_matrix = association_split.groupby('system name').sum()
    correlation_matrix.sort_index(key=lambda x: x.str.lower(), ascending=False, inplace=True)
    correlation_matrix = correlation_matrix[sorted(correlation_matrix.columns.tolist(), key=lambda x: int(x))]
    source = ColumnDataSource(data=dict(
        modules=list(correlation_matrix.columns) * len(correlation_matrix.index),
        systems=list(correlation_matrix.index) * len(correlation_matrix.columns),
        corr=correlation_matrix.to_numpy().flatten(),
    ))

    high_corr = correlation_matrix.values.max()
    if high_corr == 1:
        color_palette = ["#ffffff", "#000000"]
    elif high_corr == 2:
        color_palette = ["#ffffff"] + list(reversed(list(Colorblind[3])))[1:]
    elif high_corr <= 8:
        color_palette = ["#ffffff"] + list(Colorblind[high_corr])
    else:
        color_palette = ["#ffffff"] + list(reversed(linear_palette(Reds256, high_corr)))
    # color_mapper = LinearColorMapper(palette=color_palette, low=0, high=high_corr)

    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"
    p = figure(x_range=list(correlation_matrix.columns), y_range=list(correlation_matrix.index),
               width=1920, height=1080, tools=tools, toolbar_location='below')
    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.grid.grid_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_standoff = 1
    p.xaxis.axis_label = 'Module ID'
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = 'System name'
    p.yaxis.major_label_orientation = 1 / 3

    r = p.rect('modules', 'systems', 1, 1, source=source, line_color="white",
               fill_color=linear_cmap('corr', palette=color_palette, low=0, high=high_corr))

    p.add_layout(r.construct_color_bar(label_standoff=12,
                                       ticker=BasicTicker(desired_num_ticks=len(color_palette)-1),
                                       border_line_color=None,
                                       location=(0, 0)),
                 'right')

    if "html" in out_format:
        output_path = Path.cwd() / output / name / f"{name}_correlation_module.html"
        output_file(output_path)
        save(p)
        logging.getLogger("PANORAMA").debug(f"Saved partition heatmap in HTML format to {output_path}")
    if "png" in out_format:
        output_path = Path.cwd() / output / name / f"{name}_correlation_module.png"
        export_png(p, width=1920, height=1080)
        logging.getLogger("PANORAMA").debug(f"Saved partition heatmap in PNG format to {output_path}")


def association_pangenome_systems(pangenome: Pangenome, association: List[str], output: Path,
                                  out_format: List[str] = None):
    """Write association between systems and pangenome object
    """
    out_format = out_format if out_format is not None else ['html']

    association_df = get_association_df(pangenome, association)
    association_df.to_csv(output / 'association.csv', sep='\t')
    write_correlation_matrix(association_df, output, pangenome.name, out_format)
