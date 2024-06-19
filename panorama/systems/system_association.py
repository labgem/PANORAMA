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
from bokeh.transform import linear_cmap
from bokeh.palettes import Colorblind, Reds256, linear_palette
from bokeh.models import BasicTicker

# local libraries
from panorama.pangenomes import Pangenome


def get_association_df(pangenome: Pangenome, association: List[str]) -> pd.DataFrame:
    columns = ['system number', 'system_name', 'families']
    if 'RGPs' in association:
        columns.append('RGPs')
    if 'spots' in association:
        columns.append('spots')
    if 'modules' in association:
        columns.append('modules')
    association_list = {}
    for system in pangenome.systems:
        association_list[system.ID] = [system.name, ",".join([fam.name for fam in system.families])]
        if 'RGPs' in association:
            association_list[system.ID].append(",".join(map(str, system.regions)))
        if 'spots' in association:
            association_list[system.ID].append(",".join(map(lambda x: str(x.ID), system.spots)))
        if 'modules' in association:
            association_list[system.ID].append(",".join(map(lambda x: str(x.ID), system.modules)))

    association_df = pd.DataFrame.from_dict(association_list, orient='index', columns=columns[1:])
    association_df.index.name = columns[0]
    return association_df


def write_correlation_matrix(df: pd.DataFrame, association: str, output: Path, name: str, out_format: List[str] = None):
    out_format = out_format if out_format is not None else ['html']

    association_split = df.drop(columns=['families']).join(df[association].str.get_dummies(sep=','))
    association_split = association_split.drop(columns=[association])

    correlation_matrix = association_split.groupby('system_name').sum()
    correlation_matrix.sort_index(key=lambda x: x.str.lower(), ascending=False, inplace=True)
    correlation_matrix.columns.name = association
    if association == 'RGPs':
        correlation_matrix = correlation_matrix[sorted(correlation_matrix.columns.tolist())]
        xlabel = "RGP name"
    else:
        correlation_matrix = correlation_matrix[sorted(correlation_matrix.columns.tolist(), key=lambda x: int(x))]
        xlabel = f"{association} ID"

    high_corr = correlation_matrix.values.max()
    if high_corr == 1:
        color_palette = ["#ffffff", "#000000"]
    elif high_corr == 2:
        color_palette = ["#ffffff"] + list(reversed(list(Colorblind[3])))[1:]
    elif high_corr <= 8:
        color_palette = ["#ffffff"] + list(Colorblind[high_corr])
    else:
        color_palette = ["#ffffff"] + list(reversed(linear_palette(Reds256, high_corr + 4)))[4:]

    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"
    p = figure(x_range=list(correlation_matrix.columns), y_range=list(correlation_matrix.index),
               width=1720, height=960, tools=tools, toolbar_location='below',
               tooltips=[(association, f'@{association}'), ('system', '@system_name')])
    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.grid.grid_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_standoff = 1
    p.xaxis.axis_label = xlabel
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = 'System name'
    p.yaxis.major_label_orientation = 1 / 3

    source = pd.DataFrame(correlation_matrix.stack(), columns=['corr']).reset_index()

    r = p.rect(association, 'system_name', 1, 1, source=source, line_color="white",
               fill_color=linear_cmap('corr', palette=color_palette, low=0, high=high_corr + 1))

    p.add_layout(r.construct_color_bar(label_standoff=12,
                                       ticker=BasicTicker(desired_num_ticks=len(color_palette)),
                                       border_line_color=None,
                                       location=(0, 0)),
                 'right')

    if "html" in out_format:
        output_path = Path.cwd() / output / name / f"{name}_correlation_{association}.html"
        output_file(output_path)
        save(p)
        logging.getLogger("PANORAMA").debug(f"Saved partition heatmap in HTML format to {output_path}")
    if "png" in out_format:
        output_path = Path.cwd() / output / name / f"{name}_correlation_{association}.png"
        export_png(p, width=1920, height=1080)
        logging.getLogger("PANORAMA").debug(f"Saved partition heatmap in PNG format to {output_path}")


def association_pangenome_systems(pangenome: Pangenome, association: List[str], output: Path,
                                  out_format: List[str] = None):
    """Write association between systems and pangenome object
    """
    out_format = out_format if out_format is not None else ['html']

    association_df = get_association_df(pangenome, association)
    association_df.to_csv(output / 'association.csv', sep='\t')
    for asso in association:
        write_correlation_matrix(association_df.drop([other for other in association if other != asso], axis=1),
                                 asso, output, pangenome.name, out_format)
