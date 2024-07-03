#!/usr/bin/env python3
# coding:utf-8

"""
This module allows the association of systems to other pangenome elements.
"""

# default libraries
from __future__ import annotations
import logging
from collections import defaultdict, namedtuple
from pathlib import Path
from typing import Dict, List, Set, Tuple, Union

# installed libraries
import pandas as pd
from bokeh.io import output_file, save, export_png
from bokeh.layouts import gridplot, row
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.palettes import Colorblind, Reds256, Blues256, linear_palette
from bokeh.models import BasicTicker, ColumnDataSource, LinearColorMapper, ColorBar, GlyphRenderer, FactorRange

# local libraries
from panorama.pangenomes import Pangenome
from panorama.region import Region, Spot, Module
from panorama.systems.system import System

total_width, total_height = 1780, 920
left_width, center_width, right_width = [int(x * total_width) for x in [0.15, 0.75, 0.1]]
top_height, middle_height, below_height = [int(y * total_height) for y in [0.15, 0.7, 0.15]]


def get_coverage_df(asso2sys: Dict[Union[Region, Spot, Module], Set[System]],
                    pangenome: Pangenome = None) -> pd.DataFrame:
    field = ["name", "systems_ID", "systems_name", "coverage"]
    if pangenome is not None:
        field.append("frequency")
    elem_out = namedtuple("Elem", field)
    out = []
    for element, systems in asso2sys.items():
        element_fam = set(element.families)
        sys_fam, sys_org, sys_id, sys_name = set(), set(), set(), set()
        for sys in systems:
            sys_fam |= set(sys.families)
            sys_org |= set(sys.organisms)
            sys_id.add(sys.ID)
            sys_name.add(sys.name)
        coverage = len(element_fam.intersection(sys_fam)) / len(element_fam)
        if pangenome is not None:
            element_org = set(element.organisms)
            freq = len(element_org.intersection(sys_org)) / pangenome.number_of_organisms
            out.append(elem_out(str(element), ",".join(sys_id), ','.join(sys_name), coverage, freq))
        else:
            out.append(elem_out(str(element), ",".join(sys_id), ','.join(sys_name), coverage))
    return pd.DataFrame(out)


def get_association_df(pangenome: Pangenome, association: List[str]
                       ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Get the dataframe corresponding to the system-pangenome object association.

    Args:
        pangenome (Pangenome): Pangenome containing systems.
        association (List[str]): List of pangenome elements to associate.

    Returns:
        pd.DataFrame: DataFrame corresponding to the system-pangenome element association.
    """
    columns = ['system number', 'system_name', 'families']
    if 'RGPs' in association:
        columns.append('RGPs')
    if 'spots' in association:
        columns.append('spots')
    if 'modules' in association:
        columns.append('modules')

    association_list = {}
    rgp2sys = defaultdict(set)
    spot2sys = defaultdict(set)
    mod2sys = defaultdict(set)
    for system in pangenome.systems:
        association_list[system.ID] = [system.name, ",".join([fam.name for fam in system.families])]
        if 'RGPs' in association:
            rgps = set()
            for rgp in system.regions:
                rgps.add(rgp.name)
                rgp2sys[rgp].add(system)
            association_list[system.ID].append(",".join(rgps))
        if 'spots' in association:
            spots = set()
            for spot in system.spots:
                spots.add(str(spot.ID))
                spot2sys[spot].add(system)
            association_list[system.ID].append(",".join(spots))
        if 'modules' in association:
            modules = set()
            for mod in system.modules:
                modules.add(str(mod.ID))
                mod2sys[mod].add(system)
            association_list[system.ID].append(",".join(modules))
    association_df = pd.DataFrame.from_dict(association_list, orient='index', columns=columns[1:])
    association_df.index.name = columns[0]
    rgp2sys_df = get_coverage_df(rgp2sys)
    spot2sys_df = get_coverage_df(spot2sys, pangenome)
    mod2sys_df = get_coverage_df(mod2sys, pangenome)
    return association_df, rgp2sys_df, spot2sys_df, mod2sys_df


def preprocess_data(df: pd.DataFrame, association: str) -> pd.DataFrame:
    """
    Preprocess the data for the correlation matrix.

    Args:
        df (pd.DataFrame): Association dataframe between systems and pangenome objects.
        association (str): Pangenome object to associate systems with.

    Returns:
        pd.DataFrame: Preprocessed correlation matrix.
    """
    association_split = df.drop(columns=['families']).join(df[association].str.get_dummies(sep=','))
    association_split = association_split.drop(columns=[association])

    correlation_matrix = association_split.groupby('system_name').sum()
    correlation_matrix.sort_index(key=lambda x: x.str.lower(), ascending=False, inplace=True)
    correlation_matrix.columns.name = association

    return correlation_matrix


def create_figure(correlation_matrix: pd.DataFrame, association: str,
                  x_range: FactorRange, y_range: FactorRange) -> Tuple[figure, GlyphRenderer]:
    """
    Create the figure for the correlation matrix.

    Args:
        correlation_matrix (pd.DataFrame): Preprocessed correlation matrix.
        association (str): Pangenome object to associate systems with.
        x_range: Range of x abscissa
        y_range: Range of y abscissa

    Returns:
        figure: Bokeh figure object.
    """
    high_corr = correlation_matrix.values.max()
    if high_corr == 1:
        color_palette = ["#ffffff", "#000000"]
    elif high_corr == 2:
        color_palette = ["#ffffff"] + list(reversed(list(Colorblind[3])))[1:]
    elif high_corr <= 8:
        color_palette = ["#ffffff"] + list(Colorblind[high_corr])
    else:
        color_palette = ["#ffffff"] + list(linear_palette(Reds256[::-1], high_corr + 4))[4:]

    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"
    p = figure(x_range=x_range, y_range=y_range, width=center_width, height=middle_height, tools=tools,
               toolbar_location='below', tooltips=[(association, f'@{association}'), ('system', '@system_name')])
    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.grid.grid_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_standoff = 1
    p.xaxis.axis_label = "RGP name" if association == 'RGPs' else f"{association} ID"
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = 'System name'
    p.yaxis.major_label_orientation = 1 / 3

    source = ColumnDataSource(correlation_matrix.stack().reset_index(name='corr'))

    r = p.rect(association, 'system_name', 1, 1, source=source, line_color="white",
               fill_color=linear_cmap('corr', palette=color_palette, low=0, high=high_corr + 1))

    return p, r


def create_color_bar(r: GlyphRenderer) -> figure:
    """
    Create the color bar for the correlation matrix.

    Args:
        r (figure): Rectangle glyph object.

    Returns:
        figure: Bokeh figure object for the color bar.
    """
    color_bar = ColorBar(color_mapper=r.glyph.fill_color['transform'],
                         label_standoff=12,
                         ticker=BasicTicker(desired_num_ticks=len(r.glyph.fill_color['transform'].palette)),
                         border_line_color=None)

    color_bar_plot = figure(title="#Systems/spots", title_location="right",
                            height=middle_height, width=right_width,
                            toolbar_location=None, min_border=0,
                            outline_line_color=None)

    color_bar_plot.add_layout(color_bar, 'right')
    color_bar_plot.title.align = "center"
    color_bar_plot.title.text_font_size = '14pt'

    return color_bar_plot


def create_coverage_plot(coverage: pd.DataFrame, association: str, x_range: FactorRange) -> Tuple[figure, figure]:
    """
    Create the coverage plot for the correlation matrix.

    Args:
        association (str): Pangenome object to associate systems with.
        coverage (pd.DataFrame): Coverage dataframe.
        x_range: order of x elements

    Returns:
        figure: Bokeh figure object for the coverage plot.
    """
    coverage_mapper = LinearColorMapper(palette=Reds256, low=1, high=0)
    coverage_color_bar = ColorBar(color_mapper=coverage_mapper, label_standoff=12,
                                  ticker=BasicTicker(desired_num_ticks=5),
                                  border_line_color=None,
                                  location=(0, 0))

    coverage_color_bar_plot = figure(title="Coverage", title_location="below",
                                     height=int(below_height / 3), width=int(center_width / 2.5),
                                     toolbar_location=None, min_border=0,
                                     outline_line_color=None)

    coverage_color_bar_plot.add_layout(coverage_color_bar, 'below')
    coverage_color_bar_plot.title.align = "center"
    coverage_color_bar_plot.title.text_font_size = '12pt'
    coverage_source = ColumnDataSource(pd.DataFrame({
        association: x_range.factors,
        'coverage': coverage.loc[x_range.factors]["coverage"]
    }))

    coverage_p = figure(x_range=x_range, height=int(below_height / 4), width=center_width,
                        tooltips=[(association, f'@{association}'), ('coverage', '@coverage')])
    coverage_p.xaxis.visible = False
    coverage_p.yaxis.visible = False
    coverage_p.grid.grid_line_color = None
    coverage_p.outline_line_color = None

    coverage_p.rect(association, 0.5, 1, 1, source=coverage_source,
                    fill_color={'field': 'coverage', 'transform': coverage_mapper})

    return coverage_p, coverage_color_bar_plot


def create_frequency_plot(frequency: pd.DataFrame, association: str, x_range: FactorRange) -> Tuple[figure, figure]:
    """
    Create the frequency plot for the correlation matrix.

    Args:
        association (str): Pangenome object to associate systems with.
        frequency (pd.DataFrame): Coverage dataframe.
        x_range: order of x elements


    Returns:
        figure: Bokeh figure object for the coverage plot.
    """
    frequency_mapper = LinearColorMapper(palette=Blues256, low=1, high=0)
    frequency_color_bar = ColorBar(color_mapper=frequency_mapper, label_standoff=12,
                                   ticker=BasicTicker(desired_num_ticks=5),
                                   border_line_color=None,
                                   location=(0, 0))

    frequency_color_bar_plot = figure(title="Genome frequencies", title_location="below",
                                      height=int(below_height / 3), width=int(center_width / 2.5),
                                      toolbar_location=None, min_border=0,
                                      outline_line_color=None)

    frequency_color_bar_plot.add_layout(frequency_color_bar, 'below')
    frequency_color_bar_plot.title.align = "center"
    frequency_color_bar_plot.title.text_font_size = '12pt'

    frequency_source = ColumnDataSource(pd.DataFrame({
        association: x_range.factors,
        'frequency': frequency.loc[x_range.factors]["frequency"]
    }))

    frequency_p = figure(x_range=x_range, height=int(below_height / 4), width=center_width,
                         tooltips=[(association, f'@{association}'), ('frequency', '@frequency')])
    frequency_p.xaxis.visible = False
    frequency_p.yaxis.visible = False
    frequency_p.grid.grid_line_color = None
    frequency_p.outline_line_color = None

    frequency_p.rect(association, 0.5, 1, 1, source=frequency_source,
                     fill_color={'field': 'frequency', 'transform': frequency_mapper})

    return frequency_p, frequency_color_bar_plot


def create_bar_plots(correlation_matrix: pd.DataFrame, association: str) -> Tuple[figure, figure, List[str]]:
    """
    Create the bar plots for the correlation matrix.

    Args:
        correlation_matrix (pd.DataFrame): Preprocessed correlation matrix.
        association (str): Pangenome object to associate systems with.

    Returns:
        figure: Bokeh figure objects for the bar plots.
    """
    system_counts = correlation_matrix.sum(axis=1)
    left_bar_source = ColumnDataSource(pd.DataFrame({
        f"left_bar_{association}": correlation_matrix.index.values,
        'count': system_counts
    }))

    left_bar = figure(y_range=list(correlation_matrix.index), width=left_width, height=middle_height,
                      toolbar_location=None, tools="")
    left_bar.hbar(y='system_name', right='count', height=0.9, source=left_bar_source,
                  color='navy', alpha=0.6)
    left_bar.yaxis.visible = False
    left_bar.xaxis.axis_label = 'Count'
    left_bar.x_range.flipped = True
    left_bar.grid.grid_line_color = None
    left_bar.outline_line_color = None

    individual_counts = correlation_matrix.sum(axis=0)
    top_bar_source = ColumnDataSource(pd.DataFrame({
        f"top_bar_{association}": correlation_matrix.columns.values,
        'count': individual_counts
    }))
    x_ord = list(sorted(correlation_matrix.columns, key=lambda x: individual_counts.loc[x], reverse=True))
    top_bar = figure(x_range=x_ord, height=top_height, width=center_width,
                     toolbar_location=None, tools="")
    top_bar.vbar(x=association, top='count', width=0.9, source=top_bar_source,
                 color='green', alpha=0.6)
    top_bar.xaxis.visible = False
    top_bar.yaxis.axis_label = 'Count'
    top_bar.grid.grid_line_color = None
    top_bar.outline_line_color = None

    return left_bar, top_bar, x_ord


def write_correlation_matrix(df: pd.DataFrame, association: str, coverage: pd.DataFrame, output: Path,
                             frequency: pd.DataFrame = None, out_format: List[str] = None):
    """
    Write the correlation matrix.

    Args:
        df (pd.DataFrame): Association dataframe between systems and pangenome objects.
        association (str): Pangenome object to associate systems with.
        coverage (pd.DataFrame): Coverage dataframe between systems and pangenome objects.
        output (Path): Path to the output directory.
        out_format (List[str], optional): Formats of the output file (default is ['html']).
    """
    out_format = out_format if out_format is not None else ['html']

    correlation_matrix = preprocess_data(df, association)
    left_bar, top_bar, xord = create_bar_plots(correlation_matrix, association)
    p, r = create_figure(correlation_matrix, association, top_bar.x_range, left_bar.y_range)
    color_bar_plot = create_color_bar(r)
    coverage_p, coverage_color_bar_plot = create_coverage_plot(coverage, association, top_bar.x_range)
    if frequency is not None:
        frequency_p, frequency_color_bar_plot = create_frequency_plot(frequency, association, top_bar.x_range)
        p = gridplot([[None, top_bar, None, None],
                      [left_bar, p, color_bar_plot],
                      [None, frequency_p, None],
                      [None, coverage_p, None],
                      [None, row([frequency_color_bar_plot, coverage_color_bar_plot], align='end', spacing=100), None]],
                     toolbar_location='above')
    else:
        p = gridplot([[None, top_bar, None],
                      [left_bar, p, color_bar_plot],
                      [None, coverage_p, None],
                      [None, row(coverage_color_bar_plot, align="center"), None]],
                     toolbar_location='above')

    if "html" in out_format:
        output_path = output / f"correlation_{association}.html"
        output_file(output_path)
        save(p)
        logging.getLogger("PANORAMA").debug(f"Saved partition heatmap in HTML format to {output_path}")
    if "png" in out_format:
        output_path = output / f"correlation_{association}.png"
        export_png(p, filename=output_path, width=total_width, height=total_height)
        logging.getLogger("PANORAMA").debug(f"Saved partition heatmap in PNG format to {output_path}")


def association_pangenome_systems(pangenome: Pangenome, association: List[str], output: Path,
                                  out_format: List[str] = None):
    """
    Write the association between systems and pangenome objects.

    Args:
        pangenome (Pangenome): The pangenome containing systems and other elements.
        association (List[str]): List of pangenome elements to associate.
        output (Path): Path to the output directory.
        out_format (List[str], optional): Formats of the output files (default is ['html']).

    Returns:
        None
    """
    out_format = out_format if out_format is not None else ['html']

    association_df, rgp2sys_df, spot2sys_df, mod2sys_df = get_association_df(pangenome, association)
    association_df.to_csv(output / 'association.tsv', sep='\t')
    rgp2sys_df.to_csv(output / 'rgp_to_systems.tsv', sep='\t')
    spot2sys_df.to_csv(output / 'spot_to_systems.tsv', sep='\t')
    mod2sys_df.to_csv(output / 'module_to_systems.tsv', sep='\t')
    for asso in association:
        if asso == "RGPs":
            coverage = rgp2sys_df.set_index("name")
            frequency = None
        elif asso == "spots":
            coverage = spot2sys_df.set_index("name")
            coverage.index = coverage.index.str.replace("spot_", "")
            frequency = coverage.loc[:, coverage.columns != "coverage"]
            coverage = coverage.loc[:, coverage.columns != "frequency"]
        elif asso == "modules":
            coverage = mod2sys_df.set_index("name")
            coverage.index = coverage.index.str.replace("module_", "")
            frequency = coverage.loc[:, coverage.columns != "coverage"]
            coverage = coverage.loc[:, coverage.columns != "frequency"]
        else:
            raise Exception("Unexpected error")
        write_correlation_matrix(association_df.drop([other for other in association if other != asso], axis=1),
                                 asso, coverage, output, frequency, out_format)
