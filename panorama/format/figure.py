#!/usr/bin/env python3
# coding:utf-8

import numpy as np
from typing import List, Tuple
from bokeh.models import (BasicTicker, ColorBar, LinearColorMapper, CategoricalColorMapper,
                          PrintfTickFormatter, ColumnDataSource)
from bokeh.palettes import Magma256
import pandas as pd
from bokeh.io import output_file, export_png
from bokeh.plotting import figure, save
from bokeh.layouts import gridplot
from bokeh.palettes import turbo
from pathlib import Path
from bokeh.models import FactorRange, CustomJS, Button


def per_pan_heatmap(name: str, system_projection: pd.DataFrame, output: Path):
    """ Call functions to draw heatmap figures for the pangenome

    :param name: Name of the pangenome
    :param system_projection: Systems in the pangenome
    :param output: Path to output directory
    """

    data_systems = system_projection[['system number', 'system name', 'organism', 'partition']]
    system_names = data_systems['system name'].unique().tolist()
    system_names.sort(key=str.casefold)
    # We have a sorted list with system names without duplicates
    organism_names = data_systems['organism'].unique().tolist()
    organism_names.sort(key=str.casefold)
    organism_names.sort(key=lambda x: ('assembled' in x or 'MAGS' in x, x), reverse=True)
    # We have a sorted list with organism names without duplicates

    matrix_genomes_systems = np.zeros((len(organism_names), len(system_names)))
    for i, organism in enumerate(organism_names):
        data_genome = data_systems[(data_systems['organism'] == organism)]
        dict_defense_genome = pd.Series(data_genome['system name'].values, index=data_genome['system number']).to_dict()
        for j, system in enumerate(system_names):
            matrix_genomes_systems[i][j] = sum(
                x == system for x in dict_defense_genome.values())  # Sum of system number for each system

    figure_partition_heatmap(name=name, data=data_systems, list_systems=system_names, list_organisms=organism_names,
                             output=output)
    figure_count_heatmap(name=name, data=matrix_genomes_systems, list_system=system_names, list_organism=organism_names,
                         output=output)


def figure_partition_heatmap(name: str, data: pd.DataFrame, list_systems: List[str],
                             list_organisms: List[str], output: Path):
    """ Draw partition heatmap figure for the pangenome

    :param name: Name of the pangenome
    :param data: Data used to produce the heatmap
    :param list_systems: List of systems in the pangenome
    :param list_organisms: List of organisms in the pangenome
    :param output: Path to output directory
    """

    output_file(Path.cwd() / output / name / f"{name}_partition.html")
    data_partition = data.pivot_table(index='organism', columns='system name', values='partition',
                                      fill_value='Not_found', aggfunc=lambda x: ','.join(set(x)))
    df_stack_partition = pd.DataFrame(data_partition.stack(), columns=['partition']).reset_index()

    colors = ["#FF9933", "#3DCA2E", "#B1F1EB", "#FFFFFF"]
    partitions = ["persistent", "shell", "cloud", "Not_found"]
    mapper = CategoricalColorMapper(palette=colors, factors=partitions)
    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

    p = figure(title="Partition of systems for {0}".format(name), x_range=list_systems, y_range=list_organisms,
               x_axis_location="above", tools=tools, toolbar_location='below',
               tooltips=[('Partition', '@partition'), ('Organism', '@organism'), ('System', '@{system name}')])

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Systems name'
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = 'Genomes name'
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.xgrid.grid_line_color = None

    p.rect(x="system name", y="organism", width=1, height=1, source=df_stack_partition,
           fill_color={'field': 'partition', 'transform': mapper}, line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="14px", label_standoff=6,
                         border_line_color=None,
                         major_label_overrides={0: "Persistent", 1: "Shell", 2: "Cloud", 3: "Not_found"},
                         major_tick_line_color='black', bar_line_color='black', bar_line_width=0.2,
                         border_line_width=0.2)

    p.add_layout(color_bar, 'right')
    save(p)


def figure_count_heatmap(name: str, data: np.ndarray, list_system: List, list_organism: List, output: Path):
    """ Draw count heatmap figure for the pangenome

    :param name: Name of the pangenome
    :param data: Data used to produce the heatmap
    :param list_system: List of systems in the pangenome
    :param list_organism: List of organisms in the pangenome
    :param output: Path to output directory
    """

    output_file(Path.cwd() / output / name / "{0}_count.html".format(name))
    df_gcount = pd.DataFrame(data, index=list_organism, columns=list_system)
    df_stack_gcount = pd.DataFrame(df_gcount.stack(), columns=['number']).reset_index()

    mapper = LinearColorMapper(palette=list(reversed(Magma256)), low=1, low_color='#FFFFFF',
                               high=df_stack_gcount.number.max())
    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

    p = figure(title="Systems for {0}".format(name), x_range=list_system, y_range=list_organism,
               x_axis_location="above", tools=tools, toolbar_location='below',
               tooltips=[('Count', '@number'), ('Organism', '@level_0'), ('System', '@level_1')])

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Systems name'
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = 'Genomes name'
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.xgrid.grid_line_color = None

    p.rect(x="level_1", y="level_0", width=1, height=1, source=df_stack_gcount,
           fill_color={'field': 'number', 'transform': mapper}, line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="14px",
                         ticker=BasicTicker(desired_num_ticks=int(data.max())),
                         formatter=PrintfTickFormatter(format="%s"), label_standoff=6, border_line_color=None)
    p.add_layout(color_bar, 'right')
    save(p)


def heatmap(global_systems_proj: pd.DataFrame, output: Path):
    """ Transform data for heatmap and histogram figures used for all pangenomes

    :param global_systems_proj: Systems detected in all pangenomes
    :param output: Path to output directory
    """

    names = global_systems_proj['pangenome name'].unique().tolist()
    names.sort(key=str.casefold)
    systems_type = global_systems_proj['system name'].unique().tolist()
    systems_type.sort(key=str.casefold)

    save_sum_pans_systems = []
    save_sum_pans_ratio = []
    for i, name in enumerate(names):
        sum_pans_systems = {}
        sum_pans_ratio = {}
        pan_df = global_systems_proj[global_systems_proj['pangenome name'] == name]
        print("pika", pan_df['organism'].unique().shape)
        number_org_tot = pan_df['organism'].unique().shape[0]
        for system in systems_type:
            system_df = pan_df[pan_df['system name'] == system]
            print("cara", system_df['organism'].unique().shape[0])
            number_org = system_df['organism'].unique().shape[0]
            sum_pans_systems[system] = number_org
            sum_pans_ratio[system] = number_org / number_org_tot

        save_sum_pans_systems.append(sum_pans_systems)
        save_sum_pans_ratio.append(sum_pans_ratio)

    df_pan_count = pd.DataFrame(save_sum_pans_systems, columns=systems_type, index=names).fillna(0)
    df_pan_ratio = pd.DataFrame(save_sum_pans_ratio, columns=systems_type, index=names).fillna(0)

    figure_histogram(names=names, df_pan_count=df_pan_count, systems_type=systems_type, output=output)
    figure_heatmap(names=names, df_pan_ratio=df_pan_ratio, output=output)


def figure_histogram(names: List[str], df_pan_count: pd.DataFrame, systems_type: List[str], output: Path):
    """ Draw histogram figure for all pangenomes

    :param names: Names of all pangenomes
    :param df_pan_count: Count of all types of systems in all pangenomes
    :param systems_type: Systems type in all pangenomes
    :param output: Path to output directory
    """

    output_file(Path.cwd() / output / "bar_stacked_pangenomes.html")

    dict_graph_pans = [(key, value) for i, (key, value) in enumerate(zip(names, df_pan_count.values.tolist()))]
    dict_pan_names = {'systems_names': systems_type}
    dict_pan_names.update(dict_graph_pans)
    col_pan_names = ColumnDataSource(data=dict_pan_names)

    p = figure(title="Systems count", background_fill_color="#FAFAFA", x_range=systems_type,
               sizing_mode='scale_both', x_axis_location="below", tools="hover", toolbar_location=None,
               tooltips=[('Pangenome', '$name'), ('system', '@systems_names'), ('Count', '@$name')])

    p.vbar_stack(names, x='systems_names', width=0.5, color=turbo(len(names)),
                 source=col_pan_names, legend_label=names, line_color='black')

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Systems'
    p.xaxis.major_label_orientation = 0.7
    p.x_range.range_padding = 0.1
    p.yaxis.axis_label = 'Count'
    p.y_range.start = 0
    p.xgrid.grid_line_color = None
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.axis.minor_tick_line_color = None
    p.legend.location = "top_center"
    p.legend.orientation = "horizontal"
    p.outline_line_color = None
    p.legend.click_policy = "hide"

    save(p)


def figure_heatmap(names: List[str], df_pan_ratio: pd.DataFrame, output: Path, threshold: float = 0.1):
    """ Draw heatmap figure for all pangenomes

    :param names: Names of all pangenomes
    :param df_pan_ratio: Ratio (count types of systems/count organisms) in all pangenomes
    :param output: Path to output directory
    """
    df_pan_ratio_f = df_pan_ratio.loc[:, (df_pan_ratio >= threshold).any()]  # remove systems present in less that 10% of each pangenome
    print(df_pan_ratio, df_pan_ratio_f)
    systems_list = df_pan_ratio_f.columns.tolist()
    systems_list.sort()

    output_file(Path.cwd() / output / "pangenomes_count.html")
    df_stack_pan_ratio = pd.DataFrame(df_pan_ratio_f.stack(), columns=['ratio']).reset_index()
    mapper = LinearColorMapper(palette=(list(reversed(Magma256))), low=0.000001, low_color='#FFFFFF',
                               high=df_stack_pan_ratio.ratio.max())

    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

    len_pangenomes = len(names)
    len_systems = len(systems_list)
    height = len_pangenomes * 100
    width = len_systems * 50

    p = figure(title="Systems in pangenomes", x_range=systems_list, y_range=names,
               x_axis_location="above", height=height, width=width, tools=tools, toolbar_location='below',
               tooltips=[('pangenome', '@level_0'), ('System', '@level_1'), ('Ratio', '@ratio{0.00}')])

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Systems name'
    p.xaxis.major_label_orientation = 0.7
    p.yaxis.axis_label = 'Pangenomes name'
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.xgrid.grid_line_color = None

    p.rect(x='level_1', y='level_0', width=1, height=1, source=df_stack_pan_ratio,
           fill_color={'field': 'ratio', 'transform': mapper}, line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="14px",
                         ticker=BasicTicker(desired_num_ticks=int(df_stack_pan_ratio.ratio.max())),
                         formatter=PrintfTickFormatter(format="%s"), label_standoff=6, border_line_color=None)
    p.add_layout(color_bar, 'right')
    save(p)


def upsetplot(name: str, systems_projection: pd.DataFrame, system_to_feature: pd.DataFrame, output: Path):
    """ Draw upsetplot figure (presence/absence of systems in spots or modules predicted by PPanGGOLiN)

    :param name: Name of the pangenome
    :param systems_projection: Dataframe with systems detected in the pangenome
    :param system_to_feature: Dataframe with modules, spots and RGPs associated to each system ID of pangenome
    :param output: Path to output directory
    """

    systems_projection = systems_projection[['system number', 'system name', 'organism']]
    system_type_list = systems_projection['system name'].unique().tolist()
    system_type_list.sort(key=str.lower)
    organism_name_list = systems_projection['organism'].unique().tolist()
    organism_name_list.sort(key=lambda x: ('assembled' in x or 'MAGS' in x, x), reverse=True)
    organism_id_pan = [i + 1 for i in range(len(organism_name_list))]
    dict_organism_id_pan = {name: ID + 1 for ID, name in enumerate(organism_name_list)}
    system_to_feature = system_to_feature[['system_name', 'mod_organism', 'module', 'rgp_organism', 'spot', 'rgp']]

    systems_count, organisms_count = count(systems_projection=systems_projection, system_type_list=system_type_list,
                                           organism_name_list=organism_name_list)

    spot_set = set()
    dict_spot_system = {}
    dict_spot_organism_id = {}

    module_set = set()
    dict_module_system = {}
    dict_module_organism_id = {}

    x_axis_systems_present = []
    y_axis_organisms_id_present = []

    for system in system_type_list:
        systems_projection_filter = systems_projection[systems_projection['system name'] == system]
        organism_id = 1
        for organism in dict_organism_id_pan.keys():
            checking_presence = (systems_projection_filter['organism'].eq(organism)).to_string()
            if 'True' in checking_presence:
                x_axis_systems_present.append(system)
                y_axis_organisms_id_present.append(organism_id)
            organism_id += 1

        if system_to_feature['system_name'].str.contains(system).any():
            data_features_filter = system_to_feature[system_to_feature['system_name'] == system]
            for org_index, row in data_features_filter.iterrows():
                rgp_organisms = row.rgp_organism
                mod_organisms = row.mod_organism
                spots = row.spot
                modules = row.module
                if spots is not None:
                    spot_set |= set(spots)
                    spot_set.discard("")
                if modules is not None:
                    module_set |= set(modules)
                if spots is not None:
                    for i, spot in enumerate(spots):
                        dict_spot_system.setdefault(spot, [])
                        dict_spot_system[spot].append(system)
                        organism_id_rgp = dict_organism_id_pan[rgp_organisms[i]]
                        dict_spot_organism_id.setdefault(spot, [])
                        dict_spot_organism_id[spot].append(organism_id_rgp)

                if modules is not None:
                    for j, module in enumerate(modules):
                        dict_module_system.setdefault(module, [])
                        dict_module_system[module].append(system)
                        organism_id_mod = dict_organism_id_pan[mod_organisms[j]]
                        dict_module_organism_id.setdefault(module, [])
                        dict_module_organism_id[module].append(organism_id_mod)

    tooltips = [('System name', '@x'), ('Features', '$name')]

    # Set size for plot to show the whole legend
    len_sys = len(system_type_list)
    len_figure = len(spot_set) + len(module_set) + len(organism_name_list)
    width_center = len_sys * 40

    if len_figure < 30:
        height_center = len_figure * 25
        height_top = 400
        width_left = 350
        size_spot = 20
        size_module = 17
    elif 30 <= len_figure <= 100:
        height_center = len_figure * 28
        height_top = 450
        width_left = 400
        size_spot = 21
        size_module = 18
    elif 100 < len_figure <= 200:
        height_center = len_figure * 30
        height_top = 475
        width_left = 425
        size_spot = 23
        size_module = 20
    else:
        height_center = len_figure * 32
        height_top = 500
        width_left = 450
        size_spot = 25
        size_module = 21

    p_features = figure(width=width_center, height=height_center, x_range=system_type_list, tooltips=tooltips,
                        x_axis_location="above")
    p_features.square(x=x_axis_systems_present, y=y_axis_organisms_id_present, fill_color="white",
                      line_color="black", legend_label="Alone", size=size_spot, name="")

    spot_list = sorted(list(spot_set), key=float)
    color_spots = turbo(len(spot_list))
    for i, spot in enumerate(spot_list):
        x_axis_systems_spot = dict_spot_system[spot]
        y_axis_organisms_id_spot = dict_spot_organism_id[spot]
        p_features.square(x=x_axis_systems_spot, y=y_axis_organisms_id_spot, line_color="black",
                          fill_color=color_spots[i], size=size_spot,
                          legend_label="Spot " + str(spot), name=str(spot))

    module_list = sorted(list(module_set), key=float)
    color_modules = turbo(len(module_list))
    for i, module in enumerate(module_list):
        x_axis_systems_module = dict_module_system[module]
        y_axis_organisms_id_module = dict_module_organism_id[module]
        p_features.circle(x=x_axis_systems_module, y=y_axis_organisms_id_module, line_color="black",
                          fill_color=color_modules[i], size=size_module,
                          legend_label="Module " + str(module), name=str(module))

    # Presence Absence of spot and module

    p_features.yaxis.ticker.desired_num_ticks = len(organism_id_pan)
    p_features.yaxis.minor_tick_line_color = None
    dict_organism_id_pan_reverse = dict((v, k) for k, v in dict_organism_id_pan.items())
    p_features.yaxis.major_label_overrides = dict_organism_id_pan_reverse
    p_features.x_range.range_padding = 0.1
    p_features.y_range.range_padding = 0.1
    p_features.y_range.start = 0.5
    p_features.y_range.end = len(organism_id_pan) + 0.5
    p_features.xaxis.major_label_orientation = 0.8
    p_features.outline_line_width = 1
    p_features.outline_line_color = "black"
    p_features.axis.major_label_text_font_size = "14px"
    p_features.yaxis.major_label_text_align = 'center'
    p_features.min_border_bottom = 100
    p_features.add_layout(p_features.legend[0], 'right')
    p_features.legend.click_policy = "hide"

    # Button to show/hide all spots and modules
    button = Button(label='Hide all', button_type="success", width=100, height=30, margin=(100, 10, 10, -115))
    cb = CustomJS(args=dict(fig=p_features, btn=button),
                  code='''
                  if (btn.label=='Hide all'){
                      for (var i=0; i<fig.renderers.length; i++){
                              fig.renderers[i].visible=false}
                      btn.label = 'Show all';
                      btn.button_type = "primary"
                      }
                  else {for (var i=0; i<fig.renderers.length; i++){
                          fig.renderers[i].visible=true}
                  btn.label = 'Hide all';
                  btn.button_type = "success"
                  }
                  ''')
    button.js_on_click(cb)

    # Histogram of organism
    tooltips_org = [('Organism', '@y'), ('Count', '@left')]

    p_org = figure(x_range=[max(organisms_count) + 1, 0], y_range=organism_name_list, width=width_left,
                   height=height_center,
                   tooltips=tooltips_org)
    p_org.hbar(y=organism_name_list, left=organisms_count, right=0, height=0.8, fill_color="#458B74",
               line_color="black")
    p_org.yaxis.visible = False
    p_org.outline_line_width = 1
    p_org.outline_line_color = "black"
    p_org.axis.major_label_text_font_size = "14px"
    p_org.xgrid.grid_line_color = "black"
    p_org.min_border_bottom = 100

    tooltips_sys = [('System name', '@x'), ('Count', '@top')]
    p_sys = figure(title="{0}".format(name), x_range=system_type_list, width=width_center, height=height_top,
                   tooltips=tooltips_sys)
    p_sys.vbar(x=system_type_list, top=systems_count, width=0.8, fill_color="#B22222", line_color="black")

    p_sys.title.align = "center"
    p_sys.title.text_font_size = "20pt"
    p_sys.outline_line_width = 1
    p_sys.outline_line_color = "black"
    p_sys.x_range.range_padding = 0.1
    p_sys.ygrid.grid_line_color = "black"
    p_sys.y_range.start = 0
    p_sys.xaxis.visible = False
    p_sys.axis.major_label_text_font_size = "14px"

    output_file(Path.cwd() / output / name / f"upsetplot_{name}.html")

    grid = gridplot([[None, p_sys, None], [p_org, p_features, button]])
    save(grid)
    export_png(grid, filename=output / name / f"upsetplot_{name}.png")


def count(systems_projection: pd.DataFrame, system_type_list: List, organism_name_list: List) -> Tuple[List, List]:
    """ Count number of system per type or per organism

    :param systems_projection: Dataframe with systems detected in the pangenome
    :param system_type_list: All type of systems detected in the pangenome
    :param organism_name_list: All organisms name in the pangenome
    """

    systems_count = [None] * len(system_type_list)
    for i, system in enumerate(system_type_list):
        data_pangenome_filter = systems_projection[systems_projection['system name'] == system]
        data_pangenome_filter_drop = data_pangenome_filter.drop_duplicates(
            subset=['system number', 'system name', 'organism'])
        systems_count[i] = data_pangenome_filter_drop.shape[0]

    organisms_count = [None] * len(organism_name_list)
    for j, organism in enumerate(organism_name_list):
        data_pangenome_filter_2 = systems_projection[systems_projection['organism'] == organism]
        data_pangenome_filter_drop_2 = data_pangenome_filter_2.drop_duplicates(
            subset=['system number', 'system name', 'organism'])
        organisms_count[j] = data_pangenome_filter_drop_2.shape[0]
    return systems_count, organisms_count


def hbar_id_total(name: str, dataframe_id: pd.DataFrame, dataframe_total: pd.DataFrame, output: Path):
    """ Draw histogram figure with number of ID per type and total number of system per type

    :param name: Name of the pangenome
    :param dataframe_id: ID data used to produce the histogram
    :param dataframe_total: Total count data used to produce the histogram
    :param output: Path to output directory
    """
    data_systems = dataframe_id[['system name', 'Number of IDs']]
    system_names = data_systems['system name'].tolist()
    system_names.sort(key=str.casefold)
    number_id = data_systems['Number of IDs'].tolist()

    data_systems_org = dataframe_total[['system name', 'Total number of systems']]
    number_total = list(data_systems_org['Total number of systems'])

    category = ['ID', 'Total']
    data_graph = {'ID': number_id, 'Total': number_total}
    x = [(system, cat) for system in system_names for cat in category]
    count = sum(zip(data_graph['ID'], data_graph['Total']), ())
    source = ColumnDataSource(data=dict(x=x, counts=count))

    output_file(Path.cwd() / output / name / "hbar_ID_Total_{0}.html".format(name))

    # Set size for plot
    len_sys = len(system_names)
    max_value = max(number_total)
    width = max_value * 5
    height = len_sys * 35

    if width < 500:
        width = 500

    p = figure(y_range=FactorRange(*x), title="Systems count for {0}".format(name),
               height=height, width=width)

    p.hbar(y='x', right='counts', source=source, height=0.8, fill_color="#458B74", line_color="black")

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.x_range.range_padding = 0.01
    p.yaxis.major_label_orientation = 0
    p.yaxis.group_label_orientation = 0
    p.ygrid.grid_line_color = None

    save(p)
