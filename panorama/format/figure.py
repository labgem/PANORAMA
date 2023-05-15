#!/usr/bin/env python3
# coding:utf-8

import numpy as np
from bokeh.models import (BasicTicker, ColorBar, LinearColorMapper, CategoricalColorMapper, PrintfTickFormatter, ColumnDataSource)
from bokeh.palettes import Magma256, magma
import pandas as pd
from bokeh.io import output_file
from bokeh.plotting import figure, save
from bokeh.layouts import gridplot
from bokeh.palettes import turbo
from pathlib import Path
from bokeh.models import Title, FactorRange

def per_pan_heatmap(pan_name: str, system_proj, output: str):
    data_systems = system_proj[['system number', 'system name', 'organism', 'partition']]
    system_names = list(data_systems['system name'].unique())  # We have a list with system names without duplicates
    system_names.sort(key=str.casefold)
    organism_names = list(data_systems['organism'].unique())  # We have a list with organism names without duplicates
    organism_names.sort(key=str.casefold)
    organism_names.sort(key=lambda x: ('assembled' in x or 'MAGS' in x, x), reverse=True)

    matrix_genomes_systems = np.zeros((len(organism_names), len(system_names)))
    for i, organism in enumerate(organism_names):
        data_genome = data_systems[(data_systems['organism'] == organism)]
        dict_defense_genome = pd.Series(data_genome['system name'].values, index=data_genome['system number']).to_dict()
        for j, system in enumerate(system_names):
            matrix_genomes_systems[i][j] = sum(x == system for x in dict_defense_genome.values())  # Sum of system number for each system

    figure_partition_heatmap(pan_name, data_systems, system_names, organism_names, output)
    figure_count_heatmap(pan_name, system_names, organism_names, matrix_genomes_systems, output)


def figure_partition_heatmap(pan_name: str, data, list_system, list_organism, output: str):

    output_file(Path.cwd() / output / pan_name / "{0}_partition.html".format(pan_name))
    data_partition = data.pivot_table(index='organism', columns='system name', values='partition',
                                              fill_value='Not_found', aggfunc=lambda x: ','.join(set(x)))
    df_stack_partition = pd.DataFrame(data_partition.stack(), columns=['partition']).reset_index()

    colors = ["#FF9933", "#3DCA2E", "#B1F1EB", "#FFFFFF"]
    partitions = ["persistent", "shell", "cloud", "Not_found"]
    mapper = CategoricalColorMapper(palette=colors, factors=partitions)
    Tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

    p = figure(title="Partition of defense systems for {0}".format(pan_name), x_range=list_system, y_range=list_organism,
               x_axis_location="above", sizing_mode='scale_both', tools=Tools, toolbar_location='below',
               tooltips=[('Partition', '@partition'), ('Organism', '@organism'), ('Defense system', '@{system name}')])

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Defense systems name'
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = 'Genomes name'
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.xgrid.grid_line_color = None

    p.rect(x="system name", y="organism", width=1, height=1, source=df_stack_partition,
           fill_color={'field': 'partition', 'transform': mapper}, line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="14px", label_standoff=6, border_line_color=None,
                         major_label_overrides={0: "Persistent", 1: "Shell", 2: "Cloud", 3: "Not_found"},
                         major_tick_line_color='black', bar_line_color='black', bar_line_width=0.2, border_line_width=0.2)

    p.add_layout(color_bar, 'right')
    save(p)

def figure_count_heatmap(pan_name: str, list_system, list_organism, matrix_genomes_systems, output: str):

    output_file(Path.cwd() / output / pan_name /"{0}_count.html".format(pan_name))
    df_gcount = pd.DataFrame(matrix_genomes_systems, index=list_organism, columns=list_system)
    df_stack_gcount = pd.DataFrame(df_gcount.stack(), columns=['number']).reset_index()

    mapper = LinearColorMapper(palette=list(reversed(Magma256)), low=1, low_color='#FFFFFF',
                               high=df_stack_gcount.number.max())
    Tools = "hover,save,pan,box_zoom,reset,wheel_zoom"
    p = figure(title="Defense systems for {0}".format(pan_name),
                x_range=list_system, y_range=list_organism,
                x_axis_location="above", sizing_mode='scale_both',
                tools=Tools, toolbar_location='below',
                tooltips=[('Count', '@number'), ('Organism', '@level_0'), ('Defense system', '@level_1')])

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Defense systems name'
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
                         ticker=BasicTicker(desired_num_ticks=int(matrix_genomes_systems.max())),
                         formatter=PrintfTickFormatter(format="%s"), label_standoff=6, border_line_color=None)
    p.add_layout(color_bar, 'right')
    save(p)

def heatmap(global_df, output):
    pan_names = list(global_df['pangenome name'].unique())
    pan_names.sort(key=str.casefold)
    system_names = list(global_df['system name'].unique())
    system_names.sort(key=str.casefold)

    save_sum_pans_systems = []
    save_sum_pans_ratio = []
    for i, pan in enumerate(pan_names):
        sum_pans_systems = {}
        sum_pans_ratio = {}
        pan_df = global_df[global_df['pangenome name'] == pan]
        number_org_tot = len(list(pan_df['organism'].unique()))
        for system in system_names:
            system_df = pan_df[pan_df['system name'] == system]
            number_org = len(list(system_df['organism'].unique()))
            sum_pans_systems[system] = number_org
            sum_pans_ratio[system] = number_org/number_org_tot

        save_sum_pans_systems.append(sum_pans_systems)
        save_sum_pans_ratio.append(sum_pans_ratio)

    system_names = sorted(list(system_names), key=str.casefold)
    df_pan_count = pd.DataFrame(save_sum_pans_systems, columns=system_names, index=pan_names).fillna(0)
    df_pan_ratio = pd.DataFrame(save_sum_pans_ratio, columns=system_names, index=pan_names).fillna(0)

    figure_histogram(pan_names, df_pan_count, system_names, output)
    figure_heatmap(pan_names, df_pan_count, df_pan_ratio, system_names, output)


def figure_histogram(pan_names, df_pan_count, system_names, output):

    output_file(Path.cwd() / output / "bar_stacked_pangenomes.html")

    dict_graph_pans = [(key, value) for i, (key, value) in enumerate(zip(pan_names, df_pan_count.values.tolist()))]
    dict_pan_names = {'systems_names' : system_names}
    dict_pan_names.update(dict_graph_pans)
    col_pan_names = ColumnDataSource(data=dict_pan_names)

    p = figure(title="Defense systems Count", background_fill_color="#FAFAFA", x_range=system_names,
               sizing_mode='scale_both', x_axis_location="below", tools="hover", toolbar_location=None,
               tooltips=[('Pangenome', '$name'),('Defense system', '@systems_names'), ('Count', '@$name')])

    p.vbar_stack(pan_names, x='systems_names', width=0.5, color=magma(len(pan_names)),
                 source=col_pan_names, legend_label=pan_names, line_color='black')

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Defense systems'
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

def figure_heatmap(pan_names, df_pan_count, df_pan_ratio, system_names, output):

    output_file(Path.cwd() / output / "pangenomes_count.html")
    df_stack_pan_ratio = pd.DataFrame(df_pan_ratio.stack(), columns=['ratio']).reset_index()
    mapper = LinearColorMapper(palette=(list(reversed(Magma256))), low=0.000001, low_color='#FFFFFF', high=df_stack_pan_ratio.ratio.max())

    Tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

    p = figure(title="Defense systems in pangenomes", x_range=system_names, y_range=pan_names,
               x_axis_location="above", sizing_mode='scale_both', tools=Tools, toolbar_location='below',
               tooltips=[('pangenome', '@level_0'), ('Defense system', '@level_1'), ('Ratio', '@ratio{0.00}'), ('Count', f"{df_pan_count}")])

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = 'Defense systems name'
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

def upsetplot(pan_name, pan_df, features_df, output):
    pan_df = pan_df[['system number', 'system name', 'organism']]
    system_name_pan = list(pan_df['system name'].unique())
    system_name_pan.sort(key=str.lower)
    organism_name_pan = list(pan_df['organism'].unique())
    organism_name_pan.sort(key=lambda x: ('assembled' in x or 'MAGS' in x, x), reverse=True)
    organism_id_pan = [i + 1 for i in range(len(organism_name_pan))]
    dict_organism_id_pan = {name: ID+1 for ID, name in enumerate(organism_name_pan)}
    features_df = features_df[['system_name', 'mod_organism', 'module', 'rgp_organism', 'spot', 'rgp']]
    #system_name_features = list(features_df['system_name'].unique())

    systems_count, organisms_count = count(pan_df, system_name_pan, organism_name_pan)

    #systems_count_features_list = [None] * (len(system_name_features) + 1)
    #count_system = 0

    spot_set = set()
    dict_spot_system = {}
    dict_spot_organism_id = {}

    module_set = set()
    dict_module_system = {}
    dict_module_organism_id = {}

    x_axis_systems_present = []
    y_axis_organisms_id_present = []


    for system in system_name_pan:
        pan_df_filter = pan_df[pan_df['system name'] == system]
        organism_id = 1
        for organism in dict_organism_id_pan.keys():
            checking_presence = (pan_df_filter['organism'].eq(organism)).to_string()
            if 'True' in checking_presence:
                x_axis_systems_present.append(system)
                y_axis_organisms_id_present.append(organism_id)
            organism_id += 1

        if features_df['system_name'].str.contains(system).any():
            data_features_filter = features_df[features_df['system_name'] == system]
            for org_index, row in data_features_filter.iterrows():
                rgp_organisms = row.rgp_organism
                mod_organisms = row.mod_organism
                spots = row.spot
                modules = row.module
                if spots is not None:
                    spot_set |= set(spots)
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

    width =len(system_name_pan) * 40
    height = len(organism_name_pan) * 40

    p_features = figure(width=width, height=height, x_range=system_name_pan, tooltips=tooltips,
                                            x_axis_location="above")
    p_features.square(x=x_axis_systems_present, y=y_axis_organisms_id_present, fill_color="white",
                                          line_color="black", legend_label="Alone", size=15, name="")

    if "no_spot" in spot_set:
        x_axis_systems_spot = dict_spot_system["no_spot"]
        y_axis_organisms_id_spot = dict_spot_organism_id["no_spot"]
        p_features.asterisk(x=x_axis_systems_spot, y=y_axis_organisms_id_spot, line_color="black",
                                                fill_color="white", size=20, legend_label="no_spot", name="no_spot")
        spot_set.remove("no_spot")

    spot_list = sorted(list(spot_set), key=float)
    color_spots = turbo(len(spot_list))
    for i, spot in enumerate(spot_list):
        x_axis_systems_spot = dict_spot_system[spot]
        y_axis_organisms_id_spot = dict_spot_organism_id[spot]
        p_features.square(x=x_axis_systems_spot, y=y_axis_organisms_id_spot, line_color="black",
                                                fill_color=color_spots[i], size=15,
                                                legend_label="Spot " + str(spot), name=str(spot))

    if "no_module" in module_set:
        x_axis_systems_module = dict_module_system["no_module"]
        y_axis_organisms_id_module = dict_module_organism_id["no_module"]
        p_features.asterisk(x=x_axis_systems_module, y=y_axis_organisms_id_module, line_color="black",
                                                fill_color="white", size=10, legend_label="no_module")
        module_set.remove("no_module")

    module_list = sorted(list(module_set), key=float)
    color_modules = turbo(len(module_list))
    for i, module in enumerate(module_list):
        x_axis_systems_module = dict_module_system[module]
        y_axis_organisms_id_module = dict_module_organism_id[module]
        p_features.circle(x=x_axis_systems_module, y=y_axis_organisms_id_module, line_color="black",
                                              fill_color=color_modules[i], size=10,
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
    p_features.add_layout(p_features.legend[0], 'right')
    p_features.legend.click_policy = "hide"
    p_features.outline_line_width = 1
    p_features.outline_line_color = "black"
    p_features.axis.major_label_text_font_size = "14px"
    p_features.yaxis.major_label_text_align = 'center'
    p_features.min_border_bottom = 100

    # Histogram of organism

    tooltips_org = [('Organism', '@y'), ('Count', '@left')]

    p_org = figure(x_range=[max(organisms_count) + 1, 0], y_range=organism_name_pan, width=width, height=height,
                   tooltips=tooltips_org)
    p_org.hbar(y=organism_name_pan, left=organisms_count, right=0, height=0.8, fill_color="#458B74",
               line_color="black")
    p_org.yaxis.visible = False
    p_org.outline_line_width = 1
    p_org.outline_line_color = "black"
    p_org.axis.major_label_text_font_size = "14px"
    p_org.xgrid.grid_line_color = "black"
    p_org.min_border_bottom = 100


    tooltips_sys = [('System name', '@x'), ('Count', '@top')]
    p_sys = figure(title="{0}".format(pan_name), x_range=system_name_pan, width=width, height=height, tooltips=tooltips_sys)
    p_sys.vbar(x=system_name_pan, top=systems_count, width=0.8, fill_color="#B22222", line_color="black")

    p_sys.title.align = "center"
    p_sys.title.text_font_size = "20pt"
    p_sys.outline_line_width = 1
    p_sys.outline_line_color = "black"
    p_sys.x_range.range_padding = 0.1
    p_sys.ygrid.grid_line_color = "black"
    p_sys.y_range.start = 0
    p_sys.xaxis.visible = False
    p_sys.axis.major_label_text_font_size = "14px"


    output_file(Path.cwd() / output / pan_name / f"upsetplot_{pan_name}.html")

    grid = gridplot([[None, p_sys], [p_org, p_features]])
    save(grid)

def count(pan_df, system_name_pan, organism_name_pan):

    systems_count = [None] * len(system_name_pan)
    for i, system in enumerate(system_name_pan):
        data_pangenome_filter = pan_df[pan_df['system name'] == system]
        data_pangenome_filter_drop = data_pangenome_filter.drop_duplicates(
            subset=['system number', 'system name', 'organism'])
        systems_count[i] = data_pangenome_filter_drop.shape[0]

    organisms_count = [None] * len(organism_name_pan)
    for j, organism in enumerate(organism_name_pan):
        data_pangenome_filter_2 = pan_df[pan_df['organism'] == organism]
        data_pangenome_filter_drop_2 = data_pangenome_filter_2.drop_duplicates(
            subset=['system number', 'system name', 'organism'])
        organisms_count[j] = data_pangenome_filter_drop_2.shape[0]
    return systems_count, organisms_count

def hbar_ID_type(pangenome_name, dataframe_ID, dataframe_ID_org, output):
    data_systems = dataframe_ID[['system name', 'Number of IDs']]
    system_names = list(data_systems['system name'])
    system_names.sort(key=str.casefold)
    number_ID = list(data_systems['Number of IDs'])

    data_systems_org = dataframe_ID_org[['system name', 'Number of IDs per org']]
    number_ID_org = list(data_systems_org['Number of IDs per org'])

    category = ['ID', 'Type']
    data_graph = {'ID' : number_ID, 'Type' : number_ID_org}
    x = [ (system, cat) for system in system_names for cat in category ]
    count = sum(zip(data_graph['ID'], data_graph['Type']), ())
    source = ColumnDataSource(data=dict(x=x, counts=count))

    output_file(Path.cwd() / output / pangenome_name / "hbar_ID_Type_{0}.html".format(pangenome_name))
    p = figure(y_range=FactorRange(*x), title="Systems count per ID or Type for {0}".format(pangenome_name),
               sizing_mode='scale_both')

    p.hbar(y='x', right='counts', source=source, height=0.8,
                   fill_color="#458B74", line_color="black")

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.x_range.range_padding = 0.01
    p.yaxis.major_label_orientation = 0
    p.yaxis.group_label_orientation = 0
    p.ygrid.grid_line_color = None

    save(p)
