#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import List, Set
from pathlib import Path

# installed libraries
from tqdm import tqdm
import numpy as np
import pandas as pd
from bokeh.models import (
    BasicTicker,
    ColorBar,
    LinearColorMapper,
    CategoricalColorMapper,
    PrintfTickFormatter,
)
from bokeh.palettes import Magma256
from bokeh.io import output_file, export_png
from bokeh.plotting import figure, save


def systems_partition(
    name: str, system_projection: pd.DataFrame, output: Path, disable_bar: bool = False
) -> None:
    """
    Call functions to draw heatmap figures for the pangenome.

    Args:
        name (str): Name of the pangenome.
        system_projection (pd.DataFrame): Systems in the pangenome.
        output (Path): Path to the output directory.
        disable_bar (bool, optional): If True, disable the progress bar. Defaults to False.
    """
    data_systems = system_projection[
        ["system number", "system name", "organism", "partition"]
    ]
    system_names = data_systems["system name"].unique().tolist()
    system_names.sort(key=str.casefold)
    organism_names = data_systems["organism"].unique().tolist()
    organism_names.sort(key=str.casefold)

    matrix_genomes_systems = np.zeros((len(organism_names), len(system_names)))
    for i, organism in tqdm(
        enumerate(organism_names),
        unit="organisms",
        desc="Write system partition",
        disable=disable_bar,
    ):
        data_genome = data_systems[data_systems["organism"] == organism]
        dict_defense_genome = pd.Series(
            data_genome["system name"].values, index=data_genome["system number"]
        ).to_dict()
        for j, system in enumerate(system_names):
            matrix_genomes_systems[i][j] = sum(
                x == system for x in dict_defense_genome.values()
            )

    figure_partition_heatmap(
        name=name,
        data=data_systems,
        list_systems=system_names,
        list_organisms=organism_names,
        output=output,
    )
    figure_count_heatmap(
        name=name,
        data=matrix_genomes_systems,
        list_system=system_names,
        list_organism=organism_names,
        output=output,
    )


def figure_partition_heatmap(
    name: str,
    data: pd.DataFrame,
    list_systems: List[str],
    list_organisms: List[str],
    output: Path,
    out_format: List[str] = None,
) -> None:
    """
    Draw a partition heatmap figure for the pangenome.

    Args:
        name (str): Name of the pangenome.
        data (pd.DataFrame): Data used to produce the heatmap.
        list_systems (List[str]): List of systems in the pangenome.
        list_organisms (List[str]): List of organisms in the pangenome.
        output (Path): Path to the output directory.
        out_format (List[str], optional): Format of the output file. Defaults to None.
    """

    def conciliate_system_partition(system_partition: Set[str]) -> str:
        """
        Conciliate the partition of the system.

        Args:
            system_partition (Set[str]): All found partitions for gene that code the system.

        Returns:
            str: Reconciled system partition.

        Raises:
            Exception: If no partition is found. Could happen if a partition is not loaded or computed.
        """
        if len(system_partition) == 1:
            return system_partition.pop()
        else:
            if "persistent" in system_partition:
                return "persistent|accessory"
            else:
                return "accessory"

    out_format = out_format if out_format is not None else ["html"]
    data_partition = data.pivot_table(
        index="organism",
        columns="system name",
        values="partition",
        fill_value="Not_found",
        aggfunc=lambda x: conciliate_system_partition(set(x)),
    )
    df_stack_partition = pd.DataFrame(
        data_partition.stack(), columns=["partition"]
    ).reset_index()

    colors = ["#79DEFF", "#00D860", "#EB37ED", "#F7A507", "#FF2828"]
    partitions = ["cloud", "shell", "accessory", "persistent", "persistent|accessory"]
    mapper = CategoricalColorMapper(
        palette=colors, factors=partitions, nan_color="white"
    )
    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"
    p = figure(
        title="Partition of systems for {0}".format(name),
        x_range=list_systems,
        y_range=list_organisms,
        x_axis_location="above",
        tools=tools,
        toolbar_location="below",
        tooltips=[
            ("Partition", "@partition"),
            ("Organism", "@organism"),
            ("System", "@{system name}"),
        ],
        width=1200,
        height=900,
    )
    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = "Systems name"
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = "Genomes name"
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.xgrid.grid_line_color = None

    p.rect(
        x="system name",
        y="organism",
        width=1,
        height=1,
        source=df_stack_partition,
        fill_color={"field": "partition", "transform": mapper},
        line_color=None,
    )

    color_bar = ColorBar(
        color_mapper=mapper,
        major_label_text_font_size="14px",
        label_standoff=6,
        border_line_color=None,
        major_label_overrides={
            0: "Persistent",
            1: "Shell",
            2: "Cloud",
            3: "Accessory",
            4: "Persistent|accessory",
        },
        major_tick_line_color="black",
        bar_line_color="black",
        bar_line_width=0.2,
        border_line_width=0.2,
    )
    p.add_layout(color_bar, "right")
    if "html" in out_format:
        output_path = output / "partition.html"
        output_file(output_path)
        save(p)
        logging.getLogger("PANORAMA").debug(
            f"Saved partition heatmap in HTML format to {output_path}"
        )
    if "png" in out_format:
        output_path = output / "partition.png"
        output_file(output_path)
        export_png(p, width=1920, height=1080)
        logging.getLogger("PANORAMA").debug(
            f"Saved partition heatmap in PNG format to {output_path}"
        )


def figure_count_heatmap(
    name: str,
    data: np.ndarray,
    list_system: List[str],
    list_organism: List[str],
    output: Path,
) -> None:
    """
    Draw count heatmap figure for the pangenome.

    Args:
        name (str): Name of the pangenome.
        data (np.ndarray): Data used to produce the heatmap.
        list_system (List[str]): List of systems in the pangenome.
        list_organism (List[str]): List of organisms in the pangenome.
        output (Path): Path to the output directory.
    """
    output_path = output / "count.html"
    output_file(output_path)
    df_gcount = pd.DataFrame(data, index=list_organism, columns=list_system)
    df_stack_gcount = pd.DataFrame(df_gcount.stack(), columns=["number"]).reset_index()

    mapper = LinearColorMapper(
        palette=list(reversed(Magma256)),
        low=1,
        low_color="#FFFFFF",
        high=df_stack_gcount.number.max(),
    )
    tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

    p = figure(
        title="Systems for {0}".format(name),
        x_range=list_system,
        y_range=list_organism,
        x_axis_location="above",
        tools=tools,
        toolbar_location="below",
        tooltips=[
            ("Count", "@number"),
            ("Organism", "@level_0"),
            ("System", "@level_1"),
        ],
        width=1200,
        height=900,
    )

    p.title.align = "center"
    p.title.text_font_size = "20pt"
    p.xaxis.axis_label = "Systems name"
    p.xaxis.major_label_orientation = 1
    p.yaxis.axis_label = "Genomes name"
    p.axis.axis_label_text_font_size = "16pt"
    p.axis.axis_line_color = None
    p.axis.major_label_text_font_size = "12px"
    p.axis.major_label_standoff = 2
    p.xgrid.grid_line_color = None

    p.rect(
        x="level_1",
        y="level_0",
        width=1,
        height=1,
        source=df_stack_gcount,
        fill_color={"field": "number", "transform": mapper},
        line_color=None,
    )

    color_bar = ColorBar(
        color_mapper=mapper,
        major_label_text_font_size="14px",
        ticker=BasicTicker(desired_num_ticks=int(data.max())),
        formatter=PrintfTickFormatter(format="%s"),
        label_standoff=6,
        border_line_color=None,
    )
    p.add_layout(color_bar, "right")
    save(p)
    logging.getLogger("PANORAMA").debug(
        f"Saved count heatmap in HTML format to {output_path}"
    )
