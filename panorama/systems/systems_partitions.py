#!/usr/bin/env python3
# coding:utf-8
"""
Systems partitions visualization module for pangenome analysis.

This module provides functionality to create heatmap visualizations for
pangenome systems partitions and system counts across organisms.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Set

import numpy as np
import pandas as pd
from bokeh.models import (
    BasicTicker,
    CategoricalColorMapper,
    ColorBar,
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    PrintfTickFormatter,
    FactorRange,
)
from bokeh.layouts import gridplot, row
from bokeh.transform import linear_cmap, transform
from bokeh.palettes import Magma256
from bokeh.plotting import figure, show
from tqdm import tqdm

from panorama.systems.utils import VisualizationBuilder

# Configuration constants
DEFAULT_COLORS = ["#79DEFF", "#00D860", "#EB37ED", "#F7A507", "#FF2828"]
PARTITION_TYPES = ["cloud", "shell", "accessory", "persistent", "persistent|accessory"]
DEFAULT_FIGURE_WIDTH = 1200
DEFAULT_FIGURE_HEIGHT = 900

# Logger instance
logger = logging.getLogger("PANORAMA")


class SystemsPartitionVisualizer(VisualizationBuilder):
    """
    A class to create heatmap visualizations for pangenome systems partitions.

    This class handles the generation of both partition and count heatmaps
    for pangenome analysis, providing methods to visualize how systems are
    distributed across different organisms.
    """

    def __init__(self, name: str, output_dir: Path, formats: List[str] = None):
        """
        Initialize the SystemsPartitionVisualizer.

        Args:
            name: Name of the pangenome for visualization titles.
            output_dir: Directory path where output files will be saved.
        """
        super().__init__(name, output_dir, formats)
        self.color_palette = ["#79DEFF", "#00D860", "#EB37ED", "#F7A507", "#FF2828"]
        self.partitions = [
            "cloud",
            "shell",
            "accessory",
            "persistent",
            "persistent|accessory",
        ]
        self.partition2color = dict(zip(self.partitions, self.color_palette))
        self.mapper = CategoricalColorMapper(
            palette=self.color_palette, factors=self.partitions, nan_color="white"
        )
        self.top_bar = None

    def create_left_bar(self, source: ColumnDataSource):
        # Create the plot
        self.left_bar = figure(
            y_range=source.data["systems"],
            width=VisualizationBuilder.LEFT_WIDTH,
            height=VisualizationBuilder.MIDDLE_HEIGHT,
            toolbar_location=None,
            tools="",
        )

        # Add bars for each partition
        self.left_bar.hbar_stack(
            self.partitions,
            y="systems",
            color=[self.partition2color[part] for part in self.partitions],
            source=source,
            height=0.9,
        )
        #
        # # Customize the plot
        self.left_bar.yaxis.visible = False
        self.left_bar.xaxis.axis_label = "Count"
        self.left_bar.x_range.flipped = True
        self.left_bar.grid.grid_line_color = None
        self.left_bar.outline_line_color = None
        #
        # # Add hover tool
        hover = HoverTool(
            tooltips=[
                ("System", "@systems"),
                ("Partition", "$name"),
                ("Count", "@$name{count}"),
            ],
        )
        self.left_bar.add_tools(hover)

        # show(self.left_bar)

    def create_top_bar_plot(self, source: ColumnDataSource, data: pd.DataFrame):
        x_order = data.sort_values("count", ascending=False)["organism"].tolist()
        self.top_bar = figure(
            x_range=x_order,
            height=VisualizationBuilder.TOP_HEIGHT,
            width=VisualizationBuilder.CENTER_WIDTH,
            toolbar_location=None,
            tools="",
        )

        # Add bars for each partition
        self.top_bar.vbar(
            x="organism",
            top="count",
            width=0.9,
            source=source,
            color="green",
            alpha=0.6,
        )
        #
        self.top_bar.add_tools(
            HoverTool(
                tooltips=[
                    ("organism", "@{organism}"),
                    ("count", "@count"),
                ],
            )
        )

        # Configure the top bar plot
        self.top_bar.xaxis.visible = False
        self.top_bar.yaxis.axis_label = "Count"
        self.top_bar.grid.grid_line_color = None
        self.top_bar.outline_line_color = None

        # show(self.top_bar)

    def create_bar_plot(self, partition_matrix: pd.DataFrame):
        filtered_partition_matrix = partition_matrix[
            partition_matrix["partition"] != "Not_found"
        ]
        pivot_df = (
            filtered_partition_matrix.groupby(["system name", "partition"])
            .size()
            .unstack(fill_value=0)
        )

        # Prepare data for stacked bars using dictionary comprehension
        source_data = {
            "systems": pivot_df.index.tolist(),
            **{
                partition: pivot_df[partition].tolist()
                for partition in self.partitions
                if partition in pivot_df.columns
            },
        }

        source = ColumnDataSource(source_data)
        self.create_left_bar(source)

        organism_count = pd.DataFrame(
            filtered_partition_matrix["organism"].value_counts()
        ).reset_index()
        source = ColumnDataSource(organism_count)
        self.create_top_bar_plot(source, organism_count)

    def create_color_bar(self, title: str):
        """
        Create a color bar for the correlation matrix.

        Args:


        Returns:
            Bokeh figure containing the color bar.
        """
        color_bar = ColorBar(
            color_mapper=self.mapper,
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

        self.color_bar = figure(
            title=title,
            title_location="right",
            height=VisualizationBuilder.MIDDLE_HEIGHT,
            width=VisualizationBuilder.RIGHT_WIDTH,
            toolbar_location=None,
            min_border=0,
            outline_line_color=None,
        )

        self.color_bar.add_layout(color_bar, "right")
        self.color_bar.title.align = "center"
        self.color_bar.title.text_font_size = "14pt"

    def create_main_figure(
        self,
        partition_matrix: pd.DataFrame,
        x_range: FactorRange,
        y_range: FactorRange,
    ):
        """
        Create the main correlation matrix heatmap figure.

        Args:
            partition_matrix: Preprocessed partition matrix.
            x_range: X-axis range for the plot.
            y_range: Y-axis range for the plot.

        Returns:
            Tuple of the Bokeh figure and GlyphRenderer objects.
        """
        tooltips = [
            ("Partition", "@partition"),
            ("Organism", "@organism"),
            ("System", "@{system name}"),
        ]
        self._create_main_figure(partition_matrix, x_range, y_range, tooltips)
        source = ColumnDataSource(partition_matrix.reset_index())

        self.glyph_renderer = self.main_plot.rect(
            "organism",
            "system name",
            1,
            1,
            source=source,
            line_color="white",
            fill_color={"field": "partition", "transform": self.mapper},
        )
        self.main_plot.xaxis.major_label_orientation = 0.5

    def plot(self):
        grid_layout_matrix = [
            [None, self.top_bar, None, None],
            [self.left_bar, self.main_plot, self.color_bar],
        ]

        grid_layout = gridplot(grid_layout_matrix, toolbar_location="above")

        self._save_figure(grid_layout, "partition")


def preprocess_data(data: pd.DataFrame, disable_bar: bool = False) -> pd.DataFrame:
    """
    Preprocess data to draw a partition heatmap figure for the pangenome.

    Args:
        data (pd.DataFrame): Projection of pangenome systems
        disable_bar (bool, optional): If True, disable the progress bar. Defaults to False.
    """
    system_names = data["system name"].unique().tolist()
    system_names.sort(key=str.casefold)
    organism_names = data["organism"].unique().tolist()
    organism_names.sort(key=str.casefold)

    matrix_genomes_systems = np.zeros((len(system_names), len(organism_names)))
    for j, organism in tqdm(
        enumerate(organism_names),
        unit="organisms",
        desc="Write system partition",
        disable=disable_bar,
    ):
        data_genome = data[data["organism"] == organism]
        dict_defense_genome = pd.Series(
            data_genome["system name"].values, index=data_genome["system number"]
        ).to_dict()
        for i, system in enumerate(system_names):
            matrix_genomes_systems[i][j] = sum(
                x == system for x in dict_defense_genome.values()
            )
    data_count = pd.DataFrame(
        matrix_genomes_systems, columns=organism_names, index=system_names
    )
    data_count.index.name = "system name"
    data_count.columns.name = "organism"
    return data_count


def preprocess_partition_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Preprocess data to draw a partition heatmap figure for the pangenome.

    Args:
        data (pd.DataFrame): Data used to produce the heatmap.
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

    return df_stack_partition


def systems_partition(
    name: str,
    system_projection: pd.DataFrame,
    output: Path,
    output_formats: List[str] = None,
    disable_bar: bool = False,
) -> None:
    """
    Create heatmap visualizations for pangenome systems partitions.

    This function serves as the main entry point for generating both partition
    and count heatmap visualizations from pangenome system projection data.

    Args:
        name: Name of the pangenome for visualization titles.
        system_projection: DataFrame containing system projection data with columns:
            ['system number', 'system name', 'organism', 'partition'].
        output: Path to the directory where output files will be saved.
        output_formats: List of output formats for visualizations.
            Valid options: ['html', 'png']. Defaults to ['html'].
        disable_bar: If True, disable the progress bar during processing.

    Raises:
        ValueError: If required columns are missing from system_projection DataFrame.
    """
    if output_formats is None:
        output_formats = [VisualizationBuilder.DEFAULT_FORMAT]

    # Validate output formats

    for fmt in output_formats:
        if fmt not in VisualizationBuilder.OUTPUT_FORMATS:
            raise ValueError(
                f"Unsupported output format: {fmt}. "
                f"Supported formats: {VisualizationBuilder.OUTPUT_FORMATS}"
            )

    # data_genomes_systems = preprocess_data(system_projection, disable_bar=disable_bar)
    partitions_matrix = preprocess_partition_data(system_projection)
    viz_builder = SystemsPartitionVisualizer(
        name=name, output_dir=output, formats=output_formats
    )
    viz_builder.create_bar_plot(partitions_matrix)
    viz_builder.create_main_figure(
        partitions_matrix, viz_builder.top_bar.x_range, viz_builder.left_bar.y_range
    )
    viz_builder.create_color_bar("Partition")
    viz_builder.plot()
