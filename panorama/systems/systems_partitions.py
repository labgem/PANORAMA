#!/usr/bin/env python3
# coding:utf-8
"""
Systems partitions visualization module for pangenome analysis.

This module provides functionality to create heatmap visualizations for
pangenome systems partitions and system counts across organisms.
"""

import logging
from pathlib import Path
from typing import List, Optional, Set

import numpy as np
import pandas as pd
from bokeh.models import (
    CategoricalColorMapper,
    ColorBar,
    ColumnDataSource,
    HoverTool,
    FactorRange,
)
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from tqdm import tqdm

from panorama.systems.utils import VisualizationBuilder, conciliate_partition

# Configuration constants
DEFAULT_COLORS = ["#79DEFF", "#00D860", "#EB37ED", "#F7A507", "#FF2828"]
PARTITION_TYPES = ["cloud", "shell", "accessory", "persistent", "persistent|accessory"]
DEFAULT_FIGURE_WIDTH = 1200
DEFAULT_FIGURE_HEIGHT = 900

# Logger instance
logger = logging.getLogger("PANORAMA")


class SystemsPartitionVisualizer(VisualizationBuilder):
    """
    Visualizer for pangenome systems partition distributions.

    This class creates heatmap visualizations showing how pangenome systems are
    partitioned (persistent, shell, cloud, etc.) across different organisms.
    The visualization helps understand the conservation patterns of genetic systems.

    Attributes:
        partitions: List of partition categories
        partition2color: Mapping from partitions to colors
        mapper: Categorical color mapper for partitions
    """

    def __init__(self, name: str, output_dir: Path, formats: Optional[List[str]] = None):
        """
        Initialize the SystemsPartitionVisualizer.

        Args:
            name: Name of the pangenome for visualization titles
            output_dir: Directory path where output files will be saved
            formats: List of output formats to generate
        """
        super().__init__(name, output_dir, formats)

        # Partition configuration
        self.partitions = ["persistent|accessory", "persistent", "accessory", "shell", "cloud"]
        self.color_palette = ["#FF2828", "#F7A507", "#EB37ED", "#00D860", "#79DEFF"]
        self.partition2color = dict(zip(self.partitions, self.color_palette))

        # Color mapper for categorical partitions
        self.mapper = CategoricalColorMapper(
            palette=self.color_palette,
            factors=self.partitions,
            nan_color="white"
        )

    def create_left_bar(self, source: ColumnDataSource) -> None:
        """
        Create a stacked horizontal bar plot showing partition distributions by system.

        Args:
            source: ColumnDataSource containing stacked partition data
        """
        self.left_bar = figure(
            y_range=source.data["systems"],
            width=self.LEFT_WIDTH,
            height=self.MIDDLE_HEIGHT,
            toolbar_location=None,
            tools="",
        )

        # Create stacked bars for each partition
        self.left_bar.hbar_stack(
            self.partitions,
            y="systems",
            color=[self.partition2color[part] for part in self.partitions],
            source=source,
            height=0.9,
        )

        # Configure appearance
        self._configure_bar_plot_style(
            self.left_bar,
            x_label="Count",
            flip_x=True,
            hide_y_axis=True
        )

        # Add interactive hover tool
        hover = HoverTool(tooltips=[
            ("System", "@systems"),
            ("Partition", "$name"),
            ("Count", "@$name{0}"),
        ])
        self.left_bar.add_tools(hover)

    def create_color_bar(self, title: str) -> None:
        """
        Create a color bar for the partition matrix.

        Args:
            title: Title to display on the color bar
        """
        color_bar = ColorBar(
            color_mapper=self.mapper,
            major_label_text_font_size="14px",
            label_standoff=6,
            border_line_color=None,
            major_label_overrides={
                0: "Persistent",
                1: "Persistent|Accessory",
                2: "Accessory",
                3: "Shell",
                4: "Cloud",
            },
            major_tick_line_color="black",
            bar_line_color="black",
            bar_line_width=0.2,
            border_line_width=0.2,
        )

        self.color_bar = figure(
            title=title,
            title_location="right",
            height=self.MIDDLE_HEIGHT,
            width=self.RIGHT_WIDTH,
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
    ) -> None:
        """
        Create the main partition matrix heatmap figure.

        Args:
            partition_matrix: DataFrame with systems, organisms, and partition information
            x_range: X-axis range for the plot (organisms)
            y_range: Y-axis range for the plot (systems)
        """
        # Define tooltips for partition exploration
        tooltips = [
            ("Partition", "@partition"),
            ("Organism", "@organism"),
            ("System", "@{system name}"),
        ]

        # Create the base figure
        self._create_main_figure(partition_matrix, x_range, y_range, tooltips)

        # Prepare data source
        source = ColumnDataSource(partition_matrix.reset_index())

        # Create partition rectangles with categorical color mapping
        self.glyph_renderer = self.main_plot.rect(
            "organism",
            "system name",
            1,  # width
            1,  # height
            source=source,
            line_color="white",
            fill_color={"field": "partition", "transform": self.mapper},
        )

        # Configure x-axis label orientation for organism names
        self.main_plot.xaxis.major_label_orientation = 0.5
        self.main_plot.xaxis.axis_label = "Organism"

    def create_bar_plots(self, partition_matrix: pd.DataFrame) -> None:
        """
        Create bar plots showing system and organism statistics.

        Creates stacked bars for partition distributions (left) and organism
        counts (top) to provide marginal summaries of the partition matrix.

        Args:
            partition_matrix: DataFrame with partition information
        """
        # Filter out 'Not_found' entries for cleaner visualization
        filtered_matrix = partition_matrix[
            partition_matrix["partition"] != "Not_found"
            ]

        # Left bar plot: Partition distribution by system (stacked bars)
        pivot_df = (
            filtered_matrix.groupby(["system name", "partition"])
            .size()
            .unstack(fill_value=0)
        )

        # Prepare data for stacked bars - ensure all partitions are represented
        source_data = {"systems": pivot_df.index.tolist()}
        for partition in self.partitions:
            if partition in pivot_df.columns:
                source_data[partition] = pivot_df[partition].tolist()
            else:
                # Add zeros if partition not present in data
                source_data[partition] = [0] * len(pivot_df.index)

        left_bar_source = ColumnDataSource(source_data)
        self.create_left_bar(left_bar_source)

        # Top bar plot: Organism counts
        organism_counts = pd.DataFrame(
            filtered_matrix["organism"].value_counts()
        ).reset_index()
        organism_counts.columns = ["organism", "count"]

        top_bar_source = ColumnDataSource(organism_counts)

        # Sort organisms by count for better visualization
        x_order = organism_counts.sort_values("count", ascending=False)["organism"].tolist()

        self.create_top_bar_plot(
            top_bar_source,
            "organism",
            x_order=x_order,
            color="green"
        )

    def plot(self) -> None:
        """
        Create and save the complete partition visualization layout.

        Arranges the main partition heatmap with supporting bar plots and color bar
        in a clean grid layout, then saves the result in specified formats.
        """
        # Create grid layout with top bar, main components
        grid_layout_matrix = [
            [None, self.top_bar, None, None],
            [self.left_bar, self.main_plot, self.color_bar],
        ]

        # Create the final layout
        grid_layout = gridplot(grid_layout_matrix, toolbar_location="above")

        # Save the visualization
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
    data_partition = data.pivot_table(
        index="organism",
        columns="system name",
        values="partition",
        fill_value="Not_found",
        aggfunc=lambda x: conciliate_partition(set(x)),
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
    partitions_matrix = preprocess_partition_data(system_projection)
    viz_builder = SystemsPartitionVisualizer(
        name=name, output_dir=output, formats=output_formats
    )
    viz_builder.create_bar_plots(partitions_matrix)
    viz_builder.create_main_figure(
        partitions_matrix, viz_builder.top_bar.x_range, viz_builder.left_bar.y_range
    )
    viz_builder.create_color_bar("Partition")
    viz_builder.plot()
