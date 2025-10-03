#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functionality for associating systems with other pangenome elements
such as RGPs (Regions of Genomic Plasticity), spots, and modules.

The module creates correlation matrices and visualizations to analyze the relationships
between systems and various pangenome components.
"""

# default libraries
from __future__ import annotations
import logging
from collections import defaultdict, namedtuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Set, Tuple, Union, Optional
import time

# installed libraries
from tqdm import tqdm
import pandas as pd
from bokeh.layouts import gridplot, row
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.palettes import Colorblind, Reds256, Blues256, linear_palette
from bokeh.models import (
    BasicTicker,
    ColumnDataSource,
    LinearColorMapper,
    ColorBar,
    FactorRange,
    HoverTool,
)
from ppanggolin.region import Region

# local libraries
from panorama.pangenomes import Pangenome
from panorama.region import Spot, Module
from panorama.systems.system import System
from panorama.systems.utils import VisualizationBuilder

# Constants for association types
ASSOCIATION_RGPS = "RGPs"
ASSOCIATION_SPOTS = "spots"
ASSOCIATION_MODULES = "modules"

# All valid association types
VALID_ASSOCIATIONS = [ASSOCIATION_RGPS, ASSOCIATION_SPOTS, ASSOCIATION_MODULES]


class AssociationVisualizationBuilder(VisualizationBuilder):
    """
    Builder for correlation matrix visualizations between systems and genomic associations.

    This class creates comprehensive visualizations showing correlations between
    pangenome systems and various genomic elements (RGPs, modules, etc.), including
    coverage and frequency plots.

    Attributes:
        association: Type of genomic association being visualized (e.g., 'rgp', 'module')
    """

    def __init__(
            self,
            association: str,
            name: str,
            output_dir: Path,
            formats: Optional[List[str]] = None
    ):
        """
        Initialize the association visualization builder.

        Args:
            association: Type of pangenome object being visualized (e.g., 'rgp', 'module')
            name: Name of the pangenome for visualization titles
            output_dir: Directory path where output files will be saved
            formats: List of output formats to generate
        """
        super().__init__(name, output_dir, formats)
        self.association = association

        # Additional plots specific to association visualizations
        self.coverage_plot: Optional[figure] = None
        self.coverage_color_bar: Optional[figure] = None
        self.frequency_plot: Optional[figure] = None
        self.frequency_color_bar: Optional[figure] = None

    @staticmethod
    def create_color_palette(max_value: int) -> List[str]:
        """
        Create an appropriate color palette based on the maximum correlation value.

        The palette selection adapts to the data range to provide optimal visual
        discrimination between different correlation values.

        Args:
            max_value: Maximum correlation value in the matrix

        Returns:
            List of color hex codes for the palette, starting with white for zero values
        """
        if max_value == 1:
            # Binary data: white to black
            return ["#ffffff", "#000000"]
        elif max_value == 2:
            # Three categories: white + two colors from colorblind-safe palette
            return ["#ffffff"] + list(reversed(Colorblind[3]))[1:]
        elif max_value <= 8:
            # Small ranges: use discrete colorblind palette
            return ["#ffffff"] + list(Colorblind[max_value])
        else:
            # Large ranges: use continuous red palette
            n_colors = min(max_value + 4, 256)
            return ["#ffffff"] + list(linear_palette(Reds256[::-1], n_colors))[4:]

    def _configure_plot_style(self) -> None:
        """
        Configure plot styling specific to association visualizations.

        Extends the base styling with association-specific axis labels.
        """
        super()._configure_plot_style()

        # Set x-axis label based on association type
        if self.association.lower() == "rgps":
            x_label = "RGP name"
        else:
            x_label = f"{self.association} ID"

        self.main_plot.xaxis.axis_label = x_label
        self.main_plot.xaxis.major_label_orientation = 1

    def create_color_bar(self, title: str) -> None:
        """
        Create a color bar for the correlation matrix.

        Args:
            title: Title to display on the color bar
        """
        # Extract color mapper from the main plot's glyph
        color_mapper = self.glyph.fill_color["transform"]

        color_bar = ColorBar(
            color_mapper=color_mapper,
            label_standoff=12,
            ticker=BasicTicker(desired_num_ticks=len(color_mapper.palette)),
            border_line_color=None,
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
            correlation_matrix: pd.DataFrame,
            x_range: FactorRange,
            y_range: FactorRange,
    ) -> None:
        """
        Create the main correlation matrix heatmap figure.

        Args:
            correlation_matrix: Preprocessed correlation matrix with systems as rows
                               and associations as columns
            x_range: X-axis range for consistent ordering across plots
            y_range: Y-axis range for consistent ordering across plots
        """
        # Determine color palette based on data range
        max_correlation = correlation_matrix.values.max()
        color_palette = self.create_color_palette(max_correlation)

        # Define tooltips for interactive exploration
        tooltips = [
            (self.association.title(), f"@{self.association}"),
            ("System", "@system_name"),
            ("Count", "@corr"),
        ]

        # Create the base figure
        self._create_main_figure(correlation_matrix, x_range, y_range, tooltips)

        # Prepare data source from correlation matrix
        correlation_data = correlation_matrix.stack().reset_index(name="corr")
        source = ColumnDataSource(correlation_data)

        # Create the heatmap rectangles with color mapping
        self.glyph_renderer = self.main_plot.rect(
            self.association,
            "system_name",
            1,  # width
            1,  # height
            source=source,
            line_color="white",
            fill_color=linear_cmap(
                "corr",
                palette=color_palette,
                low=0,
                high=max_correlation + 1
            ),
        )

    def create_coverage_plot(self, coverage_df: pd.DataFrame, x_range: FactorRange) -> None:
        """
        Create a coverage visualization plot.

        Coverage represents how well each genomic element is covered by the systems,
        displayed as a horizontal strip below the main heatmap.

        Args:
            coverage_df: DataFrame containing coverage data with coverage values
            x_range: X-axis range for consistent ordering with main plot
        """
        self.coverage_plot, self.coverage_color_bar = self._create_metric_plot(
            coverage_df, x_range, "coverage", Reds256, "Coverage"
        )

    def create_frequency_plot(self, frequency_df: pd.DataFrame, x_range: FactorRange) -> None:
        """
        Create a frequency visualization plot.

        Frequency represents how often each genomic element appears across genomes,
        displayed as a horizontal strip below the main heatmap.

        Args:
            frequency_df: DataFrame containing frequency data with frequency values
            x_range: X-axis range for consistent ordering with main plot
        """
        self.frequency_plot, self.frequency_color_bar = self._create_metric_plot(
            frequency_df, x_range, "frequency", Blues256, "Genome frequencies"
        )

    def _create_metric_plot(
            self,
            data_df: pd.DataFrame,
            x_range: FactorRange,
            metric_name: str,
            color_palette: List[str],
            title: str,
    ) -> Tuple[figure, figure]:
        """
        Create a generic metric visualization plot (coverage or frequency).

        Creates a horizontal strip visualization with an associated color bar
        to show metric values across genomic elements.

        Args:
            data_df: DataFrame containing the metric data
            x_range: X-axis range for consistent ordering
            metric_name: Name of the metric column in the DataFrame
            color_palette: Color palette to use for the visualization
            title: Title for the color bar

        Returns:
            Tuple of (metric_plot, color_bar_plot)
        """
        # Create color mapper (inverted scale: high values = dark colors)
        mapper = LinearColorMapper(palette=color_palette, low=1, high=0)

        # Create color bar
        color_bar = ColorBar(
            color_mapper=mapper,
            label_standoff=12,
            ticker=BasicTicker(desired_num_ticks=5),
            border_line_color=None,
            location=(0, 0),
        )

        # Create color bar plot
        color_bar_plot = figure(
            title=title,
            title_location="below",
            height=int(self.BELOW_HEIGHT / 3),
            width=int(self.CENTER_WIDTH / 2.5),
            toolbar_location=None,
            min_border=0,
            outline_line_color=None,
        )
        color_bar_plot.add_layout(color_bar, "below")
        color_bar_plot.title.align = "center"
        color_bar_plot.title.text_font_size = "12pt"

        # Prepare data source aligned with x_range
        aligned_data = pd.DataFrame({
            self.association: x_range.factors,
            metric_name: data_df.loc[x_range.factors][metric_name],
        })
        data_source = ColumnDataSource(aligned_data)

        # Create main metric plot
        metric_plot = figure(
            x_range=x_range,
            height=int(self.BELOW_HEIGHT / 4),
            width=self.CENTER_WIDTH,
            tooltips=[
                (self.association.title(), f"@{self.association}"),
                (metric_name.title(), f"@{metric_name}"),
            ],
        )

        # Configure minimal appearance
        self._configure_minimal_plot(metric_plot)

        # Add rectangles for metric visualization
        metric_plot.rect(
            self.association,
            0.5,  # y-position (centered)
            1,    # width
            1,    # height
            source=data_source,
            fill_color={"field": metric_name, "transform": mapper},
        )

        return metric_plot, color_bar_plot

    def create_bar_plots(self, correlation_matrix: pd.DataFrame) -> None:
        """
        Create bar plots showing system and element counts.

        Creates both left (system counts) and top (element counts) bar plots
        to provide marginal summaries of the correlation matrix.

        Args:
            correlation_matrix: Preprocessed correlation matrix
        """
        # Left bar plot: System counts (sum across rows)
        system_counts = correlation_matrix.sum(axis=1)
        left_bar_data = pd.DataFrame({
            "system_name": list(correlation_matrix.index),
            "count": system_counts.to_list(),
        })
        left_bar_source = ColumnDataSource(left_bar_data)
        self.create_left_bar_plot(left_bar_source, correlation_matrix)

        # Top bar plot: Element counts (sum across columns)
        element_counts = correlation_matrix.sum(axis=0)

        # Sort elements by count for better visualization
        x_order = sorted(
            correlation_matrix.columns,
            key=lambda x: element_counts.loc[x],
            reverse=True,
        )

        top_bar_data = pd.DataFrame({
            self.association: list(correlation_matrix.columns),
            "count": element_counts.to_list(),
        })
        top_bar_source = ColumnDataSource(top_bar_data)
        self.create_top_bar_plot(
            top_bar_source,
            self.association,
            x_order=x_order,
            color="green"
        )

    def plot(self) -> None:
        """
        Create and save the complete association visualization layout.

        Arranges all components (main heatmap, bar plots, color bars, and metric plots)
        in a grid layout and saves the result in the specified formats.
        """
        # Base layout: top bar, main components
        grid_layout_matrix = [
            [None, self.top_bar, None, None],
            [self.left_bar, self.main_plot, self.color_bar],
            [None, self.coverage_plot, None],
        ]

        # Add frequency plot if it exists
        if self.frequency_plot is not None:
            grid_layout_matrix.insert(2, [None, self.frequency_plot, None])

            # Add color bars for both metrics
            color_bars_row = row(
                [self.frequency_color_bar, self.coverage_color_bar],
                align="end",
                spacing=100,
            )
            grid_layout_matrix.append([None, color_bars_row, None])
        else:
            # Add color bar for coverage only
            coverage_color_bar_row = row(self.coverage_color_bar, align="center")
            grid_layout_matrix.append([None, coverage_color_bar_row, None])

        # Create the final layout
        grid_layout = gridplot(grid_layout_matrix, toolbar_location="above")

        # Save the visualization
        self._save_figure(grid_layout, f"correlation_{self.association}")


def _get_region_frequency(region: Region, pangenome: Pangenome) -> float:
    """
    Calculate the frequency of a Region element across all organisms.

    Args:
        region: The Region element for frequency calculation.
        pangenome: The pangenome containing organism information.

    Returns:
        The frequency of the Region in the pangenome.

    Note:
        TODO: This implementation needs to be fixed to properly calculate
        region frequency across organisms.
    """
    return 1.0 / pangenome.number_of_organisms


def _get_element_frequency(
    element: Union[Spot, Module], system_organisms: Set, pangenome: Pangenome
) -> float:
    """
    Calculate the frequency of a Spot or Module element across organisms.

    Args:
        element: The Spot or Module element.
        system_organisms: Set of organisms associated with systems.
        pangenome: The pangenome containing organism information.

    Returns:
        The frequency of the element among organisms.
    """
    element_organisms = set(element.organisms)
    intersection = element_organisms.intersection(system_organisms)
    return len(intersection) / pangenome.number_of_organisms


def create_coverage_dataframe(
    element_to_systems: Dict[Union[Region, Spot, Module], Set[System]],
    pangenome: Optional[Pangenome] = None,
) -> pd.DataFrame:
    """
    Create a DataFrame describing coverage of systems by pangenome elements.

    Args:
        element_to_systems: Dictionary mapping pangenome elements to system sets.
        pangenome: Optional pangenome object for frequency calculations.

    Returns:
        DataFrame with coverage and frequency information.
    """
    fields = ["name", "systems_ID", "systems_name", "coverage"]
    if pangenome is not None:
        fields.append("frequency")

    ElementRecord = namedtuple("ElementRecord", fields)
    records = []

    # Return empty DataFrame if no elements provided
    if not element_to_systems:
        return pd.DataFrame(columns=fields)

    # Determine frequency calculation method based on the element type
    first_element = next(iter(element_to_systems.keys()))
    is_region = isinstance(first_element, Region)

    for element, systems in element_to_systems.items():
        element_families = set(element.families)
        system_families, system_organisms, system_ids, system_names = (
            set(),
            set(),
            set(),
            set(),
        )

        # Aggregate data from all associated systems
        for system in systems:
            system_families.update(system.families)
            system_organisms.update(system.organisms)
            system_ids.add(system.ID)
            system_names.add(system.name)

        # Calculate coverage as intersection over union
        coverage = (
            len(element_families.intersection(system_families)) / len(element_families)
            if element_families
            else 0.0
        )

        record_data = [
            str(element),
            ",".join(system_ids),
            ",".join(system_names),
            coverage,
        ]

        # Add frequency if pangenome is available
        if pangenome is not None:
            if is_region:
                frequency = _get_region_frequency(element, pangenome)
            else:
                frequency = _get_element_frequency(element, system_organisms, pangenome)
            record_data.append(frequency)

        records.append(ElementRecord(*record_data))

    return pd.DataFrame(records)


def process_system(
    system: System,
    associations: List[str],
    rgp_to_systems: defaultdict,
    spot_to_systems: defaultdict,
    module_to_systems: defaultdict,
) -> Tuple[str, List[str]]:
    """
    Process a single system and update association mappings.

    Args:
        system: The system to process.
        associations: List of association types to include.
        rgp_to_systems: Mapping from RGPs to systems (updated in-place).
        spot_to_systems: Mapping from spots to systems (updated in-place).
        module_to_systems: Mapping from modules to systems (updated in-place).

    Returns:
        Tuple of system ID and system data list.
    """
    system_data = [system.name, ",".join(fam.name for fam in system.families)]

    if ASSOCIATION_RGPS in associations:
        rgp_names = {rgp.name for rgp in system.regions}
        for rgp in system.regions:
            rgp_to_systems[rgp].add(system)
        system_data.append(",".join(rgp_names))

    if ASSOCIATION_SPOTS in associations:
        spot_ids = {str(spot.ID) for spot in system.spots}
        for spot in system.spots:
            spot_to_systems[spot].add(system)
        system_data.append(",".join(spot_ids))

    if ASSOCIATION_MODULES in associations:
        module_ids = {str(mod.ID) for mod in system.modules}
        for mod in system.modules:
            module_to_systems[mod].add(system)
        system_data.append(",".join(module_ids))

    return system.ID, system_data


def get_association_dataframes(
    pangenome: Pangenome,
    associations: List[str],
    threads: int = 1,
    disable_progress_bar: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate DataFrames for system-pangenome element associations.

    Args:
        pangenome: Pangenome containing systems and elements.
        associations: List of pangenome elements to associate with systems.
        threads: Number of threads for parallel processing.
        disable_progress_bar: Whether to disable the progress bar.

    Returns:
        Tuple containing:
            - Association DataFrame (systems to elements)
            - RGP coverage DataFrame
            - Spot coverage DataFrame
            - Module coverage DataFrame

    Raises:
        ValueError: If no systems are found in the pangenome.
    """
    if pangenome.number_of_systems() == 0:
        raise ValueError("No systems found in the pangenome")

    # Define column structure based on requested associations
    columns = ["system_name", "families"]
    has_rgps = ASSOCIATION_RGPS in associations
    has_spots = ASSOCIATION_SPOTS in associations
    has_modules = ASSOCIATION_MODULES in associations

    if has_rgps:
        columns.append(ASSOCIATION_RGPS)
    if has_spots:
        columns.append(ASSOCIATION_SPOTS)
    if has_modules:
        columns.append(ASSOCIATION_MODULES)

    # Initialize data structures
    association_data = {}
    rgp_to_systems = defaultdict(set)
    spot_to_systems = defaultdict(set)
    module_to_systems = defaultdict(set)

    # Process systems in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_system,
                system,
                associations,
                rgp_to_systems,
                spot_to_systems,
                module_to_systems,
            )
            for system in pangenome.systems
        ]

        # Collect results with progress tracking
        for future in tqdm(
            as_completed(futures),
            total=pangenome.number_of_systems(),
            unit="systems",
            desc=f"Associating systems to: {', '.join(associations)}",
            disable=disable_progress_bar,
        ):
            system_id, system_data = future.result()
            association_data[system_id] = system_data

    # Create association DataFrame
    start_time = time.time()
    association_df = pd.DataFrame.from_dict(
        association_data, orient="index", columns=columns
    )
    association_df.index.name = "system_number"

    logging.getLogger("PANORAMA").debug(
        f"Association DataFrame created in {time.time() - start_time:.2f} seconds"
    )

    # Generate coverage DataFrames
    rgp_coverage_df = (
        create_coverage_dataframe(rgp_to_systems, pangenome)
        if has_rgps
        else pd.DataFrame()
    )
    spot_coverage_df = (
        create_coverage_dataframe(spot_to_systems, pangenome)
        if has_spots
        else pd.DataFrame()
    )
    module_coverage_df = (
        create_coverage_dataframe(module_to_systems, pangenome)
        if has_modules
        else pd.DataFrame()
    )

    return association_df, rgp_coverage_df, spot_coverage_df, module_coverage_df


def preprocess_association_data(
    dataframe: pd.DataFrame, association: str
) -> pd.DataFrame:
    """
    Preprocess association data to create a correlation matrix.

    Args:
        dataframe: Association DataFrame between systems and pangenome objects.
        association: Type of pangenome object for association.

    Returns:
        Preprocessed correlation matrix DataFrame.
    """
    # Split comma-separated associations into dummy variables
    processed_df = dataframe.drop(columns=["families"]).join(
        dataframe[association].str.get_dummies(sep=",")
    )
    processed_df = processed_df.drop(columns=[association])

    # Group by system name and sum associations
    correlation_matrix = processed_df.groupby("system_name").sum()
    correlation_matrix.sort_index(
        key=lambda x: x.str.lower(), ascending=False, inplace=True
    )
    correlation_matrix.columns.name = association

    return correlation_matrix


def write_correlation_matrix_visualization(
    association_df: pd.DataFrame,
    association: str,
    coverage_df: pd.DataFrame,
    pangenome_name: str,
    output_dir: Path,
    frequency_df: Optional[pd.DataFrame] = None,
    output_formats: Optional[List[str]] = None,
):
    """
    Generate and save correlation matrix visualization.

    Args:
        association_df: Association DataFrame between systems and pangenome objects.
        association: Type of pangenome object to visualize.
        coverage_df: Coverage DataFrame for the association.
        pangenome_name (str): Name of the pangenome.
        output_dir: Directory to save output files.
        frequency_df: Optional frequency DataFrame.
        output_formats: List of output formats (default: ['html']).

    Raises:
        ValueError: If an unsupported output format is specified.
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

    # Preprocess data for correlation matrix
    correlation_matrix = preprocess_association_data(association_df, association)

    # Create a visualization builder and plot components
    viz_builder = AssociationVisualizationBuilder(
        association, name=pangenome_name, output_dir=output_dir, formats=output_formats
    )
    viz_builder.create_bar_plots(correlation_matrix)

    viz_builder.create_main_figure(
        correlation_matrix, viz_builder.top_bar.x_range, viz_builder.left_bar.y_range
    )

    viz_builder.create_color_bar("# Systems")
    viz_builder.create_coverage_plot(coverage_df, viz_builder.top_bar.x_range)

    # Create a grid layout based on whether frequency data is available
    if frequency_df is not None:
        viz_builder.create_frequency_plot(frequency_df, viz_builder.top_bar.x_range)

    viz_builder.plot()


def create_pangenome_system_associations(
    pangenome: Pangenome,
    associations: List[str],
    output_dir: Path,
    output_formats: Optional[List[str]] = None,
    threads: int = 1,
    disable_bar: bool = False,
):
    """
    Create and save associations between systems and pangenome elements.

    This function generates association matrices, coverage analysis, and
    visualizations for the relationships between systems and various
    pangenome components (RGPs, spots, modules).

    Args:
        pangenome: The pangenome containing systems and other elements.
        associations: List of pangenome elements to associate with systems.
            Valid options: ['RGPs', 'spots', 'modules']
        output_dir: Directory where output files will be saved.
        output_formats: List of output formats for visualizations.
            Valid options: ['html', 'png']. Defaults to ['html'].
        threads: Number of threads for parallel processing. Defaults to 1.
        disable_bar: Whether to disable the progress bar display. Defaults to False.

    Raises:
        ValueError: If invalid association types are provided.
        FileNotFoundError: If the output directory doesn't exist.
    """
    if output_formats is None:
        output_formats = ["html"]

    # Validate associations
    for association in associations:
        if association not in VALID_ASSOCIATIONS:
            raise ValueError(
                f"Invalid association type: {association}. "
                f"Valid options: {VALID_ASSOCIATIONS}"
            )

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("PANORAMA")

    # Generate association DataFrames
    association_df, rgp_coverage_df, spot_coverage_df, module_coverage_df = (
        get_association_dataframes(pangenome, associations, threads, disable_bar)
    )

    # Save main association DataFrame
    association_output_path = output_dir / "association.tsv"
    association_df.to_csv(association_output_path, sep="\t")
    logger.info(f"Saved association DataFrame to {association_output_path}")

    # Process each association type
    for association in tqdm(
        associations,
        unit="association",
        desc="Creating system association visualizations",
        disable=disable_bar,
    ):
        coverage_df = None
        frequency_df = None
        should_create_visualization = False

        if association == ASSOCIATION_RGPS and not rgp_coverage_df.empty:
            coverage_df = rgp_coverage_df.set_index("name")
            # Save RGP-specific data
            coverage_output_path = output_dir / "rgp_to_systems.tsv"
            coverage_df.to_csv(coverage_output_path, sep="\t")
            logger.info(f"Saved RGP-to-systems DataFrame to {coverage_output_path}")
            should_create_visualization = True

        elif association == ASSOCIATION_SPOTS and not spot_coverage_df.empty:
            # Save spot-specific data
            spot_output_path = output_dir / "spot_to_systems.tsv"
            spot_coverage_df.set_index("name").to_csv(spot_output_path, sep="\t")
            logger.info(f"Saved spot-to-systems DataFrame to {spot_output_path}")

            # Prepare coverage and frequency data
            coverage_df = spot_coverage_df.set_index("name")
            coverage_df.index = coverage_df.index.str.replace("spot_", "")

            # Separate frequency and coverage columns
            frequency_df = coverage_df.loc[:, coverage_df.columns != "coverage"]
            coverage_df = coverage_df.loc[:, coverage_df.columns != "frequency"]
            should_create_visualization = True

        elif association == ASSOCIATION_MODULES and not module_coverage_df.empty:
            # Save module-specific data
            module_output_path = output_dir / "module_to_systems.tsv"
            module_coverage_df.set_index("name").to_csv(module_output_path, sep="\t")
            logger.info(f"Saved module-to-systems DataFrame to {module_output_path}")

            # Prepare coverage and frequency data
            coverage_df = module_coverage_df.set_index("name")
            coverage_df.index = coverage_df.index.str.replace("module_", "")

            # Separate frequency and coverage columns
            frequency_df = coverage_df.loc[:, coverage_df.columns != "coverage"]
            coverage_df = coverage_df.loc[:, coverage_df.columns != "frequency"]
            should_create_visualization = True

        # Create a visualization if data is available
        if should_create_visualization and coverage_df is not None:
            # Filter association DataFrame to only include the current association type
            filtered_association_df = association_df.drop(
                columns=[other for other in associations if other != association]
            )
            write_correlation_matrix_visualization(
                association_df=filtered_association_df,
                association=association,
                coverage_df=coverage_df,
                pangenome_name=pangenome.name,
                output_dir=output_dir,
                frequency_df=frequency_df,
                output_formats=output_formats,
            )
