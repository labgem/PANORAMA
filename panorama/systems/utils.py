#!/usr/bin/env python3

"""
This module provides utility functions to detect and write biological systems in pangenomes.
"""

# default libraries
from __future__ import annotations

import logging
import warnings
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

# installed libraries
import networkx as nx
import numpy as np
import pandas as pd
from bokeh.core.validation import silence
from bokeh.core.validation.warnings import MISSING_RENDERERS
from bokeh.io import export_png, output_file, save
from bokeh.models import (
    ColumnDataSource,
    FactorRange,
    GlyphRenderer,
    HoverTool,
)
from bokeh.models.glyph import Glyph
from bokeh.plotting import figure
from ppanggolin.genome import Organism

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome
from panorama.systems.models import Family, FuncUnit, Model

# Remove when pandas3.0 available. See the caveats in the documentation:
# https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
pd.options.mode.copy_on_write = True
# Silence-specific warning
silence(MISSING_RENDERERS, True)
# Silence output file will be overwritten info from bokeh
logging.getLogger("bokeh.io.state").setLevel(logging.WARNING)
# Silence-specific UserWarning while keeping others
warnings.filterwarnings(
    "ignore",
    message="Export method called with width or height argument on a non-Plot model.*",
    category=UserWarning,
)


def filter_global_context(graph: nx.Graph, jaccard_threshold: float = 0.8) -> nx.Graph[GeneFamily]:
    """
    Filters the edges of a gene family graph based on a Jaccard gene proportion threshold.

    Copies all nodes to a new graph and retains only those edges where both connected
    GeneFamily nodes have a Jaccard gene proportion (shared genomes over unique organisms)
    greater than or equal to the specified threshold. Updates edge data with Jaccard values
    and family names.

    Args:
        graph (nx.Graph):
            The input graph with GeneFamily nodes and edge data containing 'genomes'.
        jaccard_threshold (float, optional):
            Minimum Jaccard gene proportion required for both families to retain an edge.
            Defaults to 0.8.

    Returns:
        nx.Graph[GeneFamily]:
            A new graph with filtered edges and updated edge attributes.
    """
    new_graph = nx.Graph()
    # Copy all nodes from the original graph to the new graph
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        # Calculate Jaccard gene proportions for both families
        f1_proportion = len(data["genomes"]) / len(set(f1.organisms))
        f2_proportion = len(data["genomes"]) / len(set(f2.organisms))

        # Update the local copy of the edge data
        data.update(
            {
                "f1": f1.name,
                "f2": f2.name,
                "f1_jaccard_gene": f1_proportion,
                "f2_jaccard_gene": f2_proportion,
            }
        )

        # Add the edge to the new graph if it meets the Jaccard threshold
        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            # Add the edge with the updated edge data
            new_graph.add_edge(f1, f2, **data)

    return new_graph


def filter_local_context(
    graph: nx.Graph, organisms: Set[Organism], jaccard_threshold: float = 0.8
) -> nx.Graph[GeneFamily]:
    """
    Filters a graph based on a local Jaccard index.

    Args:
        graph (nx.Graph): A sub-pangenome graph.
        organisms (Set[Organism]):
            Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional):
            Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """
    orgs = frozenset(organisms)

    # Precompute family ∩ organisms intersections
    family_orgs: dict[GeneFamily, set[Organism]] = {fam: set(fam.organisms).intersection(orgs) for fam in graph.nodes}

    new_graph = nx.Graph()
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        genomes_in_orgs = data["genomes"].intersection(orgs)

        f1_orgs = family_orgs[f1]
        f2_orgs = family_orgs[f2]

        f1_proportion = len(genomes_in_orgs) / len(f1_orgs) if f1_orgs else 0.0
        f2_proportion = len(genomes_in_orgs) / len(f2_orgs) if f2_orgs else 0.0

        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            new_graph.add_edge(
                f1,
                f2,
                **data,
                f1=f1.name,
                f2=f2.name,
                f1_jaccard_gene=f1_proportion,
                f2_jaccard_gene=f2_proportion,
            )

    return new_graph


def filter_local_context_old(
    graph: nx.Graph, organisms: Set[Organism], jaccard_threshold: float = 0.8
) -> nx.Graph[GeneFamily]:
    """
    Filters a graph based on a local Jaccard index.

    Args:
        graph (nx.Graph): A sub-pangenome graph.
        organisms (Set[Organism]):
            Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional):
            Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """

    def get_gene_proportion(gene_family: GeneFamily) -> float:
        """Returns the Jaccard gene proportion for a given gene family."""
        # Compute the proportion if not cached
        gf_orgs = set(gene_family.organisms).intersection(organisms)
        if len(gf_orgs) > 0:
            return len(data["genomes"].intersection(organisms)) / len(gf_orgs)
        else:
            return 0

    new_graph = nx.Graph()
    # Copy all nodes from the original graph to the new graph
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        # Calculate Jaccard gene proportions for both families
        f1_proportion = get_gene_proportion(f1)
        f2_proportion = get_gene_proportion(f2)

        # Update the local copy of the edge data
        data.update(
            {
                "f1": f1.name,
                "f2": f2.name,
                "f1_jaccard_gene": f1_proportion,
                "f2_jaccard_gene": f2_proportion,
            }
        )

        # Add the edge to the new graph if it meets the Jaccard threshold
        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            # Add the edge with the updated edge data
            new_graph.add_edge(f1, f2, **data)

    return new_graph


def check_for_families(
    gene_families: Set[GeneFamily],
    gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
    mod_fam2meta_source: Dict[str, str],
    func_unit: FuncUnit,
) -> Tuple[bool, Dict[GeneFamily, Tuple[str, int]]]:
    """
    Evaluate gene families against a functional unit to detect forbidden, mandatory,
    and accessory family conditions.

    Args:
        gene_families (Set[GeneFamily]): Gene families to evaluate.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]]): Map from gene families to their model families.
        mod_fam2meta_source (Dict[str, str]): Map from model family name to metadata source.
        func_unit (FuncUnit): Functional unit definition to check against.

    Returns:
        Tuple[bool, Dict[GeneFamily, Tuple[str, int]]]:
            - A boolean indicating whether the conditions are satisfied
              (False immediately if a forbidden family is found).
            - A mapping from gene families to their selected metadata (source, meta_id).
    """

    # Ranking order for family presence
    presence_priority = {"forbidden": 0, "mandatory": 1, "accessory": 2, "neutral": 3}

    def fam_sort_key(item):
        """
        Determines the sorting key for an item based on specified criteria.

        Args:
            item (Tuple): A tuple consisting of a fam object and its associated metadata.

        Returns:
            Tuple: A tuple containing the rank based on presence, the negative score (for descending order),
            and the name of the fam to determine sorting precedence.
        """
        gfam, md_info = item
        score = md_info[2]
        return presence_priority.get(gfam.presence, 4), -score, gfam.name

    gf2meta_info: Dict[GeneFamily, Tuple[str, int]] = {}
    mandatory_seen, accessory_seen = set(), set()

    # Sort gene families by how many model families they map to (fewer first)
    for gf in sorted(gene_families, key=lambda n: len(gene_fam2mod_fam[n])):
        fam2meta_info = {}

        # Collect all metadata matches for this gene family
        for family in gene_fam2mod_fam[gf]:
            source = mod_fam2meta_source[family.name]
            available_names = {family.name, *family.exchangeable}

            for meta_id, metadata in gf.get_metadata_by_source(source).items():
                if metadata.protein_name in available_names or (
                    "secondary_name" in metadata.fields
                    and any(name in available_names for name in metadata.secondary_names.split(","))
                ):
                    fam2meta_info[family] = (source, meta_id, metadata.score)

        if not fam2meta_info:
            continue  # no metadata match, skip this GF

        # Pick the best candidate (highest priority, then score, then name)
        sorted_candidates = sorted(fam2meta_info.items(), key=fam_sort_key)
        family, best_meta_info = sorted_candidates[0]

        # Skip families already used if alternatives exist
        used_names = mandatory_seen | accessory_seen
        for fam, meta in sorted_candidates:
            if fam.name not in used_names:
                family, best_meta_info = fam, meta
                break

        # Forbidden condition → stop immediately
        if family.presence == "forbidden":
            return False, {}

        # Track mandatory / accessory usage
        if family.presence == "mandatory":
            mandatory_seen.add(family.name)
        elif family.presence == "accessory":
            accessory_seen.add(family.name)

        # Store only (source, meta_id), drop score
        gf2meta_info[gf] = best_meta_info[:2]

    # Check if constraints are satisfied
    if len(mandatory_seen) >= func_unit.min_mandatory and len(mandatory_seen | accessory_seen) >= func_unit.min_total:
        return True, gf2meta_info

    return False, {}


def get_gfs_matrix_combination(
    gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]]
) -> pd.DataFrame:
    """
    Build a matrix of association between gene families and families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]): Dictionary linking gene families to model families.

    Returns:
        pd.DataFrame: Matrix of association between gene families and families.
    """

    gfs = set()
    model_fams = set()
    for gf in gene_families:
        for fam in gene_fam2mod_fam[gf]:
            if fam.presence in ["mandatory", "accessory"]:
                gfs.add(gf)
                model_fams.add(fam.name)

    gfs = list(gfs)
    model_fams = list(model_fams)
    score_matrix = np.zeros((len(model_fams), len(gfs)), dtype=int)

    for j, gf in enumerate(gfs):
        annots = [f.name for f in gene_fam2mod_fam[gf]]
        for i, fam in enumerate(model_fams):
            if fam in annots:
                score_matrix[i, j] = 1

    return pd.DataFrame(score_matrix, index=model_fams, columns=[gf.name for gf in gfs])


def check_needed_families(matrix: pd.DataFrame, func_unit: FuncUnit) -> bool:
    """
    Check if there are enough mandatory and total families to satisfy the functional unit rules.

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        Boolean: True if satisfied, False otherwise

    Notes:
        This function assumes that a family could play multiple roles to satisfy
        the model requirements if it has multiple annotations
    """

    matrix = matrix.loc[matrix.sum(axis=1) > 0]  # Remove all-zero rows

    mandatory_fams = {fam.name for fam in func_unit.mandatory}
    mandatory_count = sum(fam in mandatory_fams for fam in matrix.index)

    return mandatory_count >= func_unit.min_mandatory and len(matrix) >= func_unit.min_total


def get_metadata_to_families(pangenome: Pangenome, sources: Iterable[str]) -> Dict[str, Dict[str, Set[GeneFamily]]]:
    """
    Retrieves a mapping of metadata to sets of gene families for each metadata source.

    Args:
        pangenome (Pangenome): Pangenome object containing gene families.
        sources (iterable of str): List of metadata source names.

    Returns:
        dict: A dictionary where each metadata source maps to another dictionary of metadata to sets of gene families.
    """
    meta2fam = {source: defaultdict(set) for source in sources}
    for source in sources:
        for gf in pangenome.gene_families:
            metadata = gf.get_metadata_by_source(source)
            if metadata is not None:
                for meta in metadata.values():
                    meta2fam[source][meta.protein_name].add(gf)
                    if "secondary_names" in meta.fields and meta.secondary_names != "":
                        for secondary_name in meta.secondary_names.split(","):
                            meta2fam[source][secondary_name].add(gf)
    return meta2fam


def dict_families_context(
    model: Model, annot2fam: Dict[str, Dict[str, Set[GeneFamily]]]
) -> Tuple[Dict[GeneFamily, Set[Family]], Dict[str, str]]:
    """
    Retrieves all gene families associated with the families in the model.

    Args:
        model (Model): Model containing the families.
        annot2fam (dict): Dictionary of annotated families.

    Returns:
        tuple: A tuple containing:
            - dict: Dictionary linking gene families to their families.
            - dict: Dictionary linking families to their sources.
    """
    gf2fam = defaultdict(set)
    fam2source = {}
    for fam_model in model.families:
        exchangeable = fam_model.exchangeable | {fam_model.name}
        for exch in exchangeable:
            for source, annotation2families in annot2fam.items():
                if exch in annotation2families:
                    if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                        logging.getLogger("PANORAMA").warning(
                            f"Protein annotation {fam_model.name} is encountered in multiple sources."
                            "All sources will be used, but only first one will be associated with "
                            "the model family."
                        )
                    else:
                        fam2source[fam_model.name] = source

                    for gf in annotation2families[exch]:
                        gf2fam[gf].add(fam_model)

    return gf2fam, fam2source


def conciliate_partition(partition: Set[str]) -> str:
    """
    Conciliate  a set of partition

    Args:
        partition (Set[str]): All partitions.

    Returns:
        str: The reconciled partition.
    """
    if len(partition) == 1:
        return partition.pop()
    else:
        if "persistent" in partition:
            return "persistent|accessory"
        else:
            return "accessory"


class VisualizationBuilder(ABC):
    """
    Abstract base class for building correlation matrix and partition visualizations.

    This class maintains a common configuration and provides a cohesive interface
    for creating all components of pangenome visualizations. It handles shared
    functionality like plot dimensions, styling, and file saving.
    """

    # Constants for plot dimensions
    TOTAL_WIDTH = 1780
    """int: Total width of the complete visualization layout.
    """
    TOTAL_HEIGHT = 920
    """int: Total height of the complete visualization layout.
    """
    # Layout proportions - these define how the total space is divided
    LEFT_WIDTH = int(0.15 * TOTAL_WIDTH)
    """int: Space for left bar plots
    """
    CENTER_WIDTH = int(0.75 * TOTAL_WIDTH)
    """int: Space for main heatmap
    """
    RIGHT_WIDTH = int(0.10 * TOTAL_WIDTH)
    """int: Space for color bars
    """
    TOP_HEIGHT = int(0.15 * TOTAL_HEIGHT)
    """int Space for top bar plots
    """
    MIDDLE_HEIGHT = int(0.70 * TOTAL_HEIGHT)
    """int: Space for main heatmap
    """
    BELOW_HEIGHT = int(0.15 * TOTAL_HEIGHT)
    """int: Space for bottom plots
    """
    # Output configuration
    OUTPUT_FORMATS = ["html", "png"]
    """list: Supported output formats for saving figures
    """
    DEFAULT_FORMAT = "html"
    """str: Default output format when none specified
    """

    # Constants for PNG export dimensions
    PNG_EXPORT_WIDTH = 1920
    """int: Width of PNG exports
    """
    PNG_EXPORT_HEIGHT = 1080
    """int: Height of PNG exports
    """

    def __init__(self, name: str, output_dir: Path, formats: Optional[List[str]] = None):
        """
        Initialize the visualization builder.

        Args:
            name: Name of the pangenome for visualization titles and filenames
            output_dir: Directory path where output files will be saved
            formats: List of output formats to generate. Defaults to ["html"]
        """
        self.logger = logging.getLogger("PANORAMA")
        self.name = name
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.formats = formats if formats else [self.DEFAULT_FORMAT]

        # Validate output formats
        invalid_formats = set(self.formats) - set(self.OUTPUT_FORMATS)
        if invalid_formats:
            raise ValueError(f"Unsupported output formats: {invalid_formats}")

        # Core plot components - initialized by subclasses
        self._main_plot: Optional[figure] = None
        self._glyph_renderer: Optional[GlyphRenderer] = None
        self._color_bar: Optional[figure] = None
        self._left_bar: Optional[figure] = None
        self._top_bar: Optional[figure] = None

    # Properties with validation for core components
    @property
    def main_plot(self) -> figure:
        """Get the main heatmap plot figure."""
        if self._main_plot is None:
            raise RuntimeError("Main plot is not initialized yet.")
        return self._main_plot

    @main_plot.setter
    def main_plot(self, plot: figure) -> None:
        """Set the main heatmap plot figure with validation."""
        if not isinstance(plot, figure):
            raise TypeError("Main plot must be a Bokeh figure object.")
        self._main_plot = plot

    @property
    def glyph_renderer(self) -> GlyphRenderer:
        """Get the glyph renderer for the main plot."""
        if self._glyph_renderer is None:
            raise RuntimeError("Glyph renderer is not initialized yet.")
        return self._glyph_renderer

    @glyph_renderer.setter
    def glyph_renderer(self, glyph: GlyphRenderer) -> None:
        """Set the glyph renderer with validation."""
        if not isinstance(glyph, GlyphRenderer):
            raise TypeError("Glyph renderer must be a GlyphRenderer object.")
        self._glyph_renderer = glyph

    @property
    def glyph(self) -> Glyph:
        """Get the glyph from the renderer."""
        return self.glyph_renderer.glyph

    @property
    def color_bar(self) -> figure:
        """Get the color bar figure."""
        if self._color_bar is None:
            raise RuntimeError("Color bar is not initialized yet.")
        return self._color_bar

    @color_bar.setter
    def color_bar(self, color_bar: figure) -> None:
        """Set the color bar figure with validation."""
        if not isinstance(color_bar, figure):
            raise TypeError("Color bar must be a Bokeh figure object.")
        self._color_bar = color_bar

    @property
    def left_bar(self) -> figure:
        """Get the left bar plot figure."""
        if self._left_bar is None:
            raise RuntimeError("Left bar is not initialized yet.")
        return self._left_bar

    @left_bar.setter
    def left_bar(self, left_bar: figure) -> None:
        """Set the left bar figure with validation."""
        if not isinstance(left_bar, figure):
            raise TypeError("Left bar must be a Bokeh figure object.")
        self._left_bar = left_bar

    @property
    def top_bar(self) -> figure:
        """Get the top bar plot figure."""
        if self._top_bar is None:
            raise RuntimeError("Top bar is not initialized yet.")
        return self._top_bar

    @top_bar.setter
    def top_bar(self, top_bar: figure) -> None:
        """Set the top bar figure with validation."""
        if not isinstance(top_bar, figure):
            raise TypeError("Top bar must be a Bokeh figure object.")
        self._top_bar = top_bar

    def _save_figure(self, fig: figure, filename_base: str) -> None:
        """
        Save a Bokeh figure in the specified formats.

        Args:
            fig: The Bokeh figure object to save
            filename_base: Base filename without extension

        Raises:
            Exception: If an unsupported output format is specified
        """
        for fmt in self.formats:
            output_path = self.output_dir / f"{filename_base}.{fmt}"

            if fmt == "html":
                output_file(output_path.absolute().as_posix())
                save(fig)
                self.logger.info(f"Saved {filename_base} visualization in HTML format to {output_path}")

            elif fmt == "png":
                export_png(
                    fig,
                    filename=output_path.absolute().as_posix(),
                    width=self.PNG_EXPORT_WIDTH,
                    height=self.PNG_EXPORT_HEIGHT,
                )
                self.logger.debug(f"Saved {filename_base} visualization in PNG format to {output_path}")

            else:
                raise ValueError(f"Unsupported output format: {fmt}")

    def _configure_plot_style(self) -> None:
        """
        Configure common plot styling for the main figure.

        This method sets up consistent appearance across all visualization types,
        including fonts, colors, and axis properties.
        """
        if self._main_plot is None:
            raise RuntimeError("Cannot configure style: main plot is not initialized")

        # Title styling
        self.main_plot.title.align = "center"
        self.main_plot.title.text_font_size = "20pt"

        # Axis styling
        self.main_plot.axis.axis_label_text_font_size = "16pt"
        self.main_plot.axis.axis_line_color = None
        self.main_plot.axis.major_label_text_font_size = "12px"
        self.main_plot.axis.major_tick_line_color = None
        self.main_plot.axis.major_label_standoff = 1

        # Grid styling
        self.main_plot.grid.grid_line_color = None

        # Configure y-axis (common across all visualizations)
        self.main_plot.yaxis.axis_label = "System name"
        self.main_plot.yaxis.major_label_orientation = 1 / 3

    def _create_main_figure(
        self,
        matrix: pd.DataFrame,
        x_range: Optional[FactorRange] = None,
        y_range: Optional[FactorRange] = None,
        tooltips: Optional[List[Tuple[str, str]]] = None,
    ) -> None:
        """
        Create the main heatmap figure with common configuration.

        Args:
            matrix: Data matrix for determining ranges if not provided
            x_range: X-axis range for the plot. If None, derived from matrix columns
            y_range: Y-axis range for the plot. If None, derived from matrix index
            tooltips: List of tooltip specifications as (label, field) tuples
        """
        # Use matrix data to create ranges if not provided
        if x_range is None:
            x_range = FactorRange(factors=list(matrix.columns))
        if y_range is None:
            y_range = FactorRange(factors=list(matrix.index))

        # Define tools for interaction
        tools = "hover,save,pan,box_zoom,reset,wheel_zoom"

        # Create the main figure
        self.main_plot = figure(
            x_range=x_range,
            y_range=y_range,
            width=self.CENTER_WIDTH,
            height=self.MIDDLE_HEIGHT,
            tools=tools,
            toolbar_location="below",
            tooltips=tooltips,
        )

        # Apply common styling
        self._configure_plot_style()

    def create_left_bar_plot(
        self,
        source: ColumnDataSource,
        matrix: pd.DataFrame,
        y_field: str = "system_name",
        value_field: str = "count",
        color: str = "navy",
    ) -> None:
        """
        Create a horizontal bar plot on the left side of the visualization.

        Args:
            source: ColumnDataSource containing the data for the bars
            matrix: Data matrix for determining the y-range
            y_field: Field name for the y-axis values
            value_field: Field name for the bar values
            color: Color for the bars
        """
        self.left_bar = figure(
            y_range=list(matrix.index),
            width=self.LEFT_WIDTH,
            height=self.MIDDLE_HEIGHT,
            toolbar_location=None,
            tools="",
        )

        # Create horizontal bars
        self.left_bar.hbar(
            y=y_field,
            right=value_field,
            height=0.9,
            source=source,
            color=color,
            alpha=0.6,
        )

        # Add hover tool
        self.left_bar.add_tools(
            HoverTool(
                tooltips=[
                    ("System", f"@{y_field}"),
                    ("Count", f"@{value_field}"),
                ]
            )
        )

        # Configure appearance
        self._configure_bar_plot_style(self.left_bar, x_label="Count", flip_x=True, hide_y_axis=True)

    def create_top_bar_plot(
        self,
        source: ColumnDataSource,
        x_field: str,
        value_field: str = "count",
        color: str = "green",
        x_order: Optional[List[str]] = None,
    ) -> None:
        """
        Create a vertical bar plot on the top of the visualization.

        Args:
            source: ColumnDataSource containing the data for the bars
            x_field: Field name for the x-axis values
            value_field: Field name for the bar values
            color: Color for the bars
            x_order: Custom ordering for x-axis. If None, uses source data order
        """
        # Determine x-axis range
        if x_order is None:
            x_range = list(source.data[x_field])
        else:
            x_range = x_order

        self.top_bar = figure(
            x_range=x_range,
            height=self.TOP_HEIGHT,
            width=self.CENTER_WIDTH,
            toolbar_location=None,
            tools="",
        )

        # Create vertical bars
        self.top_bar.vbar(
            x=x_field,
            top=value_field,
            width=0.9,
            source=source,
            color=color,
            alpha=0.6,
        )

        # Add hover tool
        self.top_bar.add_tools(
            HoverTool(
                tooltips=[
                    (x_field.title(), f"@{x_field}"),
                    ("Count", f"@{value_field}"),
                ]
            )
        )

        # Configure appearance
        self._configure_bar_plot_style(self.top_bar, y_label="Count", hide_x_axis=True)

    @staticmethod
    def _configure_bar_plot_style(
        plot: figure,
        x_label: Optional[str] = None,
        y_label: Optional[str] = None,
        flip_x: bool = False,
        hide_x_axis: bool = False,
        hide_y_axis: bool = False,
    ) -> None:
        """
        Configure styling for bar plots.

        Args:
            plot: The figure to configure
            x_label: Label for x-axis
            y_label: Label for y-axis
            flip_x: Whether to flip the x-axis
            hide_x_axis: Whether to hide the x-axis
            hide_y_axis: Whether to hide the y-axis
        """
        if hide_x_axis:
            plot.xaxis.visible = False
        elif x_label:
            plot.xaxis.axis_label = x_label

        if hide_y_axis:
            plot.yaxis.visible = False
        elif y_label:
            plot.yaxis.axis_label = y_label

        if flip_x:
            plot.x_range.flipped = True

        # Common bar plot styling
        plot.grid.grid_line_color = None
        plot.outline_line_color = None

    @staticmethod
    def _configure_minimal_plot(plot: figure) -> None:
        """
        Configure a minimal plot style (no axes, grid, etc.).

        Args:
            plot: The figure to configure with minimal styling
        """
        plot.xaxis.visible = False
        plot.yaxis.visible = False
        plot.grid.grid_line_color = None
        plot.outline_line_color = None

    @abstractmethod
    def create_main_figure(self, *args, **kwargs) -> None:
        """
        Create the main visualization figure.

        This method must be implemented by subclasses to create their specific
        type of main visualization (correlation matrix, partition matrix, etc.).
        """
        pass

    @abstractmethod
    def plot(self) -> None:
        """
        Create and save the complete visualization layout.

        This method must be implemented by subclasses to define their specific
        layout arrangement and save the final visualization.
        """
        pass
