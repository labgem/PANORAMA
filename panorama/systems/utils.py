#!/usr/bin/env python3

"""
This module provides utility functions to detect and write biological systems in pangenomes.
"""

# default libraries
from __future__ import annotations
import logging
from typing import Dict, Iterable, List, Set, Tuple
from collections import defaultdict
from pathlib import Path

# installed libraries
import numpy as np
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism
from bokeh.plotting import figure
from bokeh.io import export_png, output_file, save
from bokeh.models import (
    BasicTicker,
    ColorBar,
    ColumnDataSource,
    GlyphRenderer,
    Glyph,
    FactorRange,
    HoverTool,
)

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family
from panorama.pangenomes import Pangenome

pd.options.mode.copy_on_write = True  # Remove when pandas3.0 available. See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy


def filter_global_context(
    graph: nx.Graph, jaccard_threshold: float = 0.8
) -> nx.Graph[GeneFamily]:
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
        organisms (Set[Organism]): Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """
    orgs = frozenset(organisms)

    # Precompute family ∩ organisms intersections
    family_orgs: dict[GeneFamily, set[Organism]] = {
        fam: set(fam.organisms).intersection(orgs) for fam in graph.nodes
    }

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
        organisms (Set[Organism]): Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
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
                    and any(
                        name in available_names
                        for name in metadata.secondary_names.split(",")
                    )
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
    if (
        len(mandatory_seen) >= func_unit.min_mandatory
        and len(mandatory_seen | accessory_seen) >= func_unit.min_total
    ):
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


# Note that this function assumes that a family could play multiple roles to satisfy the model requirements if it has multiple annotations.
def check_needed_families(matrix: pd.DataFrame, func_unit: FuncUnit) -> bool:
    """
    Check if there are enough mandatory and total families to satisfy the functional unit rules.

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        Boolean: True if satisfied, False otherwise
    """

    matrix = matrix.loc[matrix.sum(axis=1) > 0]  # Remove all-zero rows

    mandatory_fams = {fam.name for fam in func_unit.mandatory}
    mandatory_count = sum(fam in mandatory_fams for fam in matrix.index)

    return (
        mandatory_count >= func_unit.min_mandatory
        and len(matrix) >= func_unit.min_total
    )


def get_metadata_to_families(
    pangenome: Pangenome, sources: Iterable[str]
) -> Dict[str, Dict[str, Set[GeneFamily]]]:
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
        for exchangeable in exchangeable:
            for source, annotation2families in annot2fam.items():
                if exchangeable in annotation2families:

                    if (
                        fam_model.name in fam2source
                        and fam2source[fam_model.name] != source
                    ):
                        logging.getLogger("PANORAMA").warning(
                            f"Protein annotation {fam_model.name} is encountered in multiple sources."
                            "All sources will be used, but only first one will be associated with "
                            "the model family."
                        )
                    else:
                        fam2source[fam_model.name] = source

                    for gf in annotation2families[exchangeable]:
                        gf2fam[gf].add(fam_model)

    return gf2fam, fam2source


class VisualizationBuilder:
    """
    Builds correlation matrix visualizations with shared configuration and state.

    This class maintains a common configuration and provides a cohesive interface
    for creating all components of the correlation matrix visualization.
    """

    # Constants for plot dimensions
    TOTAL_WIDTH, TOTAL_HEIGHT = 1780, 920

    LEFT_WIDTH = int(0.15 * TOTAL_WIDTH)
    CENTER_WIDTH = int(0.75 * TOTAL_WIDTH)
    RIGHT_WIDTH = int(0.1 * TOTAL_WIDTH)

    TOP_HEIGHT = int(0.15 * TOTAL_HEIGHT)
    MIDDLE_HEIGHT = int(0.7 * TOTAL_HEIGHT)
    BELOW_HEIGHT = int(0.15 * TOTAL_HEIGHT)

    # Output supported formats
    OUTPUT_FORMATS = ["html", "png"]
    DEFAULT_FORMAT = "html"

    # Constants for PNG export dimensions
    PNG_EXPORT_WIDTH = 1920
    PNG_EXPORT_HEIGHT = 1080

    def __init__(self, name: str, output_dir: Path, formats: List[str] = None):
        """
        Initialize the visualization builder.

        Args:
            name: Name of the pangenome for visualization titles.
            output_dir: Directory path where output files will be saved.
            formats: List of output formats to generate.
        """
        self.logger = logging.getLogger("PANORAMA")
        self.name = name
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.formats = formats if formats else ["html"]
        self._main_plot = None
        self._glyph_renderer = None
        self._color_bar = None
        self._left_bar = None

    def _save_figure(self, fig, filename_base: str) -> None:
        """
        Save a Bokeh figure in the specified formats.

        Args:
            fig: The Bokeh figure object to save.
            filename_base: Base filename without extension.
        """
        for fmt in self.formats:
            if fmt == "html":
                output_path = self.output_dir / f"{filename_base}.html"
                output_file(output_path.absolute().as_posix())
                save(fig)
                self.logger.info(
                    f"Saved {filename_base} heatmap in HTML format to {output_path}"
                )

            elif fmt == "png":
                output_path = self.output_dir / f"{filename_base}.png"
                output_file(output_path.absolute().as_posix())
                export_png(
                    fig,
                    width=VisualizationBuilder.PNG_EXPORT_WIDTH,
                    height=VisualizationBuilder.PNG_EXPORT_HEIGHT,
                )
                self.logger.debug(
                    f"Saved {filename_base} heatmap in PNG format to {output_path}"
                )

            else:
                raise Exception(f"Unsupported output format: {fmt}")

    @property
    def main_plot(self) -> figure:
        if self._main_plot is None:
            raise Exception("plot is not initialized yet.")

        return self._main_plot

    @main_plot.setter
    def main_plot(self, plot: figure):
        if not isinstance(plot, figure):
            raise Exception("plot must be a Bokeh figure object.")

        self._main_plot = plot

    @property
    def glyph_renderer(self) -> GlyphRenderer:
        if self._glyph_renderer is None:
            raise Exception("Glyph renderer is not initialized yet.")

        return self._glyph_renderer

    @glyph_renderer.setter
    def glyph_renderer(self, glyph: GlyphRenderer):
        if not isinstance(glyph, GlyphRenderer):
            raise Exception("Glyph renderer must be a GlyphRenderer object.")

        self._glyph_renderer = glyph

    def _configure_plot_style(self) -> None:
        """Configure common plot styling."""
        self.main_plot.title.align = "center"
        self.main_plot.title.text_font_size = "20pt"
        self.main_plot.axis.axis_label_text_font_size = "16pt"
        self.main_plot.axis.axis_line_color = None
        self.main_plot.axis.major_label_text_font_size = "12px"
        self.main_plot.grid.grid_line_color = None
        self.main_plot.axis.major_tick_line_color = None
        self.main_plot.axis.major_label_standoff = 1

        # Set axis labels
        self.main_plot.yaxis.axis_label = "System name"
        self.main_plot.yaxis.major_label_orientation = 1 / 3
        # xaxis is configure in child class

    def _create_main_figure(
        self,
        matrix: pd.DataFrame,
        x_range: FactorRange = None,
        y_range: FactorRange = None,
        tooltips: List[Tuple[str, str]] = None,
    ) -> None:
        """
        Create the main correlation matrix heatmap figure.

        Args:
            matrix: Preprocessed correlation matrix.
            x_range: X-axis range for the plot.
            y_range: Y-axis range for the plot.
        """
        # Create the main figure
        tools = "hover,save,pan,box_zoom,reset,wheel_zoom"
        self.main_plot = figure(
            x_range=x_range,
            y_range=y_range,
            width=VisualizationBuilder.CENTER_WIDTH,
            height=VisualizationBuilder.MIDDLE_HEIGHT,
            tools=tools,
            toolbar_location="below",
            tooltips=tooltips,
        )

        # Configure plot appearance
        self._configure_plot_style()

    @property
    def glyph(self) -> Glyph:
        return self.glyph_renderer.glyph

    @property
    def color_bar(self) -> figure:
        if self._color_bar is None:
            raise Exception("Color bar is not initialized yet.")
        return self._color_bar

    @color_bar.setter
    def color_bar(self, color_bar: figure):
        if not isinstance(color_bar, figure):
            raise Exception("Color bar must be a Bokeh figure object.")
        self._color_bar = color_bar

    @property
    def left_bar(self) -> figure:
        if self._left_bar is None:
            raise Exception("Left bar is not initialized yet.")
        return self._left_bar

    @left_bar.setter
    def left_bar(self, left_bar: figure):
        if not isinstance(left_bar, figure):
            raise Exception("Left bar must be a Bokeh figure object.")
        self._left_bar = left_bar

    def create_left_bar_plot(self, source: ColumnDataSource, matrix: pd.DataFrame):
        self.left_bar = figure(
            y_range=list(matrix.index),
            width=VisualizationBuilder.LEFT_WIDTH,
            height=VisualizationBuilder.MIDDLE_HEIGHT,
            toolbar_location=None,
            tools="",
        )
        self.left_bar.hbar(
            y="system_name",
            right="count",
            height=0.9,
            source=source,
            color="navy",
            alpha=0.6,
        )
        self.left_bar.add_tools(HoverTool(tooltips=[
            ("system", "@system_name"),
            ("count", "@count"),
        ],))

        # Configure the left bar plot
        self.left_bar.yaxis.visible = False
        self.left_bar.xaxis.axis_label = "Count"
        self.left_bar.x_range.flipped = True
        self.left_bar.grid.grid_line_color = None
        self.left_bar.outline_line_color = None
