#!/usr/bin/env python3
# coding:utf-8

"""
Systems Comparison Module for PANORAMA

This module provides functionality to compare genomic systems across multiple
pangenomes, compute similarity metrics, and generate visualizations.
"""

# Standard library imports
from __future__ import annotations
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from shutil import rmtree
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock
from collections import defaultdict

# Third-party imports
from tqdm import tqdm
import networkx as nx
import pandas as pd
from bokeh.io import output_file, save, export_png
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, ColorBar
from bokeh.transform import linear_cmap
from bokeh.palettes import Reds256

# Local imports
from panorama.pangenomes import Pangenomes
from panorama.utils import mkdir, init_lock
from panorama.utility.utility import check_models
from panorama.systems.system import System, ClusterSystems
from panorama.systems.write_systems import check_pangenome_write_systems
from panorama.compare.utils import (
    parser_comparison,
    common_launch,
    cluster_on_gfrr,
    compute_gfrr,
)


logger = logging.getLogger("PANORAMA")


class SystemsComparisonError(Exception):
    """Custom exception for systems comparison errors."""

    pass


def check_compare_systems_args(args: argparse.Namespace) -> Dict:
    """
    Validate and prepare arguments for systems comparison.

    Args:
        args: Command line arguments containing sources, models, and other parameters.

    Returns:
        Dict containing required information flags and parameters for pangenome processing.

    Raises:
        argparse.ArgumentError: If the number of sources and models don't match.
    todo:
        - check if output is a directory and if empty
        - incompatibility with gfrr_metric not given & graph_formats given
    """
    # Define what information we need from pangenomes
    need_info = {
        "need_annotations": True,
        "need_families": True,
        "need_families_info": True,
        "need_metadata": True,
        "metatypes": ["families"],
        "need_systems": True,
        "systems_sources": args.sources,
        "read_canonical": args.canonical,
    }

    # Validate that sources and models lists have the same length
    if len(args.sources) != len(args.models):
        raise argparse.ArgumentError(
            argument=None,
            message=f"Number of sources ({len(args.sources)}) and models ({len(args.models)}) must match.",
        )

    return need_info


def add_system_metadata_to_graph(pangenomes: Pangenomes, graph: nx.Graph) -> None:
    """
    Add system metadata as node attributes to the systems graph.

    Args:
        pangenomes: Collection of pangenomes containing systems.
        graph: NetworkX graph to add metadata to.
    """
    for pangenome_name, pangenome in pangenomes.items():
        for system in pangenome.systems:
            # Create a comprehensive system metadata
            system_metadata = {
                "pangenome": pangenome_name,
                "system_name": system.name,
                "system_id": system.ID,
                "families_models_count": system.number_of_model_gene_families,
                "families_count": system.number_of_families,
            }

            # Generate unique hash for the system
            system_hash = hash((pangenome_name, system.name, system.ID))

            # Update node attributes if node exists
            if graph.has_node(system_hash):
                graph.nodes[system_hash].update(system_metadata)


def create_pangenome_system_graph(
    pangenome,
) -> Tuple[nx.Graph, Dict[int, str], Dict[int, System]]:
    """
    Create a graph representation of systems for a single pangenome.

    Args:
        pangenome: Pangenome object containing systems.

    Returns:
        Tuple containing:
            - NetworkX graph with system nodes
            - Dictionary mapping system hash to pangenome name
            - Dictionary mapping system hash to a system object
    """
    graph = nx.Graph()
    system_to_pangenome = {}
    system_hash_to_system = {}

    # Add each system as a node in the graph
    for system in pangenome.systems:
        system_hash = hash((pangenome.name, system.name, system.ID))

        # Add a node with basic attributes
        graph.add_node(
            system_hash,
            system_id=system.ID,
            system_name=system.name,
            pangenome=pangenome.name,
        )

        # Maintain mapping dictionaries
        system_to_pangenome[system_hash] = pangenome.name
        system_hash_to_system[system_hash] = system

    return graph, system_to_pangenome, system_hash_to_system


def create_systems_graph(
    pangenomes: Pangenomes,
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> Tuple[nx.Graph, Dict[int, str], Dict[int, System]]:
    """
    Create a comprehensive graph of all systems across pangenomes.

    Args:
        pangenomes: Collection of pangenomes to process.
        threads: Number of threads for parallel processing.
        lock: Thread lock for synchronization.
        disable_bar: Whether to disable the progress bar.

    Returns:
        Tuple containing:
            - NetworkX graph with all system nodes
            - Dictionary mapping system hash to pangenome name
            - Dictionary mapping system hash to a system object
    """

    with ThreadPoolExecutor(
        max_workers=threads, initializer=init_lock, initargs=(lock,)
    ) as executor:
        with tqdm(total=len(pangenomes), unit="Pangenome", disable=disable_bar) as pbar:
            # Submit tasks for each pangenome
            futures = []
            for pangenome in pangenomes:
                logger.debug(f"Processing systems for pangenome {pangenome.name}")
                future = executor.submit(create_pangenome_system_graph, pangenome)
                future.add_done_callback(lambda p: pbar.update())
                futures.append(future)

            # Combine results from all pangenomes
            combined_graph = nx.Graph()
            combined_system_to_pangenome = {}
            combined_system_hash_to_system = {}

            for future in futures:
                graph, sys_to_pan, sys_hash_to_sys = future.result()
                combined_graph.add_nodes_from(graph)
                combined_system_to_pangenome.update(sys_to_pan)
                combined_system_hash_to_system.update(sys_hash_to_sys)

    return combined_graph, combined_system_to_pangenome, combined_system_hash_to_system


def compute_gfrr_edges(
    graph: nx.Graph,
    system_to_pangenome: Dict[int, str],
    system_hash_to_system: Dict[int, System],
    gfrr_cutoff: Tuple[float, float] = (0.8, 0.8),
    gfrr_models_cutoff: Tuple[float, float] = (0.8, 0.8),
    disable_bar: bool = False,
) -> None:
    """
    Compute GFRR (Gene Families Repertoire Relatedness) edges between systems from different pangenomes.

    GFRR is a similarity metric that compares gene family repertoires between systems.
    Edges are only added between systems that meet both GFRR model and GFRR cutoff thresholds.

    Args:
        graph: Graph with system nodes (edges will be added).
        system_to_pangenome: Mapping from system hash to pangenome name.
        system_hash_to_system: Mapping from system hash to a system object.
        gfrr_cutoff: Minimum (min_gfrr, max_gfrr) for all gene families.
        gfrr_models_cutoff: Minimum (min_gfrr, max_gfrr) for model gene families.
        disable_bar: Whether to disable the progress bar.
    """
    # Generate all possible pairs of systems from different pangenomes
    system_pairs = [
        (sys1_hash, sys2_hash)
        for sys1_hash, sys2_hash in combinations(graph.nodes, 2)
        if system_to_pangenome[sys1_hash] != system_to_pangenome[sys2_hash]
    ]

    logger.info(f"Computing GFRR for {len(system_pairs)} system pairs")

    with tqdm(
        total=len(system_pairs),
        unit="system pair",
        desc="Computing FRR",
        disable=disable_bar,
    ) as pbar:
        for sys1_hash, sys2_hash in system_pairs:
            # Skip if edge already exists
            if graph.has_edge(sys1_hash, sys2_hash):
                pbar.update()
                continue

            sys1 = system_hash_to_system[sys1_hash]
            sys2 = system_hash_to_system[sys2_hash]

            # First, check model families GFRR (more restrictive)
            min_gfrr_models, max_gfrr_models, shared_models_families = compute_gfrr(
                set(sys1.model_families), set(sys2.model_families)
            )

            # Only proceed if model families meet cutoff
            if (
                min_gfrr_models > gfrr_models_cutoff[0]
                and max_gfrr_models > gfrr_models_cutoff[1]
            ):
                # Check all families FRR
                min_gfrr, max_gfrr, shared_families = compute_gfrr(
                    set(sys1.families), set(sys2.families)
                )

                # Add edge if both cutoffs are met
                if min_gfrr > gfrr_cutoff[0] and max_gfrr > gfrr_cutoff[1]:
                    graph.add_edge(
                        sys1_hash,
                        sys2_hash,
                        min_gfrr_models=min_gfrr_models,
                        max_gfrr_models=max_gfrr_models,
                        shared_model_families=shared_models_families,
                        min_gfrr=min_gfrr,
                        max_gfrr=max_gfrr,
                        shared_families=shared_families,
                    )

            pbar.update()


def write_conserved_systems(
    pangenomes: Pangenomes,
    output: Path,
    conserved_systems_graph: nx.Graph,
    graph_formats: List[str],
) -> None:
    """
    Write a conserved systems graph to various output formats.

    Args:
        pangenomes: Collection of pangenomes.
        output: Output directory path.
        conserved_systems_graph: Graph of conserved systems.
        graph_formats: List of output formats ('gexf', 'graphml').
    """

    # Add metadata to graph nodes
    add_system_metadata_to_graph(pangenomes, conserved_systems_graph)

    # Write in GEXF format if requested
    if "gexf" in graph_formats:
        gexf_file = output / "conserved_systems.gexf"
        logger.info(f"Writing conserved systems graph in GEXF format: {gexf_file}")
        nx.readwrite.gexf.write_gexf(conserved_systems_graph, gexf_file)

    # Write in GraphML format if requested
    if "graphml" in graph_formats:
        graphml_file = output / "conserved_systems.graphml"
        logger.info(
            f"Writing conserved systems graph in GraphML format: {graphml_file}"
        )
        nx.readwrite.graphml.write_graphml(conserved_systems_graph, graphml_file)

    # Write a summary TSV file
    tsv_file = output / "conserved_systems.tsv"
    logger.info(f"Writing conserved systems summary: {tsv_file}")
    # TODO: Implement TSV writing logic


def get_pangenomes_to_systems_data(pangenomes: Pangenomes) -> Dict[str, Dict[str, int]]:
    """
    Extract system occurrence data from pangenomes.

    Args:
        pangenomes: Collection of pangenomes.

    Returns:
        Dictionary mapping pangenome names to system occurrence counts.
    """
    data = {}

    for pangenome in pangenomes:
        data[pangenome.name] = defaultdict(int)

        # Count occurrences of each system type
        for system in pangenome.systems:
            data[pangenome.name][system.name] += 1

    return data


def compare_systems(
    pangenomes: Pangenomes,
    gfrr_metrics: str = "min_gfrr_models",
    gfrr_cutoff: Tuple[float, float] = (0.8, 0.8),
    gfrr_models_cutoff: Tuple[float, float] = (0.8, 0.8),
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> nx.Graph:
    """
    Compare systems across pangenomes and identify conserved system clusters.

    Args:
        pangenomes: Collection of pangenomes to compare.
        gfrr_metrics: Metric to use for clustering ('min_gfrr_models', 'max_gfrr_models', etc.).
        gfrr_cutoff: GFRR cutoff thresholds for all gene families.
        gfrr_models_cutoff: GFRR cutoff thresholds for model gene families.
        threads: Number of threads for parallel processing.
        lock: Thread lock for synchronization.
        disable_bar: Whether to disable progress bars.

    Returns:
        NetworkX graph containing conserved systems clusters.

    Raises:
        SystemsComparisonError: If comparison fails.
    """
    logger.info("Starting systems comparison across pangenomes")

    try:
        # Create a comprehensive systems graph
        systems_graph, system_to_pangenome, system_hash_to_system = (
            create_systems_graph(pangenomes, threads, lock, disable_bar)
        )

        logger.info(f"Created systems graph with {len(systems_graph.nodes)} nodes")

        # Compute similarity edges based on FRR
        compute_gfrr_edges(
            systems_graph,
            system_to_pangenome,
            system_hash_to_system,
            gfrr_cutoff,
            gfrr_models_cutoff,
            disable_bar,
        )

        logger.info(f"Added {len(systems_graph.edges)} similarity edges")

        # Cluster systems based on similarity
        partitions = cluster_on_gfrr(systems_graph, gfrr_metrics)
        logger.info(f"Found {len(partitions)} system clusters")

        # Process clusters and add to pangenomes
        conserved_clusters_count = 0
        for cluster_id, cluster_systems in enumerate(partitions, start=1):
            if len(cluster_systems) > 1:  # Only keep multi-system clusters
                conserved_clusters_count += 1
                cluster_system_objects = set()

                # Update node attributes for systems in this cluster
                for system_hash in cluster_systems:
                    system_obj = system_hash_to_system[system_hash]
                    pangenome_name = system_to_pangenome[system_hash]

                    # Clean up temporary clustering attributes
                    node_attrs = systems_graph.nodes[system_hash]
                    node_attrs.pop(f"{gfrr_metrics}_cluster", None)

                    # Add cluster information
                    node_attrs.update(
                        {
                            "cluster_systems_id": cluster_id,
                            "system_id": system_obj.ID,
                            "system_name": system_obj.name,
                            "pangenome": pangenome_name,
                        }
                    )

                    cluster_system_objects.add(system_obj)

                # Add a cluster to a pangenomes collection
                pangenomes.add_cluster_systems(
                    ClusterSystems(cluster_id, *cluster_system_objects)
                )
            else:
                # Remove single-system "clusters" from the graph
                systems_graph.remove_nodes_from(cluster_systems)

        logger.info(f"Identified {conserved_clusters_count} conserved system clusters")
        return systems_graph

    except Exception as e:
        logger.error(f"Systems comparison failed: {str(e)}")
        raise SystemsComparisonError(f"Failed to compare systems: {str(e)}") from e


def generate_heatmap(
    data: pd.DataFrame,
    output: Path,
    output_name: str,
    output_formats: List[str],
    figure_size: Tuple[float, float],
    title: str = "Heatmap",
    font_size: int = 18,
) -> None:
    """
    Generate a heatmap visualization using Bokeh.

    Args:
        data: Input data for the heatmap (should be a DataFrame).
        output: Directory to save the heatmap files.
        output_name: Base name for output files.
        output_formats: List of formats to save ('html', 'png').
        figure_size: Size of the figure in pixels (width, height).
        title: Title for the heatmap.
        font_size: Base font size for text elements.
    """

    # Prepare a data source for Bokeh
    melted_data = data.reset_index().melt(
        id_vars="index", var_name="columns", value_name="value"
    )
    source = ColumnDataSource(melted_data)

    # Create a figure with appropriate dimensions
    plot_width, plot_height = int(figure_size[0]), int(figure_size[1])

    heatmap = figure(
        title=title,
        x_range=list(data.columns),
        y_range=list(data.index)[::-1],  # Reverse y-axis for proper orientation
        width=plot_width,
        height=plot_height,
        tools="hover,save,pan,box_zoom,reset,wheel_zoom",
        toolbar_location="above",
    )

    # Create color mapping
    color_mapper = linear_cmap(
        field_name="value",
        palette=Reds256[::-1],  # Reverse palette for intuitive coloring
        low=data.values.min(),
        high=data.values.max(),
    )

    # Add heatmap rectangles
    heatmap.rect(
        x="columns",
        y="index",
        width=1,
        height=1,
        source=source,
        fill_color=color_mapper,
        line_color=None,
    )

    # Add color bar legend
    color_bar = ColorBar(
        color_mapper=color_mapper["transform"], width=8, location=(0, 0)
    )
    heatmap.add_layout(color_bar, "right")

    # Configure styling
    heatmap.title.text_font_size = f"{font_size + 6}px"
    heatmap.xaxis.axis_label = "Systems"
    heatmap.yaxis.axis_label = "Species"
    heatmap.xaxis.axis_label_text_font_size = f"{font_size + 2}px"
    heatmap.yaxis.axis_label_text_font_size = f"{font_size + 2}px"
    heatmap.xaxis.major_label_text_font_size = f"{font_size}px"
    heatmap.xaxis.major_label_orientation = 1  # Rotate x-axis labels
    heatmap.yaxis.major_label_text_font_size = f"{font_size}px"

    # Save in requested formats
    if "html" in output_formats:
        html_path = output / f"{output_name}.html"
        output_file(html_path)
        save(heatmap)
        logger.debug(f"Saved heatmap in HTML format: {html_path}")

    if "png" in output_formats:
        png_path = output / f"{output_name}.png"
        export_png(heatmap, filename=str(png_path), width=1920, height=1080)
        logger.debug(f"Saved heatmap in PNG format: {png_path}")


def create_pangenome_systems_heatmaps(pangenomes: Pangenomes, output: Path) -> None:
    """
    Generate heatmaps showing system distribution across pangenomes.

    Args:
        pangenomes: Collection of pangenomes to analyze.
        output: Directory to save the generated heatmaps.
    """
    logger.info("Generating pangenome systems heatmaps")

    # Prepare data matrix
    systems_data = get_pangenomes_to_systems_data(pangenomes)
    data_df = pd.DataFrame(systems_data).fillna(0)

    # Standard figure size for web display
    figure_size = (1200, 900)

    # Generate raw counts heatmap
    generate_heatmap(
        data_df.T,
        output,
        "heatmap_number_systems",
        ["html"],
        figure_size=figure_size,
        title="Number of Systems Detected in Pangenomes",
        font_size=18,
    )

    # Generate a normalized percentage heatmap
    data_normalized = data_df.div(data_df.sum(axis=0), axis=1) * 100
    generate_heatmap(
        data_normalized.T,
        output,
        "heatmap_normalized_systems",
        ["html"],
        figure_size=figure_size,
        title="Normalized Percentage of Systems Detected in Pangenomes",
        font_size=18,
    )

    logger.info("Heatmaps generation completed")


def launch(args: argparse.Namespace) -> None:
    """
    Main entry point for systems comparison analysis.

    Args:
        args: Command line arguments containing all parameters for the analysis.

    Raises:
        SystemsComparisonError: If analysis fails at any stage.
    """
    logger.info("Starting PANORAMA systems comparison analysis")

    try:
        # Validate arguments and prepare requirements
        need_info = check_compare_systems_args(args)

        # Process and validate model files
        models_list = []
        for model_path in args.models:
            validated_models = check_models(
                model_path, disable_bar=args.disable_prog_bar
            )
            models_list.append(validated_models)
        need_info["models"] = models_list

        # Load pangenomes and set up a working environment
        pangenomes, tmp_dir, _, lock = common_launch(
            args, check_pangenome_write_systems, need_info, sources=args.sources
        )

        # Create output directory
        output_dir = mkdir(args.output, force=args.force)
        logger.info(f"Output directory: {output_dir}")

        # Generate heatmaps if requested
        if args.heatmap:
            logger.info("Generating system distribution heatmaps")
            create_pangenome_systems_heatmaps(pangenomes, output_dir)

        # Perform systems comparison if GFRR metrics specified
        if args.gfrr_metrics:
            logger.info(
                f"Performing systems comparison using {args.gfrr_metrics} metric"
            )
            conserved_systems_graph = compare_systems(
                pangenomes,
                gfrr_metrics=args.gfrr_metrics,
                gfrr_cutoff=args.gfrr_cutoff,
                gfrr_models_cutoff=args.gfrr_models_cutoff,
                threads=args.cpus,
                lock=lock,
                disable_bar=args.disable_prog_bar,
            )

            # Write results to files
            write_conserved_systems(
                pangenomes, output_dir, conserved_systems_graph, args.graph_formats
            )

        logger.info("Systems comparison analysis completed successfully")

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise SystemsComparisonError(
            f"Systems comparison analysis failed: {str(e)}"
        ) from e

    finally:
        # Clean up temporary files
        if not args.keep_tmp:
            logger.debug(f"Cleaning up temporary directory: {tmp_dir}")
            rmtree(tmp_dir, ignore_errors=True)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Create a subparser for systems comparison command.

    Args:
        sub_parser: Subparser action from main argument parser.

    Returns:
        Configured argument parser for systems comparison.
    """
    parser = sub_parser.add_parser(
        "compare_systems",
        description="Compare genomic systems among pangenomes using GFRR metrics",
        help="Identify conserved systems across multiple pangenomes",
    )

    parser_comparison_systems(parser)
    return parser


def parser_comparison_systems(parser: argparse.ArgumentParser) -> None:
    """
    Configure argument parser for systems comparison command.

    Args:
        parser: Argument parser to configure.
    """
    # Get common comparison arguments
    required, compare_opt, optional = parser_comparison(parser)

    # Required arguments specific to systems comparison
    required.add_argument(
        "-m",
        "--models",
        required=True,
        type=Path,
        nargs="+",
        help="Path(s) to model list files. Multiple models can be specified "
        "corresponding to different sources. Order must match --sources.",
    )

    required.add_argument(
        "-s",
        "--sources",
        required=True,
        type=str,
        nargs="+",
        help="Name(s) of the systems sources. Multiple sources can be specified. "
        "Order must match --models argument.",
    )

    # Optional comparison-specific arguments
    compare_opt.add_argument(
        "--heatmap",
        action="store_true",
        help="Generate heatmaps showing normalized system presence distribution across pangenomes",
    )

    compare_opt.add_argument(
        "--gfrr_metrics",
        type=str,
        default=None,
        choices=["min_gfrr_models", "max_gfrr_models", "min_gfrr", "max_gfrr"],
        help="Similarity metric for clustering conserved systems. "
        "Models metrics use only model gene families, while regular metrics use all families.",
    )

    compare_opt.add_argument(
        "--gfrr_models_cutoff",
        type=float,
        nargs=2,
        default=[0.4, 0.6],
        help="GFRR cutoff thresholds for model gene families comparison. "
        "min_gfrr = shared_families / min(families1, families2), "
        "max_gfrr = shared_families / max(families1, families2). "
        "Default: 0.2 0.2",
    )

    # Additional optional arguments
    optional.add_argument(
        "--canonical",
        action="store_true",
        help="Include canonical system versions in the analysis output.",
    )
