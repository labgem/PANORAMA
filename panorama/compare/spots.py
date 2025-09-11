#!/usr/bin/env python3
# coding:utf-8

"""
Pangenome spots comparison and conserved spots detection module.

This module provides comprehensive functionality for comparing spots across multiple pangenomes,
detecting conserved genomic regions, and analyzing systems relationships through graph-based
clustering approaches. It includes utilities for building comparative graphs, computing gene
family relatedness relationships (GFRR), and generating various output formats for visualization.
"""

# TODO create a Heatmap with clustered spots and pangenomes (see code of Jupyter NB)

# Default libraries
from __future__ import annotations
import argparse
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Optional
from shutil import rmtree
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock

# Installed libraries
from tqdm import tqdm
import networkx as nx
import numpy as np
import pandas as pd

# Local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.region import Spot, ConservedSpots
from panorama.utils import mkdir, init_lock
from panorama.compare.utils import (
    parser_comparison,
    common_launch,
    cluster_on_gfrr,
    compute_gfrr,
)
from panorama.utility.utility import check_models

# Configure logger
logger = logging.getLogger("PANORAMA")

# TODO: ensure to write the canonical as expected


def check_compare_spots_args(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Validate and configure arguments for spots comparison analysis.

    This function analyzes command-line arguments to determine required data components
    for pangenome spots comparison. It validates the consistency of system-related
    arguments and configures the information dictionary for downstream processing.

    Args:
        args (argparse.Namespace): Parsed command-line arguments containing:
            - systems (bool): Whether to include systems analysis
            - models (Optional[List[Path]]): Paths to model files
            - sources (Optional[List[str]]): System source names
            - canonical (bool): Whether to include canonical systems
            - disable_prog_bar (bool): Progress bar control

    Returns:
        Dict[str, Any]: Configuration dictionary specifying required data:
            - need_annotations (bool): Whether annotations are required
            - need_families (bool): Whether gene families are required
            - need_families_info (bool): Whether family metadata is required
            - need_rgp (bool): Whether regions of genomic plasticity are required
            - need_spots (bool): Whether spots are required
            - Additional system-specific requirements if systems analysis is enabled

    Raises:
        argparse.ArgumentError: If there's a mismatch in system-related arguments:
            - Systems requested but no models provided
            - Systems requested but no sources provided
            - Models provided but no sources (or vice versa)
    """
    # Base requirements for spots comparison
    need_info = {
        "need_annotations": True,  # Required for gene family analysis
        "need_families": True,  # Required for spot border analysis
        "need_families_info": True,  # Required for family metadata
        "need_rgp": True,  # Required for spot construction
        "need_spots": True,  # Required for spot comparison
    }

    # Validate and configure systems analysis requirements
    if args.systems:
        logger.info("Systems analysis requested - validating arguments")

        if args.models is not None:
            if args.sources is None:
                raise argparse.ArgumentError(
                    argument=None,
                    message="Systems analysis requires both models and sources. "
                    "Models provided but sources missing. "
                    "Use --sources to specify system sources.",
                )

            # Validate and prepare model files
            validated_models = []
            for model_path in args.models:
                try:
                    validated_model = check_models(
                        model_path, disable_bar=args.disable_prog_bar
                    )
                    validated_models.append(validated_model)
                    logger.debug(f"Validated model file: {model_path}")
                except Exception as e:
                    raise argparse.ArgumentError(
                        argument=None, message=f"Invalid model file {model_path}: {e}"
                    )

            # Configure systems-specific requirements
            need_info.update(
                {
                    "need_metadata": True,  # Required for system metadata
                    "metatypes": ["families"],  # Type of metadata needed
                    "need_systems": True,  # Enable systems loading
                    "systems_sources": args.sources,  # System source identifiers
                    "models": validated_models,  # Validated model objects
                    "read_canonical": args.canonical,  # Include canonical systems
                }
            )

            logger.info(f"Configured systems analysis with {len(args.sources)} sources")

        else:
            # Handle cases where systems requested but models/sources missing
            if args.sources is not None:
                raise argparse.ArgumentError(
                    argument=None,
                    message="Systems analysis requires both models and sources. "
                    "Sources provided but models missing. "
                    "Use --models to specify model files.",
                )
            else:
                raise argparse.ArgumentError(
                    argument=None,
                    message="Systems analysis requires both models and sources. "
                    "Use --models and --sources to specify system parameters.",
                )

    return need_info


def check_pangenome_cs(
    pangenome: Pangenome, sources: Optional[List[str]] = None
) -> None:
    """
    Validate pangenome status for conserved spots analysis.

    This function ensures that a pangenome has the required components computed
    for conserved spots analysis, including RGPs (Regions of Genomic Plasticity)
    and spots. Optionally validates systems detection if sources are provided.

    Args:
        pangenome (Pangenome): Pangenome object to validate.
        sources (Optional[List[str]]): List of system sources to validate.
                                     If None, systems validation is skipped.

    Raises:
        ValueError: If RGPs or spots haven't been computed for the pangenome.
                   Includes guidance on required commands.
        AttributeError: If systems detection hasn't been performed when sources
                       are specified.
        KeyError: If specified sources aren't found in the pangenome's systems.
    """
    pangenome_name = pangenome.name

    # Validate RGP (Regions of Genomic Plasticity) status
    rgp_status = pangenome.status.get("predictedRGP", "NotComputed")
    if rgp_status not in ["inFile", "Computed", "Loaded"]:
        raise ValueError(
            f"RGPs have not been predicted for pangenome '{pangenome_name}'. "
            f"Current status: {rgp_status}. "
            "Please run 'ppanggolin rgp' command to predict RGPs before "
            "proceeding with conserved spots analysis."
        )

    # Validate spots status
    spots_status = pangenome.status.get("spots", "NotComputed")
    if spots_status not in ["inFile", "Computed", "Loaded"]:
        raise ValueError(
            f"Spots have not been predicted for pangenome '{pangenome_name}'. "
            f"Current status: {spots_status}. "
            "Please run 'ppanggolin spot' command to predict spots before "
            "proceeding with conserved spots analysis."
        )

    # Validate systems if sources are specified
    if sources is not None:
        systems_status = pangenome.status.get("systems", "NotDetected")
        if systems_status != "inFile":
            raise AttributeError(
                f"Systems have not been detected for pangenome '{pangenome_name}'. "
                f"Current status: {systems_status}. "
                "Use 'panorama detect' subcommand to detect systems before "
                "including systems in conserved spots analysis."
            )

        # Check if all requested sources are available
        available_sources = pangenome.status.get("systems_sources", [])
        missing_sources = set(sources) - set(available_sources)

        if missing_sources:
            logger.error(
                f"Available systems sources in pangenome '{pangenome_name}': "
                f"{available_sources}"
            )
            raise KeyError(
                f"Missing system sources in pangenome '{pangenome_name}': "
                f"{list(missing_sources)}. "
                f"Use 'panorama detect' subcommand to detect systems for "
                f"these sources in '{pangenome_name}'."
            )

    logger.debug(f"Pangenome '{pangenome_name}' validation successful")


def create_pangenome_spots_graph(
    pangenome: Pangenome, dup_margin: float = 0.05
) -> Tuple[nx.Graph, Dict[int, Set[GeneFamily]], Dict[int, str], Dict[int, Spot]]:
    """
    Create a graph representation of spots from a single pangenome.

    This function constructs a NetworkX graph where nodes represent spots from
    the given pangenome. Each spot is assigned a unique hash identifier and
    associated with its bordering gene families for subsequent comparison analysis.

    Args:
        pangenome (Pangenome): Source pangenome containing spots to process.
        dup_margin (float): Minimum ratio of organisms that must contain multiple
                          copies of a gene family for it to be considered duplicated.
                          Used for multigenic family detection. Default: 0.05 (5%).

    Returns:
        Tuple[nx.Graph, Dict[int, Set[GeneFamily]], Dict[int, str], Dict[int, Spot]]:
            - graph (nx.Graph): NetworkX graph with spot nodes (no edges)
            - spots2borders (Dict[int, Set[GeneFamily]]): Mapping of spot hashes to
              their bordering gene families
            - spots2pangenome (Dict[int, str]): Mapping of spot hashes to pangenome names
            - spothash2spot (Dict[int, Spot]): Mapping of spot hashes to spot objects
    """

    def get_borders_families() -> Set[GeneFamily]:
        """
        Extract all bordering gene families from the current spot.

        Bordering families are those located at the edges of genomic spots,
        which are crucial for determining spot similarity and conservation.

        Returns:
            Set[GeneFamily]: Set of unique gene families found at spot borders.
        """
        borders_families = set()

        # Get spot borders using pangenome-specific parameters
        spot_set_size = pangenome.parameters.get("spot", {}).get("set_size", 3)
        borders = spot.borders(spot_set_size, multigenic)

        # Extract families from both left and right borders
        for _, border_pair in borders:
            left_border, right_border = border_pair
            borders_families.update(left_border)  # Left border families
            borders_families.update(right_border)  # Right border families

        return borders_families

    # Initialize data structures
    graph = nx.Graph()
    spots2borders = {}
    spots2pangenome = {}
    spothash2spot = {}

    # Identify multigenic families (families present in multiple copies)
    # This information is used for border detection
    multigenic = pangenome.get_multigenics(dup_margin=dup_margin)

    logger.debug(
        f"Processing {pangenome.number_of_spots} spots from pangenome '{pangenome.name}'"
    )

    # Process each spot in the pangenome
    for spot in pangenome.spots:
        # Create a unique hash identifier for the spot
        spot_hash = hash((spot.ID, pangenome.name))

        # Add a spot as a node to the graph with metadata
        graph.add_node(
            spot_hash,
            spot_id=spot.ID,
            pangenome=pangenome.name,
            num_regions=len(spot),
            num_families=spot.number_of_families,
        )

        # Store spot mappings for later use
        spots2borders[spot_hash] = get_borders_families()
        spots2pangenome[spot_hash] = pangenome.name
        spothash2spot[spot_hash] = spot

    logger.debug(f"Created graph with {len(graph.nodes)} spot nodes")

    return graph, spots2borders, spots2pangenome, spothash2spot


def create_spots_graph(
    pangenomes: Pangenomes,
    dup_margin: float = 0.05,
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> Tuple[nx.Graph, Dict[int, Set[GeneFamily]], Dict[int, str], Dict[int, Spot]]:
    """
    Create a comprehensive spots graph from multiple pangenomes.

    This function processes all pangenomes in parallel to construct a unified graph
    where nodes represent spots from all pangenomes. The graph serves as the foundation
    for conserved spots detection through subsequent edge computation and clustering.

    Args:
        pangenomes (Pangenomes): Collection of pangenomes to process.
        dup_margin (float): Minimum ratio for multigenic family detection. Default: 0.05.
        threads (int): Number of threads for parallel processing. Default: 1.
        lock (Optional[Lock]): Thread synchronization lock. Default: None.
        disable_bar (bool): Whether to disable the progress bar. Default: False.

    Returns:
        Tuple[nx.Graph, Dict[int, Set[GeneFamily]], Dict[int, str], Dict[int, Spot]]:
            - spots_graph (nx.Graph): Unified graph containing all spots as nodes
            - spots2borders (Dict[int, Set[GeneFamily]]): Complete mapping of spot
              hashes to their bordering families
            - spots2pangenome (Dict[int, str]): Complete mapping of spot hashes to
              pangenome names
            - spothash2spot (Dict[int, Spot]): Complete mapping of spot hashes to
              spot objects

    Raises:
        RuntimeError: If parallel processing fails or produces inconsistent results.
    """
    # Use ThreadPoolExecutor for parallel processing of pangenomes
    with ThreadPoolExecutor(
        max_workers=threads, initializer=init_lock, initargs=(lock,)
    ) as executor:

        # Set up progress tracking
        with tqdm(
            total=len(pangenomes),
            unit="Pangenome",
            desc="Creating spots graphs",
            disable=disable_bar,
        ) as pbar:

            # Submit tasks for each pangenome
            futures = []
            for pangenome in pangenomes:
                logger.debug(
                    f"Submitting spot processing for pangenome '{pangenome.name}'"
                )
                future = executor.submit(
                    create_pangenome_spots_graph, pangenome, dup_margin
                )
                future.add_done_callback(lambda p: pbar.update())
                futures.append(future)

            # Initialize unified data structures
            spots_graph = nx.Graph()
            spots2borders = {}
            spots2pangenome = {}
            spothash2spot = {}

            # Collect and merge results from all pangenomes
            for future in futures:
                try:
                    pan_graph, pan_borders, pan_mapping, pan_spots = future.result()

                    # Merge graph nodes (with attributes)
                    spots_graph.add_nodes_from(pan_graph.nodes(data=True))

                    # Merge mappings
                    spots2borders.update(pan_borders)
                    spots2pangenome.update(pan_mapping)
                    spothash2spot.update(pan_spots)

                except Exception as e:
                    raise RuntimeError(f"Failed to process pangenome: {e}")

    total_spots = len(spots_graph.nodes)
    unique_pangenomes = len(set(spots2pangenome.values()))

    logger.info(
        f"Successfully created unified spots graph with {total_spots} spots "
        f"from {unique_pangenomes} pangenomes"
    )

    return spots_graph, spots2borders, spots2pangenome, spothash2spot


def compute_gfrr_edges(
    graph: nx.Graph,
    spots2borders: Dict[int, Set[GeneFamily]],
    spots2pangenome: Dict[int, str],
    min_gfrr_cutoff: float = 0.5,
    max_gfrr_cutoff: float = 0.8,
    disable_bar: bool = False,
) -> None:
    """
    Compute and add edges between spots based on Gene Family Relatedness Relationship (GFRR).

    This function calculates GFRR metrics between all pairs of spots from different
    pangenomes and adds edges to the graph when both minimum and maximum GFRR values
    exceed their respective cutoffs. This creates connections between potentially
    conserved genomic spots.

    Args:
        graph (nx.Graph): Spots graph to add edges to (modified in-place).
        spots2borders (Dict[int, Set[GeneFamily]]): Mapping of spot hashes to
                                                   their bordering gene families.
        spots2pangenome (Dict[int, str]): Mapping of spot hashes to pangenome names.
        min_gfrr_cutoff (float): Minimum threshold for min_gfrr metric. Default: 0.5.
        max_gfrr_cutoff (float): Minimum threshold for max_gfrr metric. Default: 0.8.
        disable_bar (bool): Whether to disable the progress bar. Default: False.
    """
    # Generate all possible pairs of spots from different pangenomes
    # This avoids unnecessary intra-pangenome comparisons
    spots_pairs = [
        (spot1, spot2)
        for spot1, spot2 in combinations(graph.nodes, 2)
        if spots2pangenome[spot1] != spots2pangenome[spot2]
    ]

    logger.info(
        f"Computing GFRR edges for {len(spots_pairs)} inter-pangenome spot pairs"
    )

    # Track statistics for reporting
    edges_added = 0
    total_comparisons = len(spots_pairs)

    # Process each spot pair with progress tracking
    with tqdm(
        total=total_comparisons,
        unit="spot pairs",
        desc="Computing GFRR edges",
        disable=disable_bar,
    ) as pbar:

        for spot1, spot2 in spots_pairs:
            # Skip if the edge already exists (shouldn't happen in normal flow)
            if graph.has_edge(spot1, spot2):
                logger.debug(f"Edge already exists between {spot1} and {spot2}")
                pbar.update()
                continue

            # Extract bordering families for both spots
            border1_families = spots2borders[spot1]
            border2_families = spots2borders[spot2]

            # Skip comparison if either spot has no bordering families
            if not border1_families or not border2_families:
                pbar.update()
                continue

            # Compute GFRR metrics between the two spots
            try:
                min_gfrr, max_gfrr, shared_count = compute_gfrr(
                    border1_families, border2_families
                )
            except ValueError as e:
                logger.warning(
                    f"GFRR computation failed for spots {spot1}, {spot2}: {e}"
                )
                pbar.update()
                continue

            # Add edge if both GFRR thresholds are satisfied
            if min_gfrr >= min_gfrr_cutoff and max_gfrr >= max_gfrr_cutoff:
                graph.add_edge(
                    spot1,
                    spot2,
                    min_gfrr=min_gfrr,
                    max_gfrr=max_gfrr,
                    shared_families=shared_count,
                    weight=max_gfrr,  # Use max_gfrr as the default edge weight
                )
                edges_added += 1

                logger.debug(
                    f"Added edge: {spots2pangenome[spot1]} spot {spot1} <-> "
                    f"{spots2pangenome[spot2]} spot {spot2} "
                    f"(min_gfrr={min_gfrr:.3f}, max_gfrr={max_gfrr:.3f})"
                )

            pbar.update()

    logger.info(
        f"Added {edges_added} edges out of {total_comparisons} comparisons "
        f"({edges_added/total_comparisons*100:.1f}% pass thresholds)"
    )


def add_systems_info(pangenomes: Pangenomes, cs_graph: nx.Graph) -> None:
    """
    Annotate conserved spots graph with systems information.

    This function enriches the conserved spots graph by adding system-related
    attributes to nodes. It identifies which systems are associated with each
    conserved spot and adds boolean attributes for system presence as well as
    a count of total systems per spot.

    Args:
        pangenomes (Pangenomes): Collection of pangenomes containing systems data.
        cs_graph (nx.Graph): Conserved spots' graph to annotate (modified in-place).
    """
    # Collect all unique system names across pangenomes
    all_system_names = set()
    node2system_count = defaultdict(int)

    logger.info("Adding systems information to conserved spots graph")

    # Process each pangenome to extract system-spot relationships
    for pangenome in pangenomes:
        logger.debug(f"Processing systems from pangenome '{pangenome.name}'")

        for system in pangenome.systems:
            all_system_names.add(system.name)

            # Find all spots associated with this system through model families
            associated_spots = set()
            for gene_family in system.model_families:
                # Gene families can be associated with multiple spots
                associated_spots.update(gene_family.spots)

            # Update graph nodes for spots associated with this system
            for spot in associated_spots:
                spot_hash = hash((spot.ID, pangenome.name))

                # Only process if the spot exists in the graph
                if cs_graph.has_node(spot_hash):
                    node_attributes = cs_graph.nodes[spot_hash]
                    node_attributes[system.name] = True
                    node2system_count[spot_hash] += 1

                    logger.debug(
                        f"Associated system '{system.name}' with spot {spot.ID} "
                        f"from pangenome '{pangenome.name}'"
                    )

    # Ensure all nodes have all system attributes (set False for missing systems)
    logger.info(
        f"Standardizing {len(all_system_names)} system attributes across all nodes"
    )

    for node_hash in cs_graph.nodes:
        node_attributes = cs_graph.nodes[node_hash]

        # Add missing system attributes as False
        for system_name in all_system_names:
            if system_name not in node_attributes:
                node_attributes[system_name] = False

        # Add total systems count
        node_attributes["#systems"] = node2system_count[node_hash]

    systems_count = len(all_system_names)
    nodes_with_systems = sum(1 for count in node2system_count.values() if count > 0)

    logger.info(
        f"Added {systems_count} system attributes to {len(cs_graph.nodes)} nodes. "
        f"{nodes_with_systems} nodes have associated systems."
    )


def create_pangenome_system_graph(
    pangenome: Pangenome, canonical: bool = False
) -> Tuple[nx.Graph, Dict[int, str], Dict[int, Any], defaultdict]:
    """
    Create a system graph for a single pangenome.

    This function constructs a graph representation of systems within a pangenome,
    where each node represents a unique system-spot combination. It computes
    system coverage statistics and tracks conserved spots associations.

    Args:
        pangenome (Pangenome): Source pangenome containing systems to process.
        canonical (bool): Whether to use canonical systems (True) or not (False). Default: False.

    Returns:
        Tuple containing:
            - graph (nx.Graph): NetworkX graph with system nodes
            - sys2pangenome (Dict[int, str]): Mapping of system hashes to pangenome names
            - syshash2sys (Dict[int, Any]): Mapping of system hashes to system objects
            - systemhash2conserved_spots (defaultdict): Mapping of system hashes to
              sets of conserved spot IDs
    """
    # Initialize data structures
    graph = nx.Graph()
    sys2pangenome = {}
    syshash2sys = {}
    systemhash2conserved_spots = defaultdict(set)
    syshash2organisms = defaultdict(set)  # Track organisms per system-spot combination

    logger.debug(
        f"Processing {pangenome.number_of_systems(with_canonical=canonical)} systems from pangenome '{pangenome.name}'"
    )

    # Process each system in the pangenome
    for system in pangenome.systems:
        # Process each spot associated with this system
        for spot in system.spots:
            # Create a unique hash for a system-spot combination
            system_hash = hash((pangenome.name, system.name, spot.ID))

            # Add a node only once per system-spot combination
            if system_hash not in syshash2sys:
                # Track organisms associated with this system
                syshash2organisms[system_hash].update(system.organisms)

                # Add a node with comprehensive attributes
                graph.add_node(
                    system_hash,
                    system_name=system.name,
                    pangenome=pangenome.name,
                    spot=spot.ID,
                    system_id=system.ID if hasattr(system, "ID") else None,
                    num_model_families=(
                        system.number_of_model_gene_families
                        if hasattr(system, "model_families")
                        else 0
                    ),
                )

                # Store mappings for later use
                sys2pangenome[system_hash] = pangenome.name
                syshash2sys[system_hash] = system

            # Track conserved spots if available
            if hasattr(spot, "conserved_id") and spot.conserved_id:
                systemhash2conserved_spots[system_hash].add(spot.conserved_id)

    # Calculate organism coverage percentages for each system-spot combination
    organism_percentages = {
        sys_hash: (len(organisms) / pangenome.number_of_organisms) * 100
        for sys_hash, organisms in syshash2organisms.items()
    }

    # Set organism percentage as node attribute
    nx.set_node_attributes(graph, organism_percentages, "percent_org")

    logger.debug(
        f"Created system graph with {len(graph.nodes)} system-spot nodes "
        f"from pangenome '{pangenome.name}'"
    )

    return graph, sys2pangenome, syshash2sys, systemhash2conserved_spots


def create_systems_graph(
    pangenomes: Pangenomes,
    canonical: bool = False,
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> Tuple[nx.Graph, Dict[int, str], Dict[int, Any], defaultdict]:
    """
    Create a unified systems graph from multiple pangenomes.

    This function processes all pangenomes in parallel to construct a comprehensive
    graph where nodes represent system-spot combinations across all pangenomes.
    The resulting graph serves as input for systems clustering and conserved spots
    analysis.

    Args:
        pangenomes (Pangenomes): Collection of pangenomes to process.
        canonical (bool): Whether to use canonical systems (True) or not (False). Default: False.
        threads (int): Number of threads for parallel processing. Default: 1.
        lock (Optional[Lock]): Thread synchronization lock. Default: None.
        disable_bar (bool): Whether to disable the progress bar. Default: False.

    Returns:
        Tuple containing:
            - systems_graph (nx.Graph): Unified graph with all system-spot nodes
            - systems2pangenome (Dict[int, str]): Mapping of system hashes to pangenome names
            - systemhash2system (Dict[int, Any]): Mapping of system hashes to system objects
            - systemhash2conserved_spots (defaultdict): Mapping of system hashes to
              sets of conserved spot IDs

    Raises:
        RuntimeError: If parallel processing fails for any pangenome.
    """
    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(
        max_workers=threads, initializer=init_lock, initargs=(lock,)
    ) as executor:

        with tqdm(
            total=len(pangenomes),
            unit="Pangenome",
            desc="Creating systems graphs",
            disable=disable_bar,
        ) as pbar:

            # Submit processing tasks for each pangenome
            futures = []
            for pangenome in pangenomes:
                logger.debug(
                    f"Submitting system processing for pangenome '{pangenome.name}'"
                )
                future = executor.submit(
                    create_pangenome_system_graph, pangenome, canonical
                )
                future.add_done_callback(lambda p: pbar.update())
                futures.append(future)

            # Initialize unified data structures
            systems_graph = nx.Graph()
            systems2pangenome = {}
            systemhash2system = {}
            systemhash2conserved_spots = defaultdict(set)

            # Collect and merge results from all pangenomes
            for future in futures:
                try:
                    pan_graph, pan_sys_mapping, pan_sys_objects, pan_conserved = (
                        future.result()
                    )

                    # Merge graph nodes with their attributes
                    systems_graph.add_nodes_from(pan_graph.nodes(data=True))

                    # Merge all mappings
                    systems2pangenome.update(pan_sys_mapping)
                    systemhash2system.update(pan_sys_objects)
                    systemhash2conserved_spots.update(pan_conserved)

                except Exception as e:
                    raise RuntimeError(f"Failed to process pangenome systems: {e}")

    total_system_nodes = len(systems_graph.nodes)
    unique_pangenomes = len(set(systems2pangenome.values()))

    logger.info(
        f"Created unified systems graph with {total_system_nodes} system-spot nodes "
        f"from {unique_pangenomes} pangenomes"
    )

    return (
        systems_graph,
        systems2pangenome,
        systemhash2system,
        systemhash2conserved_spots,
    )


def graph_systems_link_with_conserved_spots(
    pangenomes: Pangenomes,
    output: Path,
    graph_formats: Optional[List[str]] = None,
    canonical: bool = False,
    community: str = "Louvain",
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> None:
    """
    Generate and analyze systems linkage graphs based on conserved spots.

    This function creates comprehensive graphs linking systems through shared conserved
    spots and performs dual clustering analysis using both Louvain community detection
    or Minimum Spanning Tree (MST) approaches. The analysis identifies system clusters
    that span multiple pangenomes and share conserved genomic locations.

    Args:
        community:
        pangenomes (Pangenomes): Collection of pangenomes with systems and conserved spots.
        output (Path): Output directory for generated graphs and results.
        graph_formats (Optional[List[str]]): Output formats ['gexf', 'graphml']. Default: None.
        canonical: Whether to use the canonical systems (True) or not (False). Default: False.
        threads (int): Number of threads for parallel processing. Default: 1.
        lock (Optional[Lock]): Thread synchronization lock. Default: None.
        disable_bar (bool): Whether to disable progress bars. Default: False.
    """
    logger.info("Starting systems linkage analysis with conserved spots")

    # Create a unified systems graph from all pangenomes
    systems_graph, system2pangenome, systemhash2system, systemhash2conserved_spots = (
        create_systems_graph(pangenomes, canonical, threads, lock, disable_bar)
    )

    if len(systems_graph.nodes) == 0:
        logger.warning(
            "No systems found in pangenomes - skipping systems linkage analysis"
        )
        return

    # Calculate total possible edges for progress tracking
    num_nodes = len(systems_graph.nodes)
    total_combinations = num_nodes * (num_nodes - 1) // 2

    logger.info(
        f"Computing edges between {num_nodes} system nodes ({total_combinations} combinations)"
    )

    # Compute edges based on shared conserved spots
    edges_added = 0
    with tqdm(
        combinations(systems_graph.nodes, 2),
        total=total_combinations,
        desc="Computing system linkages",
        disable=disable_bar,
    ) as pbar:

        for sys_hash1, sys_hash2 in pbar:
            # Find conserved spots shared between these two systems
            common_conserved_spots = (
                systemhash2conserved_spots[sys_hash1]
                & systemhash2conserved_spots[sys_hash2]
            )

            # Add edge if systems share conserved spots
            if common_conserved_spots:
                edge_weight = len(common_conserved_spots)
                systems_graph.add_edge(
                    sys_hash1,
                    sys_hash2,
                    cluster_spots=",".join(map(str, sorted(common_conserved_spots))),
                    weight=edge_weight,
                    num_shared_spots=edge_weight,
                )
                edges_added += 1

                # Enrich node attributes with system and pangenome information
                for sys_hash in [sys_hash1, sys_hash2]:
                    node_attr = systems_graph.nodes[sys_hash]
                    pangenome_name = system2pangenome[sys_hash]
                    system_obj = systemhash2system[sys_hash]

                    # Update node attributes
                    node_attr.update(
                        {"system_name": system_obj.name, "pangenome": pangenome_name}
                    )

                    # Collect associated spot IDs for this system in shared conserved spots
                    associated_spots = set()
                    for cs_id in common_conserved_spots:
                        conserved_spot = pangenomes.get_conserved_spots(cs_id)
                        for spot in conserved_spot.spots:
                            if spot.pangenome.name == pangenome_name:
                                associated_spots.add(spot.ID)

                    # Update or initialize spots attribute
                    if "spots" in node_attr:
                        node_attr["spots"].update(associated_spots)
                    else:
                        node_attr["spots"] = associated_spots.copy()

    logger.info(f"Added {edges_added} edges between systems sharing conserved spots")

    # Remove isolated nodes (systems with no shared conserved spots)
    isolated_nodes = list(nx.isolates(systems_graph))
    if isolated_nodes:
        systems_graph.remove_nodes_from(isolated_nodes)
        logger.info(f"Removed {len(isolated_nodes)} isolated system nodes")

    # Convert spot sets to comma-separated strings for output
    for node in systems_graph.nodes:
        node_attr = systems_graph.nodes[node]
        if "spots" in node_attr and isinstance(node_attr["spots"], set):
            node_attr["spots"] = ",".join(map(str, sorted(node_attr["spots"])))

    # LOUVAIN CLUSTERING ANALYSIS
    if community == "Louvain":
        logger.info("Performing Louvain community detection")
        louvain_graph = systems_graph.copy()
        partitions = nx.algorithms.community.louvain_communities(
            louvain_graph, weight="weight"
        )

        cluster_id = 1
        clusters_retained = 0

        for cluster_systems in partitions:
            if len(cluster_systems) > 1:  # Only process multi-node clusters
                # Check if the cluster spans multiple pangenomes
                pangenomes_in_cluster = {
                    system2pangenome[sys_hash] for sys_hash in cluster_systems
                }

                if len(pangenomes_in_cluster) > 1:
                    # Multi-pangenome cluster - retain and label
                    for sys_hash in cluster_systems:
                        node_attr = louvain_graph.nodes[sys_hash]
                        system_obj = systemhash2system[sys_hash]
                        pangenome_name = system2pangenome[sys_hash]

                        node_attr.update(
                            {
                                "system_ID": getattr(system_obj, "ID", None),
                                "system_name": system_obj.name,
                                "pangenome": pangenome_name,
                                "cluster_id": cluster_id,
                            }
                        )

                    cluster_id += 1
                    clusters_retained += 1
                else:
                    # Single-pangenome cluster - remove
                    louvain_graph.remove_nodes_from(cluster_systems)
            else:
                # Single-node cluster - remove
                louvain_graph.remove_nodes_from(cluster_systems)

        logger.info(
            f"Louvain clustering: retained {clusters_retained} inter-pangenome clusters"
        )

        # Save Louvain clustering results
        if graph_formats is not None:
            for fmt in graph_formats:
                if fmt == "gexf":
                    graph_file = output / "systems_link_with_conserved_spots_louvain.gexf"
                    logger.info(f"Writing Louvain graph in GEXF format: {graph_file}")
                    nx.readwrite.gexf.write_gexf(louvain_graph, graph_file)
                elif fmt == "graphml":
                    graph_file = (
                        output / "systems_link_with_conserved_spots_louvain.graphml"
                    )
                    logger.info(f"Writing Louvain graph in GraphML format: {graph_file}")
                    nx.readwrite.graphml.write_graphml(louvain_graph, graph_file)

    if community == "MST":
        # MINIMUM SPANNING TREE (MST) CLUSTERING ANALYSIS
        logger.info("Performing MST-based clustering with automatic threshold detection")

        if len(systems_graph.edges) == 0:
            logger.warning("No edges in systems graph - skipping MST analysis")
            return

        # Build Minimum Spanning Tree
        mst = nx.minimum_spanning_tree(systems_graph, weight="weight")

        # Extract edge weights for threshold analysis
        edge_weights = np.array([d["weight"] for _, _, d in mst.edges(data=True)])

        if len(edge_weights) == 0:
            logger.warning("MST has no edges - skipping threshold analysis")
            return

        # Automatic threshold detection using weight distribution analysis
        sorted_weights = np.sort(edge_weights)

        if len(sorted_weights) > 1:
            # Find the largest jump in consecutive weights
            weight_differences = np.diff(sorted_weights)
            max_jump_index = np.argmax(weight_differences)
            threshold = sorted_weights[max_jump_index]

            logger.info(
                f"MST threshold detection: identified threshold = {threshold} "
                f"(max jump at index {max_jump_index})"
            )
        else:
            # Fallback for the single edge case
            threshold = sorted_weights[0]
            logger.info(f"MST threshold detection: single edge, using weight = {threshold}")

        # Remove edges with weights above the threshold
        heavy_edges = [
            (u, v) for u, v, d in mst.edges(data=True) if d["weight"] > threshold
        ]
        if heavy_edges:
            mst.remove_edges_from(heavy_edges)
            logger.info(f"Removed {len(heavy_edges)} edges above threshold from MST")

        # Identify connected components as clusters
        mst_clusters_retained = 0
        for cluster_id, cluster_systems in enumerate(nx.connected_components(mst), start=1):
            if len(cluster_systems) > 1:
                # Assign cluster labels to nodes
                for sys_hash in cluster_systems:
                    node_attr = mst.nodes[sys_hash]
                    system_obj = systemhash2system[sys_hash]
                    pangenome_name = system2pangenome[sys_hash]

                    node_attr.update(
                        {
                            "system_name": system_obj.name,
                            "pangenome": pangenome_name,
                            "cluster_id": cluster_id,
                        }
                    )

                mst_clusters_retained += 1
            else:
                # Remove single-node components
                mst.remove_nodes_from(cluster_systems)

        logger.info(f"MST clustering: identified {mst_clusters_retained} clusters")

        # Save MST clustering results
        if graph_formats is not None:
            for fmt in graph_formats:
                if fmt == "gexf":
                    graph_file = output / "systems_link_with_conserved_spots_mst.gexf"
                    logger.info(f"Writing MST graph in GEXF format: {graph_file}")
                    nx.readwrite.gexf.write_gexf(mst, graph_file)
                elif fmt == "graphml":
                    graph_file = output / "systems_link_with_conserved_spots_mst.graphml"
                    logger.info(f"Writing MST graph in GraphML format: {graph_file}")
                    nx.readwrite.graphml.write_graphml(mst, graph_file)


def write_conserved_spots(
    pangenomes: Pangenomes,
    output: Path,
    graph_formats: Optional[List[str]] = None,
    cs_graph: Optional[nx.Graph] = None,
    force: bool = False,
    disable_bar: bool = False,
) -> None:
    """
    Write conserved spots data to files and optionally export graphs.

    This function generates comprehensive output files documenting conserved spots
    across pangenomes, including individual spot details and summary statistics.
    It also integrates systems information and exports graph representations when requested.

    Args:
        pangenomes (Pangenomes): Collection of pangenomes with conserved spots.
        output (Path): Output directory for generated files.
        graph_formats (Optional[List[str]]): Graph export formats ['gexf', 'graphml']. Default: None.
        cs_graph (Optional[nx.Graph]): Conserved spots' graph to export. Default: None.
        force (bool): Whether to overwrite existing files. Default: False.
        disable_bar (bool): Whether to disable progress bars. Default: False.

    Side Effects:
        - Creates output directory structure
        - Generates individual TSV files for each conserved spot
        - Creates a summary TSV file with all conserved spots
        - Exports graph files in specified formats
        - Modifies cs_graph by adding systems information

    Raises:
        AssertionError: If graph_formats specified but cs_graph is None.
        IOError: If file writing fails.
    """
    logger.info("Writing conserved spots data to files")

    # Create output directory structure
    cs_dir = mkdir(output / "conserved_spots", force=force, erase=force)

    # Data collection for a summary file
    all_conserved_spots_data = []

    # Process each conserved spot
    with tqdm(
        pangenomes.conserved_spots,
        total=pangenomes.number_of_conserved_spots,
        desc="Processing conserved spots",
        disable=disable_bar,
    ) as pbar:

        for conserved_spot in pbar:
            # Data for an individual conserved spot file
            individual_spot_data = []

            # Process each spot within this conserved spot
            for spot in conserved_spot.spots:
                # Add to summary data
                all_conserved_spots_data.append(
                    [
                        conserved_spot.ID,
                        spot.ID,
                        spot.pangenome.name,
                        len(spot),  # Number of RGPs in spot
                        spot.number_of_families,
                    ]
                )

                # Add detailed information for an individual file
                for rgp in spot.regions:
                    family_names = [family.name for family in rgp.families]
                    individual_spot_data.append(
                        [spot.ID, spot.pangenome.name, rgp.name, ",".join(family_names)]
                    )

            # Write an individual conserved spot file
            if individual_spot_data:
                individual_df = pd.DataFrame(
                    individual_spot_data,
                    columns=["Spot_ID", "Pangenome", "RGP_Name", "Gene_Families"],
                )
                individual_df = individual_df.sort_values(
                    ["Spot_ID", "Pangenome", "RGP_Name"]
                )

                individual_file = cs_dir / f"conserved_spot_{conserved_spot.ID}.tsv"
                individual_df.to_csv(
                    individual_file, sep="\t", header=True, index=False
                )

                logger.debug(f"Written individual file: {individual_file}")

    # Write a summary file with all conserved spots
    if all_conserved_spots_data:
        summary_df = pd.DataFrame(
            all_conserved_spots_data,
            columns=[
                "Conserved_Spot_ID",
                "Spot_ID",
                "Pangenome",
                "Num_RGPs",
                "Num_Gene_Families",
            ],
        )
        summary_df = summary_df.sort_values(
            [
                "Conserved_Spot_ID",
                "Spot_ID",
                "Pangenome",
                "Num_RGPs",
                "Num_Gene_Families",
            ]
        )

        summary_file = output / "all_conserved_spots.tsv"
        summary_df.to_csv(summary_file, sep="\t", header=True, index=False)

        logger.info(
            f"Written summary file: {summary_file} "
            f"({len(all_conserved_spots_data)} spot entries)"
        )

    # Add systems information to the graph if available
    if cs_graph is not None:
        logger.info("Integrating systems information into conserved spots graph")
        add_systems_info(pangenomes, cs_graph)

    # Export graph files if requested
    if graph_formats is not None:
        if cs_graph is None:
            raise AssertionError(
                "Graph formats specified but no conserved spots graph provided. "
                "Pass cs_graph parameter to enable graph export."
            )

        logger.info(
            f"Exporting conserved spots graph in {len(graph_formats)} format(s)"
        )

        for fmt in graph_formats:
            if fmt == "gexf":
                graph_file = output / "conserved_spots.gexf"
                logger.info(
                    f"Writing conserved spots graph in GEXF format: {graph_file}"
                )
                try:
                    nx.readwrite.gexf.write_gexf(cs_graph, graph_file)
                except Exception as e:
                    logger.error(f"Failed to write GEXF file: {e}")

            elif fmt == "graphml":
                graph_file = output / "conserved_spots.graphml"
                logger.info(
                    f"Writing conserved spots graph in GraphML format: {graph_file}"
                )
                try:
                    nx.readwrite.graphml.write_graphml(cs_graph, graph_file)
                except Exception as e:
                    logger.error(f"Failed to write GraphML file: {e}")
            else:
                logger.warning(f"Unsupported graph format: {fmt}")


def compare_spots(
    pangenomes: Pangenomes,
    dup_margin: float = 0.05,
    gfrr_metrics: str = "min_gfrr",
    gfrr_cutoff: Tuple[float, float] = (0.8, 0.8),
    threads: int = 1,
    lock: Optional[Lock] = None,
    disable_bar: bool = False,
) -> nx.Graph:
    """
    Identify and cluster conserved spots across multiple pangenomes.

    This is the main function for conserved spots detection. It creates a comprehensive
    graph of spots from all pangenomes, computes similarity edges based on Gene Family
    Relatedness Relationship (GFRR) metrics, performs clustering to identify conserved
    spots, and integrates the results into the pangenomes object.

    Args:
        pangenomes (Pangenomes): Collection of pangenomes to analyze.
        dup_margin (float): Minimum ratio for multigenic family detection. Default: 0.05.
        gfrr_metrics (str): GFRR metric for clustering ('min_gfrr' or 'max_gfrr'). Default: 'min_gfrr'.
        gfrr_cutoff (Tuple[float, float]): Thresholds for (min_gfrr, max_gfrr). Default: (0.8, 0.8).
        threads (int): Number of threads for parallel processing. Default: 1.
        lock (Optional[Lock]): Thread synchronization lock. Default: None.
        disable_bar (bool): Whether to disable progress bars. Default: False.

    Returns:
        nx.Graph: Final spots graph with clustering information and conserved spots annotations.

    Side Effects:
        - Adds ConservedSpots objects to the pangenomes collection
        - Modifies the returned graph by adding cluster assignments and removing isolated nodes

    Raises:
        ValueError: If gfrr_metrics is not 'min_gfrr' or 'max_gfrr'.
        RuntimeError: If clustering fails or produces no valid clusters.
    """
    if gfrr_metrics not in ["min_gfrr", "max_gfrr"]:
        raise ValueError(
            f"Invalid gfrr_metrics: {gfrr_metrics}. Must be 'min_gfrr' or 'max_gfrr'"
        )

    logger.info(
        f"Starting conserved spots comparison with {len(pangenomes)} pangenomes "
        f"using {gfrr_metrics} metric and cutoffs {gfrr_cutoff}"
    )

    # Step 1: Create a comprehensive spots graph from all pangenomes
    logger.info("Step 1: Creating spots graph from all pangenomes")
    spots_graph, spots2borders, spots2pangenome, spothash2spot = create_spots_graph(
        pangenomes, dup_margin, threads, lock, disable_bar
    )

    initial_nodes = len(spots_graph.nodes)
    logger.info(f"Created initial spots graph with {initial_nodes} nodes")

    # Step 2: Compute GFRR-based edges between spots
    logger.info("Step 2: Computing GFRR-based similarity edges")
    compute_gfrr_edges(spots_graph, spots2borders, spots2pangenome, min_gfrr_cutoff=gfrr_cutoff[0],
                       max_gfrr_cutoff=gfrr_cutoff[1], disable_bar=disable_bar)

    edges_count = len(spots_graph.edges)
    logger.info(f"Added {edges_count} similarity edges to spots graph")

    # Step 3: Perform clustering to identify conserved spots
    logger.info(f"Step 3: Clustering spots using {gfrr_metrics} metric")
    partitions = cluster_on_gfrr(spots_graph, gfrr_metrics)

    # Step 4: Process clusters and create ConservedSpots objects
    logger.info("Step 4: Processing clusters and creating conserved spots")
    conserved_spots_created = 0
    nodes_in_clusters = 0

    for cs_id, cluster_spots in enumerate(partitions, start=1):
        if len(cluster_spots) > 1:  # Only process multi-spot clusters
            # Collect actual spot objects for this cluster
            conserved_spot_members = set()

            for spot_hash in cluster_spots:
                spot = spothash2spot[spot_hash]
                pangenome_name = spots2pangenome[spot_hash]

                # Update node attributes with clustering information
                node_attributes = spots_graph.nodes[spot_hash]

                # Remove the temporary clustering attribute
                cluster_attr_name = f"{gfrr_metrics}_cluster"
                if cluster_attr_name in node_attributes:
                    del node_attributes[cluster_attr_name]

                # Add final clustering information
                node_attributes.update(
                    {
                        "conserved_spot_id": cs_id,
                        "spot_id": spot.ID,
                        "pangenome": pangenome_name,
                        "num_families": spot.number_of_families,
                        "num_regions": len(spot),
                    }
                )

                conserved_spot_members.add(spot)
                nodes_in_clusters += 1

            # Create and add ConservedSpots object to pangenomes
            conserved_spot = ConservedSpots(cs_id, *conserved_spot_members)
            pangenomes.add_conserved_spots(conserved_spot)
            conserved_spots_created += 1

            logger.debug(
                f"Created conserved spot {cs_id} with {len(conserved_spot_members)} spots "
                f"spanning {len(set(s.pangenome.name for s in conserved_spot_members))} pangenomes"
            )

        else:
            # Remove single-spot clusters (not conserved)
            spots_graph.remove_nodes_from(cluster_spots)

    # Final statistics and logging
    final_nodes = len(spots_graph.nodes)
    removed_nodes = initial_nodes - final_nodes

    logger.info(
        f"Conserved spots analysis completed:\n"
        f"  - Initial spots: {initial_nodes}\n"
        f"  - Similarity edges: {edges_count}\n"
        f"  - Conserved spots created: {conserved_spots_created}\n"
        f"  - Spots in conserved groups: {nodes_in_clusters}\n"
        f"  - Isolated spots removed: {removed_nodes}\n"
        f"  - Final graph nodes: {final_nodes}"
    )

    if conserved_spots_created == 0:
        logger.warning(
            "No conserved spots were identified. Consider:\n"
            "  - Lowering GFRR cutoff thresholds\n"
            "  - Checking spot quality and distribution\n"
            "  - Verifying pangenome compatibility"
        )

    return spots_graph


def launch(args: argparse.Namespace) -> None:
    """
    Main entry point for conserved spots comparison analysis.

    This function orchestrates the complete workflow for identifying and analyzing
    conserved spots across multiple pangenomes. It handles argument validation,
    resource setup, analysis execution, and results output.

    Args:
        args (argparse.Namespace): Parsed command-line arguments containing all
                                 configuration parameters for the analysis.

    Raises:
        Various exceptions from underlying functions for validation, processing, or I/O errors.
    todo:
        - add a community argument
    """
    logger.info("Starting conserved spots comparison analysis")

    # Step 1: Validate arguments and determine data requirements
    try:
        need_info = check_compare_spots_args(args)
        logger.info("Arguments validation successful")
    except argparse.ArgumentError as e:
        logger.error(f"Argument validation failed: {e}")
        raise

    # Step 2: Load pangenomes and set up resources
    logger.info("Loading pangenomes and setting up analysis resources")
    pangenomes, tmpdir, _, lock = common_launch(args, check_pangenome_cs, need_info)

    # Step 3: Create an output directory
    output = mkdir(args.output, force=args.force)
    logger.info(f"Analysis results will be written to: {output}")

    # Step 4: Perform conserved spots analysis
    logger.info("Performing conserved spots identification and clustering")
    spots_graph = compare_spots(
        pangenomes=pangenomes,
        dup_margin=args.dup_margin,
        gfrr_metrics=args.gfrr_metrics,
        gfrr_cutoff=args.gfrr_cutoff,
        threads=args.cpus,
        lock=lock,
        disable_bar=args.disable_prog_bar
    )

    # Step 5: Write conserved spots results
    logger.info("Writing conserved spots results to files")
    write_conserved_spots(
        pangenomes,
        output,
        cs_graph=spots_graph,
        graph_formats=args.graph_formats,
        force=args.force,
        disable_bar=args.disable_prog_bar
    )

    # Step 6: Perform systems linkage analysis if requested
    if args.systems:
        logger.info("Performing systems linkage analysis with conserved spots")
        graph_systems_link_with_conserved_spots(pangenomes=pangenomes, output=output, graph_formats=args.graph_formats,
                                                threads=args.cpus, lock=lock, disable_bar=args.disable_prog_bar)
    else:
        logger.info("Systems analysis not requested - skipping systems linkage")

    # Step 7: Clean up temporary files
    if not args.keep_tmp:
        logger.info(f"Cleaning up temporary directory: {tmpdir}")
        rmtree(tmpdir, ignore_errors=True)
    else:
        logger.info(f"Temporary files preserved in: {tmpdir}")

    logger.info("Conserved spots comparison analysis completed successfully")


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Create a subparser for conserved spots comparison command.

    This function configures the command-line interface for the conserved spots
    comparison functionality, setting up the argument parser with appropriate
    description and configuration.

    Args:
        sub_parser (argparse._SubParsersAction): Parent subparser to add this command to.

    Returns:
        argparse.ArgumentParser: Configured parser for the compare_spots command.

    Example:
        >>> main_parser = argparse.ArgumentParser()
        >>> subparsers = main_parser.add_subparsers()
        >>> compare_parser = subparser(subparsers)
    """
    parser = sub_parser.add_parser(
        "compare_spots",
        description="Compare and identify conserved spots across multiple pangenomes. "
                    "This analysis identifies genomic regions that are conserved "
                    "across different pangenomes based on gene family similarity "
                    "and optionally analyzes systems relationships within these regions.",
        help="Compare spots across pangenomes to identify conserved regions"
    )

    # Configure all arguments for this subcommand
    parser_comparison_spots(parser)
    return parser


def parser_comparison_spots(parser: argparse.ArgumentParser) -> None:
    """
    Configure argument parser for conserved spots comparison command.

    This function adds all necessary command-line arguments for conserved spots
    comparison, including core comparison parameters, systems analysis options,
    and various output configurations.

    Args:
        parser (argparse.ArgumentParser): Parser to configure with comparison arguments.
    """
    # Get base comparison arguments (required, compare_opt, optional)
    _, compare_opt, optional = parser_comparison(parser)

    # Add spots-specific comparison options
    compare_opt.add_argument(
        '--gfrr_metrics',
        required=False,
        type=str,
        default="min_gfrr",
        choices=["min_gfrr", "max_gfrr"],
        help="GFRR metric used for spots clustering. "
             "'min_gfrr': conservative metric (shared/smaller_set), "
             "'max_gfrr': liberal metric (shared/larger_set). "
             "Default: min_gfrr"
    )

    # Add general optional arguments
    optional.add_argument(
        "--dup_margin",
        required=False,
        type=float,
        default=0.05,
        help="Minimum ratio of genomes in which a gene family must have "
             "multiple copies to be considered 'duplicated'. This affects "
             "multigenic family detection for spot border analysis. "
             "Range: 0.0-1.0. Default: 0.05 (5%%)"
    )

    # Systems analysis argument group
    systems = parser.add_argument_group(
        title="Systems analysis options",
        description="Optional analysis of systems relationships within conserved spots"
    )

    systems.add_argument(
        '--systems',
        required=False,
        action='store_true',
        default=False,
        help="Enable systems analysis to examine relationships between "
             "conserved spots and detected biological systems. This adds "
             "systems linkage graphs and enriched annotations to the output."
    )

    systems.add_argument(
        '-m', '--models',
        required=False,
        type=Path,
        nargs="+",
        default=None,
        metavar="MODEL_FILE",
        help="Path(s) to system model files. Multiple model files can be "
             "specified (space-separated) for different system sources. "
             "Must be provided in the same order as --sources. "
             "Required if --systems is used."
    )

    systems.add_argument(
        "-s", "--sources",
        required=False,
        type=str,
        nargs="+",
        default=None,
        metavar="SOURCE_NAME",
        help="Name(s) of systems sources corresponding to model files. "
             "Multiple sources can be specified (space-separated). "
             "Must be provided in the same order as --models. "
             "Required if --systems is used. "
    )

    systems.add_argument(
        "--canonical",
        required=False,
        action="store_true",
        default=False,
        help="Include canonical versions of systems in the analysis. "
             "This provides additional system representations that may "
             "be useful for comprehensive systems analysis."
    )
