#!/usr/bin/env python3
# coding:utf-8

"""
A collection of utility functions for pangenome comparison.

This module contains functionalities to compute gene family relatedness, cluster gene families
based on graph metrics, and launch comparative workflows for pangenome analysis. It also
includes utilities for building argument parsers related to these workflows.
"""

# default libraries
from __future__ import annotations
import logging
import tempfile
from pathlib import Path
from typing import Any, Dict, Callable, List, Set, Tuple
from multiprocessing import Manager, Lock

# installed libraries
import networkx as nx

# local libraries
from panorama.pangenomes import Pangenomes
from panorama.geneFamily import GeneFamily
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.cluster import (
    cluster_gene_families,
    write_clustering,
    parser_mmseqs2_cluster,
)


def compute_gfrr(
    queries: Set[GeneFamily], targets: Set[GeneFamily]
) -> Tuple[float, float, int]:
    """
    Compute Gene Family Repertoire Relatedness (GFRR) metrics between query and target gene families.

    This function evaluates the overlap of 'akin' relationships between two sets of gene families
    and calculates both minimum and maximum Family Relatedness Relationship (FRR) values. The FRR
    is a metric to assess similarity between gene family sets based on their shared relationships.

    Args:
        queries (Set[GeneFamily]): Set of query gene families to analyze.
        targets (Set[GeneFamily]): Set of target gene families to compare against.

    Returns:
        Tuple[float, float, int]: A tuple containing:
            - min_frr (float): Minimum FRR value (shared akins / min set size)
            - max_frr (float): Maximum FRR value (shared akins / max set size)
            - num_reciprocal (int): Number of reciprocal akin relationships found

    Raises:
        ValueError: If either queries or targets set are empty.

    Note:
        - min_frr represents conservative similarity (harder to achieve high values)
        - max_frr represents liberal similarity (easier to achieve high values)
        - Akin relationships represent evolutionary or functional similarities between gene families
    """
    if not queries or not targets:
        raise ValueError("Both queries and targets sets must be non-empty")

    # Find all unique akin IDs that are shared between query and target gene families
    # This represents the "common ground" between the two sets
    shared_akins = {
        query_gf.akin.ID
        for query_gf in queries
        for target_gf in targets
        if query_gf.akin == target_gf.akin  # Check if akin relationships match
    }

    # Calculate conservative similarity metric (normalized by a smaller set)
    min_frr = len(shared_akins) / min(len(queries), len(targets))

    # Calculate liberal similarity metric (normalized by a larger set)
    max_frr = len(shared_akins) / max(len(queries), len(targets))

    return min_frr, max_frr, len(shared_akins)


def cluster_on_gfrr(graph: nx.Graph, gfrr_metric: str) -> List[Set[Any]]:
    """
    Cluster graph nodes using Louvain community detection based on GFRR metrics.

    This function applies the Louvain algorithm for community detection to partition the input
    graph into clusters. Each node is assigned a cluster identifier as a node attribute, and
    the function returns the cluster partitions.

    Args:
        graph (nx.Graph): Input graph with nodes and weighted edges to cluster.
        gfrr_metric (str): Name of the edge weight attribute to use for clustering
                          (e.g., 'min_frr', 'max_frr').

    Returns:
        List[Set[Any]]: List where each element is a set containing nodes belonging
                       to the same cluster.

    Raises:
        KeyError: If the specified gfrr_metric is not found in graph edge attributes.

    TODO:
        - Make the resolution parameter configurable
    """
    # Apply Louvain community detection algorithm
    # This algorithm optimizes modularity to find densely connected communities
    partitions = nx.algorithms.community.louvain_communities(
        graph, weight=gfrr_metric, resolution=1.0  # Default resolution for modularity
    )

    # Assign cluster labels to each node as graph attributes
    # This enables easy retrieval of cluster assignments later
    cluster_attribute_name = f"{gfrr_metric}_cluster"
    for cluster_idx, cluster_nodes in enumerate(partitions):
        cluster_label = f"cluster_{cluster_idx}"
        # Create an attribute dictionary for all nodes in this cluster
        node_attributes = dict.fromkeys(cluster_nodes, cluster_label)
        # Set the attributes on the graph
        nx.set_node_attributes(graph, node_attributes, name=cluster_attribute_name)

    # Log clustering results for monitoring and debugging
    logging.getLogger("PANORAMA").info(
        f"Graph clustered into {len(partitions)} communities using '{gfrr_metric}' metric"
    )

    return partitions


def common_launch(
    args: Any, check_func: Callable, need_info: Dict[str, Any], **kwargs: Any
) -> Tuple[Pangenomes, Path, Manager, Lock]:
    """
    Launch a common setup for comparative pangenome workflows.

    This function handles the common initialization steps for pangenome comparison workflows,
    including loading pangenomes, setting up multiprocessing resources, creating temporary
    directories, and optionally performing gene family clustering.

    Args:
        args (Any): Parsed command-line arguments containing workflow parameters such as
                   - pangenomes: Path to pangenome list file
                   - cpus: Number of CPU threads to use
                   - cluster: Optional path to existing clustering results
                   - tmpdir: Temporary directory path
                   - Various MMSeqs2 clustering parameters
        check_func (Callable): Validation function to check pangenome integrity during loading.
        need_info (Dict[str, Any]): Dictionary specifying required data components:
                                   - 'need_families_sequences': bool, whether sequences are needed
                                   - Other data requirements as key-value pairs
        **kwargs (Any): Additional keyword arguments passed to pangenome loading function.

    Returns:
        Tuple[Pangenomes, Path, Manager, Lock]: Setup resources containing:
            - pangenomes (Pangenomes): Loaded and processed pangenomes object
            - tmpdir (Path): Path to temporary directory for intermediate files
            - manager (Manager): Multiprocessing manager for shared resources
            - lock (Lock): Thread-safe lock for concurrent operations

    Raises:
        FileNotFoundError: If pangenome files or clustering file cannot be found.
        ValueError: If clustering parameters are invalid.
    """
    # Initialize multiprocessing resources for thread-safe operations
    manager = Manager()
    lock = manager.Lock()

    # Determine data requirements based on clustering availability
    # If no clustering is provided, we need sequences to perform clustering
    if args.cluster is None:
        need_info["need_families_sequences"] = True
        logging.getLogger("PANORAMA").info(
            "No clustering provided - will perform gene family clustering"
        )
    else:
        logging.getLogger("PANORAMA").info(
            f"Using existing clustering from: {args.cluster}"
        )

    # Load pangenomes with specified requirements and validation
    logging.getLogger("PANORAMA").info("Loading pangenomes...")
    pangenomes = load_pangenomes(
        pangenome_list=args.pangenomes,
        need_info=need_info,
        check_function=check_func,
        max_workers=args.cpus,
        lock=lock,
        disable_bar=args.disable_prog_bar,
        **kwargs,
    )

    # Create a temporary directory for intermediate processing files
    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    logging.getLogger("PANORAMA").info(f"Created temporary directory: {tmpdir}")

    # Handle gene family clustering workflow
    if args.cluster is None:
        # Configure MMSeqs2 clustering parameters from command-line arguments
        mmseqs2_options = {
            "max_seqs": args.max_seqs,  # Maximum sequences per query
            "min_ungapped": args.min_ungapped,  # Minimum ungapped alignment score
            "comp_bias_corr": args.comp_bias_corr,  # Composition bias correction
            "sensitivity": args.sensitivity,  # Search sensitivity level
            "kmer_per_seq": args.kmer_per_seq,  # K-mers per sequence
            "identity": args.clust_identity,  # Sequence identity threshold
            "coverage": args.clust_coverage,  # Coverage threshold
            "cov_mode": args.clust_cov_mode,  # Coverage calculation mode
            "eval": args.eval,  # E-value threshold
            "max_seq_len": args.max_seq_len,  # Maximum sequence length
            "max_reject": args.max_reject,  # Maximum rejections per query
            "align_mode": args.align_mode,  # Alignment mode
            "clust_mode": args.clust_mode,  # Clustering mode
            "reassign": args.reassign,  # Whether to reassign sequences
        }

        logging.getLogger("PANORAMA").info(
            f"Starting gene family clustering with method: {args.method}"
        )
        # Perform gene family clustering using MMSeqs2
        clustering_results = cluster_gene_families(
            pangenomes=pangenomes,
            method=args.method,
            mmseqs2_opt=mmseqs2_options,
            tmpdir=tmpdir,
            keep_tmp=args.keep_tmp,
            threads=args.cpus,
            lock=lock,
            disable_bar=args.disable_prog_bar,
        )

        # Write clustering results to a file for future use
        cluster_file = tmpdir / "pangenome_gf_clustering.tsv"
        write_clustering(clustering_results, cluster_file)
        logging.getLogger("PANORAMA").info(
            f"Clustering results written to: {cluster_file}"
        )

    else:
        # Use existing clustering file
        cluster_file = args.cluster

    # Load clustering results into a pangenomes object
    logging.getLogger("PANORAMA").info("Loading clustering results into pangenomes...")
    pangenomes.read_clustering(cluster_file, args.disable_prog_bar)

    return pangenomes, tmpdir, manager, lock


def parser_comparison(parser: Any) -> Tuple[Any, Any, Any]:
    """
    Configure an argument parser for pangenome comparison commands.

    This function sets up command-line argument groups and options specific to pangenome
    comparison workflows, including required arguments, comparison options, and MMSeqs2
    clustering parameters.

    Args:
        parser (Any): ArgumentParser instance to configure with comparison-specific arguments.

    Returns:
        Tuple[Any, Any, Any]: Tuple containing three argument groups:
            - required (ArgumentGroup): Required arguments group
            - compare_opt (ArgumentGroup): Comparison optional arguments group
            - optional (ArgumentGroup): General optional arguments group

    Side Effects:
        - Modifies the input parser by adding argument groups and options
        - Configures MMSeqs2 clustering arguments through parser_mmseqs2_cluster

    Note:
        This function defines the command-line interface for comparison workflows,
        making it easier to maintain consistent argument handling across different
        comparison subcommands.
    """
    # Required arguments that must be provided by the user
    required = parser.add_argument_group(
        title="Required arguments",
        description="All of the following arguments are required:",
    )
    required.add_argument(
        "-p",
        "--pangenomes",
        required=True,
        type=Path,
        help="Path to TSV file containing list of pangenome .h5 files to compare",
    )
    required.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output directory where result files will be written",
    )

    # Optional arguments specific to comparison workflows
    compare_opt = parser.add_argument_group(title="Comparison optional arguments")

    compare_opt.add_argument(
        "--cluster",
        required=False,
        type=Path,
        default=None,
        help="Path to tab-separated file with pre-computed clustering results "
        "(cluster_name\\tfamily_id format). If not provided, clustering will be performed.",
    )

    compare_opt.add_argument(
        "--gfrr_cutoff",
        required=False,
        type=float,
        nargs=2,
        default=(0.5, 0.8),
        metavar=("MIN_FRR", "MAX_FRR"),
        help="FRR (Family Relatedness Relationship) cutoff values for similarity assessment. "
        "min_gfrr = shared_families / min(families1, families2), "
        "max_gfrr = shared_families / max(families1, families2)"
        "Default: 0.5 0.8",
    )

    # Configure MMSeqs2 clustering arguments
    # This adds a comprehensive set of clustering parameters
    cluster_group = parser_mmseqs2_cluster(parser)
    cluster_group.description = (
        "MMSeqs2 clustering arguments (used only if --cluster is not provided)"
    )

    cluster_group.add_argument(
        "--method",
        required=False,
        type=str,
        choices=["linclust", "cluster"],
        default="linclust",
        help="MMSeqs2 clustering method selection: "
        "'linclust' - fast linear-time clustering (less sensitive), "
        "'cluster' - slower but more sensitive clustering. "
        "Default: linclust",
    )

    # General optional arguments for workflow control
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument(
        "--graph_formats",
        required=False,
        type=str,
        choices=["gexf", "graphml"],
        nargs="+",
        default=None,
        help="Output format(s) for graph files. Multiple formats can be specified. "
        "Supported: gexf (Gephi Exchange Format), graphml (Graph Markup Language)",
    )

    optional.add_argument(
        "--tmpdir",
        required=False,
        type=Path,
        default=Path(tempfile.gettempdir()),
        help=f"Directory for temporary files. Default: {tempfile.gettempdir()}",
    )

    optional.add_argument(
        "--keep_tmp",
        required=False,
        default=False,
        action="store_true",
        help="Keep temporary files after completion (useful for debugging and inspection)",
    )

    optional.add_argument(
        "-c",
        "--cpus",
        required=False,
        type=int,
        default=1,
        help="Number of CPU threads to use for parallel processing. Default: 1",
    )

    return required, compare_opt, optional
