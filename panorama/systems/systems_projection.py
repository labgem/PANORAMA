#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to project systems onto genomes.
"""

# default libraries
from __future__ import annotations

import itertools
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import logging
from typing import Dict, List, Set, Tuple
from multiprocessing import Lock, get_context
from pathlib import Path
import time

# installed libraries
from tqdm import tqdm
import pandas as pd
import networkx as nx
from ppanggolin.genome import Gene
from ppanggolin.utils import extract_contig_window


# local libraries
from panorama.utils import mkdir, init_lock, conciliate_partition
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.systems.utils import (
    get_metadata_to_families,
    dict_families_context,
    get_gfs_matrix_combination,
    check_needed_families,
)
from panorama.systems.system import System, SystemUnit
from panorama.systems.models import Family


# TODO: if genes are within one connected component => by definition, these requirements are met
# -> this function seems unnecessary
def has_short_path(graph: nx.Graph, node_list: List[GeneFamily], n: int) -> bool:
    """
    Checks if there exists at least one path of length less than `n`
    connecting any two nodes in the given list of nodes in the graph.

    Args:
        graph (nx.Graph): the graph to search paths
        node_list (List[GeneFamily]): List of gene families to check for paths.
        n (int): The maximum length of the path to consider.

    Returns:
        bool: True if there exists at least one path of length less than `n`
              connecting any two nodes in the list, False otherwise.
    """
    path_length = defaultdict(dict)
    has_path = {node: False for node in node_list}
    for i, node1 in enumerate(node_list):
        for node2 in node_list[i + 1 :]:
            if not has_path[node2]:
                try:
                    path_length[node1][node2] = nx.shortest_path_length(
                        graph, source=node1, target=node2
                    )
                    if path_length[node1][node2] <= n:
                        has_path[node1] = True
                        has_path[node2] = True
                        break
                except nx.NetworkXNoPath:
                    continue
    return all(has_path.values())


def project_unit_on_organisms(
    components: List[List[Gene]],
    unit: SystemUnit,
    model_genes: Set[Gene],
    association: List[str] = None,
) -> List[List[str]]:
    """
    Projects a system unit onto a given organism's pangenome.

    Args:
        components (List[List[Gene]]): List of gene components to project.
        unit (SystemUnit): The unit to be projected.
        model_genes (Set[Gene]): Set of genes in one organism corresponding to model gene families.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        A list of projected system information for the organism.
    """
    association = [] if not association else association
    projection = []

    sub_id = 1  # to keep track of genes that are part of the same component
    for cc in components:
        sys_state_in_org = "strict" if model_genes <= set(cc) else "split"
        for gene in cc:
            if gene.family in unit.model_families:
                metasource, metaid = unit.get_metainfo(gene.family)
                metadata = gene.family.get_metadata(metasource, metaid)
                avail_name = {fam.name for fam in unit.functional_unit.families}
                for fam in unit.functional_unit.families:
                    avail_name |= fam.exchangeable
                fam_annot = metadata.protein_name
                fam_sec = (
                    [
                        name
                        for name in metadata.secondary_names.split(",")
                        if name in avail_name
                    ]
                    if "secondary_name" in metadata.fields
                    else ""
                )
                category = "model"
            else:
                # separate genes filtered locally at family level; could be alternatively excluded
                category = "context" if gene.family in unit.families else "filtered"
                fam_annot = fam_sec = ""

            completeness = round(
                len({g.family.name for g in model_genes}) / len({unit.model_families}),
                2,
            )  # proportion of unit families in the organism

            line_projection = [
                unit.name,
                sub_id,
                gene.organism.name,
                gene.family.name,
                gene.family.named_partition,
                fam_annot,
                fam_sec,
                gene.ID,
                gene.local_identifier,
                gene.contig.name,
                gene.start,
                gene.stop,
                gene.strand,
                gene.is_fragment,
                category,
                sys_state_in_org,
                completeness,
                gene.product,
            ]

            if "RGPs" in association:
                rgp = gene.RGP
                unit.add_region(rgp) if rgp else None
                line_projection.append(str(rgp) if rgp else "")
            if "spots" in association:
                spot = gene.spot
                unit.add_spot(spot) if spot else None
                line_projection.append(str(spot) if spot else "")

            projection.append(list(map(str, line_projection)))
        sub_id += 1

    return projection


# This function replaced `compute_genes_graph`
def compute_gene_components(
    model_genes: Set[Gene], window_size: int
) -> List[List[Gene]]:
    """
    Compute gene components within a specified window size in the contigs of an organism.

    Args:
        model_genes (Set[Gene]): Set of genes in one organism corresponding to model gene families.
        window_size (int): The size of the window to consider for grouping genes.

    Returns:
        List[List[Gene]]: A list of components, each containing genes that are within the specified window.
    """
    contig_to_genes_of_interest = defaultdict(set)
    components = []

    for gene in model_genes:
        contig = gene.contig
        contig_to_genes_of_interest[contig].add(gene)

    for contig, genes_of_interest in contig_to_genes_of_interest.items():
        genes_count = contig.number_of_genes
        genes_of_interest_positions = [g.position for g in genes_of_interest]
        contig_windows = extract_contig_window(
            genes_count,
            genes_of_interest_positions,
            window_size=window_size,
            is_circular=contig.is_circular,
        )
        for window in contig_windows:
            comp = contig.get_genes(
                begin=window[0], end=window[1] + 1, outrange_ok=True
            )
            components.append(comp)

    return components


def unit_projection(
    unit: SystemUnit,
    gf2fam: Dict[GeneFamily, set[Family]],
    fam_index: Dict[GeneFamily, int],
    association: List[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system unit onto all organisms in a pangenome.

    Args:
        unit (SystemUnit): The system unit to project.
        gf2fam (Dict[str, set[Family]]): Dictionary linking a pangenome gene family to a model family.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    association = [] if not association else association
    pangenome_projection, organisms_projection = [], []
    matrix = get_gfs_matrix_combination(set(unit.model_families), gf2fam)
    mdr_acc_gfs = {gf for gf in unit.model_families if gf.name in matrix.columns.values}

    for organism in unit.model_organisms:
        # Note that `unit.model_organisms` is the set of organisms which have unit GFs >= `min_total` requirement of the unit, but not necessarily satisfying other unit requirements
        org_fam = {
            fam for fam in unit.families if organism.bitarray[fam_index[fam]] == 1
        }
        org_mod_fam = org_fam & mdr_acc_gfs

        filtered_matrix = matrix[[gf.name for gf in org_mod_fam]]

        if check_needed_families(filtered_matrix, unit.functional_unit):
            pan_proj = [
                unit.name,
                organism.name,
                ",".join(sorted([x.name for x in org_mod_fam])),  # model families
                ",".join(
                    sorted([x.name for x in org_fam - org_mod_fam])
                ),  # context families
            ]
            model_genes = {
                gene
                for family in org_mod_fam
                for gene in family.genes
                if gene.organism == organism
            }
            components = compute_gene_components(
                model_genes, unit.functional_unit.window
            )
            org_proj = project_unit_on_organisms(
                components, unit, model_genes, association
            )
            partition = conciliate_partition(
                set(
                    line[4] for line in org_proj if line[-4] == "model"
                )  # line[4] -> named_partition; line[-4] -> category
            )
            strict_count = sum(
                1 for line in org_proj if line[-3] == "strict"
            )  # line[-3] -> sys_state_in_org
            split_count = sum(1 for line in org_proj if line[-3] == "split")
            extended_count = len(org_proj) - strict_count - split_count
            # TODO completeness model
            completeness = len(org_fam) / len(unit)
            pangenome_projection.append(
                pan_proj
                + [partition, completeness]
                + [strict_count, split_count, extended_count]
            )
            if "RGPs" in association:
                rgps = {rgp.name for rgp in unit.regions if rgp.organism == organism}
                if len(rgps) == 1:
                    pangenome_projection[-1].extend(rgps)
                elif len(rgps) > 1:
                    join_rgps = [",".join(rgps)]
                    pangenome_projection[-1].extend(join_rgps)
            if "spots" in association:
                spots = {str(spot) for spot in unit.spots if organism in spot.organisms}
                if len(spots) == 1:
                    pangenome_projection[-1].extend(spots)
                elif len(spots) > 1:
                    join_spots = [",".join(spots)]
                    pangenome_projection[-1].extend(join_spots)

            organisms_projection += org_proj

    return (
        pd.DataFrame(pangenome_projection).drop_duplicates(),
        pd.DataFrame(organisms_projection).drop_duplicates(),
    )


def system_projection(
    system: System,
    fam_index: Dict[GeneFamily, int],
    gene_family2family: Dict[GeneFamily, Set[Family]],
    association: List[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system onto all organisms in a pangenome.

    Args:
        system (System): The system to project.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        gene_family2family (Dict[GeneFamily, Set[Family]]): Dictionary linking a gene family to model families.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for systems: {system.name}")
    begin = time.time()
    association = [] if not association else association
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    gf2fam = {
        gf: fam for gf, fam in gene_family2family.items() if gf in set(system.families)
    }

    for unit in system.units:
        unit_pan_proj, unit_org_proj = unit_projection(
            unit, gf2fam, fam_index, association
        )
        pangenome_projection = pd.concat(
            [pangenome_projection, unit_pan_proj], ignore_index=True
        )
        organisms_projection = pd.concat(
            [organisms_projection, unit_org_proj], ignore_index=True
        )

    if len(system) == 1:
        pangenome_projection = pd.concat(
            [
                pd.DataFrame(
                    [
                        [system.ID] * pangenome_projection.shape[0],
                        [system.name] * pangenome_projection.shape[0],
                    ]
                ).T,
                pangenome_projection,
            ],
            axis=1,
            ignore_index=True,
        )
        organisms_projection = pd.concat(
            [
                pd.DataFrame(
                    [
                        [system.ID] * organisms_projection.shape[0],
                        [system.name] * organisms_projection.shape[0],
                    ]
                ).T,
                organisms_projection,
            ],
            axis=1,
            ignore_index=True,
        )
    logging.getLogger("PANORAMA").debug(
        f"System projection done for {system.name} in {time.time() - begin} seconds"
    )
    return (
        pangenome_projection.drop_duplicates(),
        organisms_projection.drop_duplicates(),
    )


def project_pangenome_systems(
    pangenome: Pangenome,
    system_source: str,
    association: List[str] = None,
    canonical: bool = False,
    threads: int = 1,
    lock: Lock = None,
    disable_bar: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project systems onto all organisms in a pangenome.

    Args:
        pangenome (Pangenome): The pangenome to project.
        system_source (str): Source of the systems to project.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').
        canonical (bool, optional): If True, write the canonical version of systems too. Defaults to False.
        threads (int, optional): Number of threads available (default is 1).
        lock (Lock, optional): Global lock for multiprocessing execution (default is None).
        disable_bar (bool, optional): Disable progress bar (default is False).

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projections for each organism and the pangenome.
    """
    association = [] if not association else association
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()

    fam_index = pangenome.compute_org_bitarrays()
    meta2fam = get_metadata_to_families(
        pangenome, pangenome.systems_sources_to_metadata_source()[system_source]
    )
    sys2fam_context = {}

    for system in pangenome.systems:
        # Search association now to don't repeat for the same model and different system
        if system.model.name not in sys2fam_context:
            gf2fam, _ = dict_families_context(system.model, meta2fam)
            sys2fam_context[system.model.name] = gf2fam
        if canonical:
            for canonic in system.canonical:
                if canonic.model.name not in sys2fam_context:
                    gf2fam, _ = dict_families_context(canonic.model, meta2fam)
                    sys2fam_context[canonic.model.name] = gf2fam

    with ThreadPoolExecutor(
        max_workers=threads, initializer=init_lock, initargs=(lock,)
    ) as executor:
        logging.getLogger("PANORAMA").info(
            f"Begin system projection for source : {system_source}"
        )
        with tqdm(
            total=pangenome.number_of_systems(system_source, with_canonical=canonical),
            unit="system",
            disable=disable_bar,
        ) as progress:
            futures = []
            for system in pangenome.get_system_by_source(system_source):
                gf2fam = sys2fam_context[system.model.name]
                future = executor.submit(
                    system_projection,
                    system,
                    fam_index,
                    gf2fam,
                    association,
                )
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
                if canonical:
                    for canonic in system.canonical:
                        gf2fam = sys2fam_context[canonic.model.name]
                        future = executor.submit(
                            system_projection,
                            system,
                            fam_index,
                            gf2fam,
                            association,
                        )
                        future.add_done_callback(lambda p: progress.update())
                        futures.append(future)

            for future in futures:
                result = future.result()
                pangenome_projection = pd.concat(
                    [pangenome_projection, result[0]], ignore_index=True
                )
                organisms_projection = pd.concat(
                    [organisms_projection, result[1]], ignore_index=True
                )

    pan_cols_name = [
        "system number",
        "system name",
        "functional unit name",
        "organism",
        "model_GF",
        "context_GF",
        "partition",
        "completeness",
        "strict",
        "split",
        "extended",
    ]
    org_cols_name = [
        "system number",
        "system name",
        "functional unit name",
        "subsystem number",
        "organism",
        "gene family",
        "partition",
        "annotation",
        "secondary_names",
        "gene.ID",
        "gene.name",
        "contig",
        "start",
        "stop",
        "strand",
        "is_fragment",
        "category",
        "genomic organization",
        "completeness",
        "product",
    ]
    if "RGPs" in association:
        pan_cols_name += ["RGPs"]
        org_cols_name += ["RGPs"]
    if "spots" in association:
        pan_cols_name += ["spots"]
        org_cols_name += ["spots"]
    pangenome_projection.columns = pan_cols_name
    pangenome_projection.sort_values(
        by=[
            "system number",
            "system name",
            "functional unit name",
            "organism",
            "completeness",
        ],
        key=lambda col: (
            col.apply(extract_numeric_for_sorting)
            if col.name == "system number"
            else col
        ),
        ascending=[True, True, True, True, True],
        inplace=True,
    )
    organisms_projection.columns = org_cols_name
    organisms_projection.sort_values(
        by=[
            "system name",
            "system number",
            "subsystem number",
            "functional unit name",
            "organism",
            "contig",
            "start",
            "stop",
        ],
        key=lambda col: (
            col.apply(extract_numeric_for_sorting)
            if col.name in ["system number", "subsystem number"]
            else col
        ),
        ascending=[True, True, False, True, True, True, True, True],
        inplace=True,
    )
    logging.getLogger("PANORAMA").debug("System projection done")
    return pangenome_projection, organisms_projection


def get_org_df(org_df: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    """
    Get the reformated projection dataframe for an organism

    Args:
        org_df: Dataframe for the corresponding organism

    Returns:
        pd.DataFrame: Dataframe reformated for an organism
    Todo: This function is not used anymore, should we remove it?
    """
    org_name = str(org_df["organism"].unique()[0])
    org_df = org_df.drop(columns=["organism"])
    org_df_cols = org_df.columns.tolist()

    org_df_grouped = org_df.groupby(
        ["gene family", "system name", "functional unit name", "gene.ID", "start"],
        as_index=False,
    )
    agg_dict = {
        "system number": custom_agg_unique,
        "subsystem number": custom_agg,
        "partition": custom_agg_unique,
        "annotation": custom_agg,
        "secondary_names": custom_agg,
        "contig": custom_agg_unique,
        "gene.name": custom_agg_unique,
        "stop": custom_agg_unique,
        "strand": custom_agg_unique,
        "is_fragment": custom_agg_unique,
        "category": custom_agg,
        "genomic organization": custom_agg,
        "product": "first",
        "completeness": custom_agg,
    }
    if "RGPs" in org_df_cols:
        agg_dict["RGPs"] = custom_agg_unique
    if "spots" in org_df_cols:
        agg_dict["spots"] = custom_agg_unique
    org_df_grouped = org_df_grouped.agg(agg_dict)
    org_df_grouped = org_df_grouped[org_df_cols]

    # sort columns considering "system number" numerically
    org_df_grouped_sorted = org_df_grouped.sort_values(
        by=["system number", "system name", "start", "stop"],
        key=lambda col: (
            col.apply(extract_numeric_for_sorting)
            if col.name == "system number"
            else col
        ),
        ascending=[True, True, True, True],
    )
    return org_df_grouped_sorted, org_name


# This function is an alternative for `get_org_df` to only keep one unit per family
# it is also much faster since it avoids aggregation across all columns
def get_org_df_one_unit_per_fam(
    org_df, eliminate_filtered_systems=False, eliminate_empty_systems=False
) -> Tuple[pd.DataFrame, str]:
    """
    Filters and processes a DataFrame to retain only one representative unit
    per gene family based on completeness, and optionally eliminates certain
    systems based on filtering criteria. Also calculates overlapping information.

    Args:
        org_df: The input DataFrame containing organism data with details such as
            "gene family", "system name", "functional unit name", "completeness",
            and other related columns.
        eliminate_filtered_systems: Flag indicating whether to remove systems
            where any of their model families were filtered out due to lower
            completeness.
        eliminate_empty_systems: Flag indicating whether to remove systems with
            no model families left after filtering.

    Returns:
        Tuple containing:
            - A processed and filtered DataFrame with one row per unit per gene
              family, with overlapping information added and optional system
              elimination applied.
            - The unique organism name derived from the input DataFrame.
    """
    org_name = str(org_df["organism"].unique()[0])
    org_df = org_df.drop(columns=["organism"])
    # Keep only the row with the highest completeness for each group
    group_cols = [
        "gene family",
        "system name",
        "functional unit name",
        "gene.ID",
        "start",
    ]

    # Create an overlapping_units column before filtering
    overlapping_data = []
    for group_key, group in org_df.groupby(group_cols):
        if len(group) > 1:
            # Sort by completeness (descending) to keep the highest first
            group_sorted = group.sort_values("completeness", ascending=False)

            # Get the filtered out rows (all except the first)
            filtered_rows = group_sorted.iloc[1:]

            if not filtered_rows.empty:
                # Create overlapping info with the format `unit_number:completeness`
                overlapping_info = []
                for _, row in filtered_rows.iterrows():
                    overlapping_info.append(
                        f"{row['system number']}:{row['completeness']}"
                    )
                overlapping_str = "|".join(overlapping_info)
            else:
                overlapping_str = ""
        else:
            overlapping_str = ""

        overlapping_data.append(
            {
                "gene family": group_key[0],
                "system name": group_key[1],
                "functional unit name": group_key[2],
                "gene.ID": group_key[3],
                "start": group_key[4],
                "overlapping_units": overlapping_str,
            }
        )

    # Create overlapping DataFrame
    overlapping_df = pd.DataFrame(overlapping_data)

    idx_max_completeness = org_df.groupby(group_cols)["completeness"].idxmax()
    org_df_filtered = org_df.loc[idx_max_completeness]

    # Reset index to avoid issues with duplicate indices
    org_df_filtered = org_df_filtered.reset_index(drop=True)

    # Add the overlapping_units column to the end
    org_df_filtered = org_df_filtered.merge(overlapping_df, on=group_cols, how="left")

    # Fill NaN values in overlapping column with empty string
    org_df_filtered["overlapping_units"] = org_df_filtered["overlapping_units"].fillna(
        ""
    )

    if (
        eliminate_filtered_systems
    ):  # removes all systems with any of their model families filtered out due to lower completeness
        org_df_filtered = eliminate_systems(org_df, org_df_filtered)
    elif (
        eliminate_empty_systems
    ):  # To remove systems with no model families left after filtering
        org_df_filtered = eliminate_empty(org_df_filtered)

    # sort columns considering "system number" numerically
    org_df_grouped_sorted = org_df_filtered.sort_values(
        by=["system number", "system name", "start", "stop"],
        key=lambda col: (
            col.apply(extract_numeric_for_sorting)
            if col.name == "system number"
            else col
        ),
        ascending=[True, True, True, True],
    )
    return org_df_grouped_sorted, org_name


def eliminate_systems(org_df, org_df_filtered):
    """
    Eliminates systems from a filtered DataFrame based on model family changes. This function ensures that systems
    with any eliminated model families in the filtered dataset are removed entirely.

    Args:
        org_df: pandas.DataFrame containing the original dataset with all systems and associated gene families.
        org_df_filtered: pandas.DataFrame containing the already filtered dataset, which may have excluded some
            gene families or systems.

    Returns:
        pandas.DataFrame filtered to exclude entire systems where any model families were missing after the
        initial filtering step.
    """
    # Track which systems to eliminate
    systems_to_eliminate = set()

    # For each system, check if any model families were eliminated during filtering
    for system_number, system_group in org_df.groupby(
        ["system number", "system name", "functional unit name"]
    ):
        # Get original model families for this system
        original_model_families = set(
            system_group[system_group["category"] == "model"]["gene family"].unique()
        )

        # Get model families that survived filtering
        filtered_system_data = org_df_filtered[
            (org_df_filtered["system number"] == system_number[0])
            & (org_df_filtered["system name"] == system_number[1])
            & (org_df_filtered["functional unit name"] == system_number[2])
        ]
        surviving_model_families = set(
            filtered_system_data[filtered_system_data["category"] == "model"][
                "gene family"
            ].unique()
        )

        # If any model family was eliminated, mark the entire system for elimination
        if original_model_families != surviving_model_families:
            systems_to_eliminate.add(system_number)

    # Remove eliminated systems from the filtered dataframe
    for system_number in systems_to_eliminate:
        org_df_filtered = org_df_filtered[
            ~(
                (org_df_filtered["system number"] == system_number[0])
                & (org_df_filtered["system name"] == system_number[1])
                & (org_df_filtered["functional unit name"] == system_number[2])
            )
        ]
    return org_df_filtered


def eliminate_empty(org_df):
    """
    Removes systems with no model genes left.

    Args:
        org_df (pd.DataFrame): A DataFrame with at least the columns
            "system number" and "category". The "category" column is used to
            identify "model" genes, and the "system number" column groups rows
            into systems.

    Returns:
        pd.DataFrame: A DataFrame containing only the systems that have at least
            one "model" gene. The rows are concatenated and re-indexed.
    """
    # Removes systems with no model genes left
    valid_systems = []
    for _, system_group in org_df.groupby("system number"):
        model_genes_in_system = system_group[system_group["category"] == "model"]
        if not model_genes_in_system.empty:
            valid_systems.append(system_group)
    return pd.concat(valid_systems, ignore_index=True)


def write_projection_systems(
    output: Path,
    pangenome_projection: pd.DataFrame,
    organisms_projection: pd.DataFrame,
    organisms: List[str] = None,
    threads: int = 1,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Write the projected systems to output files.

    Args:
        output (Path): Path to the output directory.
        pangenome_projection (pd.DataFrame): DataFrame containing the pangenome projection.
        organisms_projection (pd.DataFrame): DataFrame containing the organism projections.
        organisms (List[str], optional): List of organisms to project (default is all organisms).
        threads (int, optional): Number of threads to use for parallel processing. Defaults to 1.
        force (bool, optional): Force write to the output directory (default is False).
        disable_bar (bool, optional): If True, disable the progress bar. Defaults to False.

    Returns:
        None
    """

    proj_dir = mkdir(output / "projection", force=force)
    if organisms is not None:
        pangenome_projection = pangenome_projection[
            ~pangenome_projection["organism"].isin(organisms)
        ]
        organisms_projection = organisms_projection[
            ~organisms_projection["organism"].isin(organisms)
        ]
    with ProcessPoolExecutor(
        max_workers=threads, mp_context=get_context("fork")
    ) as executor:
        futures = []
        for organism_name in pangenome_projection["organism"].unique():
            org_df = organisms_projection.loc[
                organisms_projection["organism"] == organism_name
            ]
            future = executor.submit(get_org_df_one_unit_per_fam, org_df)
            futures.append(future)

        for future in tqdm(
            as_completed(futures),
            total=len(pangenome_projection["organism"].unique()),
            unit="organisms",
            disable=disable_bar,
            desc="System projection on organisms",
        ):
            organism_df, organism_name = future.result()
            organism_df.to_csv(proj_dir / f"{organism_name}.tsv", sep="\t", index=False)

    pan_df_col = pangenome_projection.columns.tolist()
    pangenome_grouped = pangenome_projection.groupby(
        by=["system number", "system name"], as_index=False
    )
    agg_dict = {
        "functional unit name": custom_agg_unique,
        "organism": custom_agg_unique,
        "model_GF": custom_agg_unique,
        "context_GF": custom_agg_unique,
        "partition": get_partition,
        "completeness": "mean",
        "strict": "sum",
        "split": "sum",
        "extended": "sum",
    }
    if "RGPs" in pan_df_col:
        agg_dict["RGPs"] = custom_agg_unique
    if "spots" in pan_df_col:
        agg_dict["spots"] = custom_agg_unique

    pangenome_grouped = pangenome_grouped.agg(agg_dict)
    pangenome_grouped = pangenome_grouped[pan_df_col]

    # Create a temporary column for sorting based on the numeric values extracted
    pangenome_grouped["sort_key"] = pangenome_grouped["system number"].apply(
        extract_numeric_for_sorting
    )

    # Sort the DataFrame using the temporary column but keep the original values
    pangenome_sorted = pangenome_grouped.sort_values(
        by=["sort_key", "system name", "organism"], ascending=[True, True, True]
    ).drop(columns=["sort_key"])

    pangenome_sorted.to_csv(output / "systems.tsv", sep="\t", index=False)


# These functions could be moved to utils
# TODO Note that since _custom_agg uses sorting, the values at the same position in different columns will not correspond to same case
# TODO The product column sometimes has commas within one product name; it must be exchanged with another seperator to avoid issues
def _custom_agg(series: pd.Series, unique: bool = False):
    """
    Aggregate a column

    Args:
        series: series to aggregate
        unique: whether to return unique values or not

    Returns:
        The aggregated series
    """
    values = list(
        itertools.chain(*list(series.replace("", pd.NA).dropna().str.split(",")))
    )
    if not values:
        return ""

    if unique:
        values = set(values)

    return ", ".join(sorted(values, key=lambda x: int(x) if x.isdigit() else x))


def custom_agg(series: pd.Series):
    """
    Aggregate a column

    Args:
        series: series to aggregate

    Returns:
        The aggregated series
    """
    return _custom_agg(series, unique=False)


def custom_agg_unique(series: pd.Series):
    """
    Aggregate a column

    Args:
        series: series to aggregate

    Returns:
        The aggregated series
    """
    return _custom_agg(series, unique=True)


def get_partition(series: pd.Series):
    """

    Args:
        series:

    Returns:

    """

    final_partition = ""
    for partition in series.unique().tolist():
        if partition != final_partition:
            if final_partition == "":
                final_partition = partition
            else:
                if final_partition == "persistent":
                    if partition == "shell":
                        final_partition = "persistent/shell"
                    elif partition == "cloud":
                        final_partition = "persistent/cloud"
                    elif partition == "accessory":
                        final_partition = "persistent/accessory"
                elif final_partition == "shell":
                    if partition == "persistent":
                        final_partition = "persistent/shell"
                    elif partition == "cloud":
                        final_partition = "accessory"
                    elif partition == "persistent/cloud":
                        final_partition = "persistent/accessory"
                elif final_partition == "cloud":
                    if partition == "persistent":
                        final_partition = "persistent/cloud"
                    elif partition == "shell":
                        final_partition = "accessory"
                    elif partition == "persistent/shell":
                        final_partition = "persistent/accessory"
                elif final_partition == "accessory":
                    if partition == "persistent":
                        final_partition = "persistent/accessory"
                elif final_partition == "persistent/shell":
                    if partition == "cloud":
                        final_partition = "persistent/accessory"
                elif final_partition == "persistent/cloud":
                    if partition == "shell":
                        final_partition = "persistent/accessory"
    return final_partition


def extract_numeric_for_sorting(val) -> float:
    """
    Function to extract the numeric value for sorting while keeping the original value

    Args:
        val: the value

    Returns:
        float: the numeric value
    """
    try:
        # If it's a simple number, return it for sorting
        return float(val)
    except ValueError:
        # If it's a list of numbers separated by commas, return the smallest number for sorting
        if "," in val:
            parts = [float(x) for x in val.replace('"', "").split(",")]
            return (
                min(parts) + len(parts) * 0.0001
            )  # Sort by minimum, then by number of parts (by adding a small penalty for more parts)
        return float("inf")  # If it cannot be converted, place it at the end


# Replaced by `compute_gene_components`; preserved for further review if needed
def compute_genes_graph(model_genes: Set[Gene], unit: SystemUnit) -> nx.Graph:
    """
    Compute the genes graph for a given genomic context in an organism.

    Args:
        model_genes (Set[Gene]): Set of genes in one organism corresponding to model gene families.
        unit (SystemUnit): The unit of interest.

    Returns:
        nx.Graph: A genomic context graph for the given organism.
    """
    genes_graph = nx.Graph()
    for gene in sorted(model_genes, key=lambda x: x.position):
        if gene.position < gene.contig.number_of_genes:
            right_genes = gene.contig.get_genes(
                begin=gene.position,
                end=gene.position + unit.functional_unit.window + 1,
                outrange_ok=True,
            )
        else:
            right_genes = [gene]

        left_genes = gene.contig.get_genes(
            begin=gene.position - unit.functional_unit.window,
            end=gene.position + 1,
            outrange_ok=True,
        )
        for l_idx, l_gene in enumerate(left_genes, start=1):
            # if l_gene in genes_graph.nodes:
            for t, t_gene in enumerate(left_genes[l_idx:]):
                # if t_gene in genes_graph.nodes:
                if unit.functional_unit.same_strand:
                    # checking for functional unit same_strand at this point does not make sense
                    # as the unit is already detected at both pangenome and genome levels regardless of strand
                    if t_gene.strand == l_gene.strand:
                        genes_graph.add_edge(t_gene, l_gene, transitivity=t)
                else:
                    genes_graph.add_edge(t_gene, l_gene, transitivity=t)

        for r_idx, r_gene in enumerate(right_genes, start=1):
            # if r_gene in genes_graph.nodes:
            for t, t_gene in enumerate(right_genes[r_idx:]):
                # if t_gene in genes_graph.nodes:
                if unit.functional_unit.same_strand:
                    if t_gene.strand == r_gene.strand:
                        genes_graph.add_edge(t_gene, r_gene, transitivity=t)
                else:
                    genes_graph.add_edge(t_gene, r_gene, transitivity=t)
    return genes_graph
