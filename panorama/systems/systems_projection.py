#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to project systems onto genomes.
"""

# default libraries
from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import logging
from typing import Dict, List, Set, Tuple
from multiprocessing import Lock
from pathlib import Path
import time

# installed libraries
from tqdm import tqdm
import numpy as np
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism, Gene

# local libraries
from panorama.utils import mkdir, init_lock, conciliate_partition
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.systems.utils import (get_metadata_to_families, dict_families_context,
                                    get_gfs_matrix_combination, check_needed_families)
from panorama.systems.system import System, SystemUnit
from panorama.systems.models import Family


def has_short_path(graph: nx.Graph, node_list, n):
    """
    Checks if there exists at least one path of length less than `n`
    connecting any two nodes in the given list of nodes in the graph.

    Args:
        graph (nx.Graph): the graph to search paths
        node_list (list): List of nodes to check for paths.
        n (int): The maximum length of the path to consider.

    Returns:
        bool: True if there exists at least one path of length less than `n`
              connecting any two nodes in the list, False otherwise.
    """
    path_length = defaultdict(dict)
    has_path = {node: False for node in node_list}
    for i, node1 in enumerate(node_list):
        for node2 in node_list[i + 1:]:
            if not has_path[node2]:
                try:
                    path_length[node1][node2] = nx.shortest_path_length(graph, source=node1, target=node2)
                    if path_length[node1][node2] <= n:
                        has_path[node1] = True
                        has_path[node2] = True
                        break
                except nx.NetworkXNoPath:
                    continue
    return all(has_path.values())


def project_unit_on_organisms(graph: nx.Graph, unit: SystemUnit, organism: Organism,
                              model_genes: Set[Gene], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                              association: List[str] = None) -> Tuple[List[List[str]], List[int], str]:
    """
    Project a system onto an organism's pangenome.

    Args:
        graph (nx.Graph): Genomic context graph of the system for the given organism.
        unit (SystemUnit): The unit to be projected.
        organism (Organism): The organism on which the system is to be projected.
        gene_fam2mod_fam (Dict[str, Set[Family]]): A dictionary mapping gene families to model families.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[List[List[str]], List[int], str]: A tuple containing:
            - A list of projected system information for the organism.
            - A list with counts of each system organization type (strict, extended, split).
            - The reconciled system partition.
    """

    def write_projection_line(gene: Gene) -> List[str]:
        """
        Write the projection information for a single gene.

        Args:
            gene (Gene): Gene to write projection for.

        Returns:
            List[str]: List of elements representing the projection for the gene.
        """
        line_projection = [gene.family.name, gene.family.named_partition, fam_annot, fam_sec, gene.ID,
                           gene.local_identifier, gene.start, gene.stop, gene.strand, gene.is_fragment,
                           sys_state_in_org, gene.product]

        if 'RGPs' in association:
            rgp = gene.RGP
            if rgp is not None:
                unit.add_region(rgp)
                line_projection.append(str(rgp))
            else:
                line_projection.append('')

        if 'spots' in association:
            spot = gene.spot
            if spot is not None:
                unit.add_spot(gene.spot)
                line_projection.append(str(spot))
            else:
                line_projection.append('')

        return list(map(str, [unit.name, sub_id, organism.name] + line_projection))

    projection = []
    partitions = set()
    counter = [0, 0, 0]  # count strict, extended, and split CC

    sub_id = 1
    for cc in nx.connected_components(graph):
        model_cc = cc.intersection(model_genes)
        if len(model_cc) > 0:
            if model_cc == model_genes:
                if len(model_cc) == 1 or has_short_path(graph, list(model_cc), unit.functional_unit.transitivity):
                    counter[0] += 1
                    sys_state_in_org = "strict"
                else:
                    counter[1] += 1
                    sys_state_in_org = "extended"
            else:
                counter[2] += 1
                sys_state_in_org = "split"
            for cc_gene in cc:
                fam_annot = ""
                fam_sec = []
                if cc_gene.family in gene_fam2mod_fam:
                    metasource, metaid = unit.get_metainfo(cc_gene.family)
                    if metaid != 0:
                        for mod_family in gene_fam2mod_fam[cc_gene.family]:
                            avail_name = {mod_family.name}.union(mod_family.exchangeable)
                            metadata = cc_gene.family.get_metadata(metasource, metaid)
                            if metadata.protein_name in avail_name:
                                fam_annot = metadata.protein_name
                            elif "secondary_name" in metadata.fields:
                                for name in metadata.secondary_name.split(","):
                                    if name in avail_name:
                                        fam_sec.append(name)
                        partitions.add(cc_gene.family.named_partition)
                fam_sec = ",".join(fam_sec)
                projection.append(write_projection_line(cc_gene))
            sub_id += 1
    return projection, counter, conciliate_partition(partitions)


def compute_genes_graph(families: Set[GeneFamily], organism: Organism,
                        unit: SystemUnit) -> Tuple[nx.Graph, Set[Gene]]:
    """
    Compute the genes graph for a given genomic context in an organism.

    Args:
        families (Set[GeneFamily]): Set of gene families.
        organism (Organism): The organism of interest.
        unit (SystemUnit): The unit of interest.

    Returns:
        nx.Graph: A genomic context graph for the given organism.
    """
    genes_graph = nx.Graph()
    for family in families:
        genes_graph.add_nodes_from({gene for gene in family.genes if gene.organism == organism})
    mod_fam = set(unit.models_families)
    model_genes = set()
    for gene in sorted(genes_graph.nodes, key=lambda x: x.position):
        if gene.family in mod_fam:
            model_genes.add(gene)
        if gene.position < gene.contig.number_of_genes:
            right_genes = gene.contig.get_genes(begin=gene.position,
                                                end=gene.position + unit.functional_unit.window + 1,
                                                outrange_ok=True)
        else:
            right_genes = [gene]

        left_genes = gene.contig.get_genes(begin=gene.position - unit.functional_unit.window, end=gene.position + 1,
                                           outrange_ok=True)
        for l_idx, l_gene in enumerate(left_genes, start=1):
            # if l_gene in genes_graph.nodes:
            for t_gene in left_genes[l_idx:unit.functional_unit.window]:
                # if t_gene in genes_graph.nodes:
                if unit.functional_unit.same_strand:
                    if t_gene.strand == l_gene.strand:
                        genes_graph.add_edge(t_gene, l_gene, transitivity=l_idx)
                else:
                    genes_graph.add_edge(t_gene, l_gene, transitivity=l_idx)

        for r_idx, r_gene in enumerate(right_genes, start=1):
            # if r_gene in genes_graph.nodes:
            for t_gene in right_genes[r_idx:unit.functional_unit.window]:
                # if t_gene in genes_graph.nodes:
                if unit.functional_unit.same_strand:
                    if t_gene.strand == r_gene.strand:
                        genes_graph.add_edge(t_gene, r_gene, transitivity=r_idx)
                else:
                    genes_graph.add_edge(t_gene, r_gene, transitivity=r_idx)
    return genes_graph, model_genes


def unit_projection(unit: SystemUnit, gf2fam: Dict[GeneFamily, set[Family]], fam2source: Dict[str, str],
                    fam_index: Dict[GeneFamily, int], association: List[str] = None
                    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system unit onto all organisms in a pangenome.

    Args:
        unit (SystemUnit): The system unit to project.
        gf2fam (Dict[str, set[Family]]): Dictionary linking a pangenome gene family to a model family.
        fam2source (Dict[str, str]): Dictionary linking a model family to his source.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    pangenome_projection, organisms_projection = [], []
    matrix, _, _ = get_gfs_matrix_combination(set(unit.models_families), gf2fam, fam2source)
    for organism in unit.models_organisms:
        org_fam = {fam for fam in unit.families if organism.bitarray[fam_index[fam]] == 1}
        org_mod_fam = org_fam & set(unit.models_families)
        filtered_matrix = matrix[list({gf.name for gf in org_mod_fam})]
        if check_needed_families(filtered_matrix, unit.functional_unit):
            pan_proj = [unit.name, organism.name, ", ".join(sorted([x.name for x in org_mod_fam])),
                        ", ".join(sorted([x.name for x in org_fam - org_mod_fam]))]
            genes_graph, model_genes = compute_genes_graph(org_mod_fam, organism, unit)
            org_proj, counter, partition = project_unit_on_organisms(genes_graph, unit, organism, model_genes,
                                                                     gf2fam, association)
            pangenome_projection.append(pan_proj + [partition, len(org_fam) / len(unit)] + counter)
            if 'RGPs' in association:
                rgps = {rgp.name for rgp in unit.regions if rgp.organism == organism}
                if len(rgps) == 1:
                    pangenome_projection[-1].extend(rgps)
                elif len(rgps) > 1:
                    join_rgps = [','.join(rgps)]
                    pangenome_projection[-1].extend(join_rgps)
            if 'spots' in association:
                spots = {str(spot) for spot in unit.spots if organism in spot.organisms}
                if len(spots) == 1:
                    pangenome_projection[-1].extend(spots)
                elif len(spots) > 1:
                    join_spots = [','.join(spots)]
                    pangenome_projection[-1].extend(join_spots)
            organisms_projection += org_proj
    return pd.DataFrame(pangenome_projection).drop_duplicates(), pd.DataFrame(organisms_projection).drop_duplicates()


def system_projection(system: System, fam_index: Dict[GeneFamily, int], gene_families: Set[GeneFamily],
                      gene_family2family: Dict[GeneFamily, Set[Family]], fam2source: Dict[str, str],
                      association: List[str] = None
                      ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system onto all organisms in a pangenome.

    Args:
        system (System): The system to project.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        gene_families (Set[GeneFamily]): The set of gene families that code for the model corresponding to system.
        gene_family2family (Dict[GeneFamily, Set[Family]]): Dictionary linking a gene family to model families.
        fam2source (Dict[str, str]): Dictionary linking a model family to his source.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for systems: {system.name}")
    begin = time.time()
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    gfs = gene_families & set(system.families)
    gf2fam = {gf: fam for gf, fam in gene_family2family.items() if gf in gfs}

    for unit in system.units:
        unit_pan_proj, unit_org_proj = unit_projection(unit, gf2fam, fam2source, fam_index, association)
        pangenome_projection = pd.concat([pangenome_projection, unit_pan_proj], ignore_index=True)
        organisms_projection = pd.concat([organisms_projection, unit_org_proj], ignore_index=True)

    if len(system) == 1:
        pangenome_projection = pd.concat([pd.DataFrame([[system.ID] * pangenome_projection.shape[0],
                                                        [system.name] * pangenome_projection.shape[0]]).T,
                                          pangenome_projection], axis=1, ignore_index=True)
        organisms_projection = pd.concat([pd.DataFrame([[system.ID] * organisms_projection.shape[0],
                                                        [system.name] * organisms_projection.shape[0]]).T,
                                          organisms_projection], axis=1, ignore_index=True)
    logging.getLogger("PANORAMA").debug(f"System projection done for {system.name} in {time.time() - begin} seconds")
    return pangenome_projection.drop_duplicates(), organisms_projection.drop_duplicates()


def project_pangenome_systems(pangenome: Pangenome, system_source: str, fam_index: Dict[GeneFamily, int],
                              association: List[str] = None, canonical: bool = False, threads: int = 1,
                              lock: Lock = None, disable_bar: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project systems onto all organisms in a pangenome.

    Args:
        pangenome (Pangenome): The pangenome to project.
        system_source (str): Source of the systems to project.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').
        threads (int, optional): Number of threads available (default is 1).
        lock (Lock, optional): Global lock for multiprocessing execution (default is None).
        disable_bar (bool, optional): Disable progress bar (default is False).

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projections for each organism and the pangenome.
    """
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    meta2fam = get_metadata_to_families(pangenome, pangenome.systems_sources_to_metadata_source()[system_source])
    sys2fam_context = {}

    for system in pangenome.systems:  # Search association now to don't repeat for same model and different system
        if system.model.name not in sys2fam_context:
            gene_families, gf2fam, fam2source = dict_families_context(system.model, meta2fam)
            sys2fam_context[system.model.name] = (gene_families, gf2fam, fam2source)
        if canonical:
            for canonic in system.canonical:
                if canonic.model.name not in sys2fam_context:
                    gene_families, gf2fam, fam2source = dict_families_context(canonic.model, meta2fam)
                    sys2fam_context[canonic.model.name] = (gene_families, gf2fam, fam2source)

    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        logging.getLogger("PANORAMA").info(f'Begin system projection for source : {system_source}')
        with tqdm(total=pangenome.number_of_systems(system_source, with_canonical=canonical), unit='system',
                  disable=disable_bar) as progress:
            futures = []
            for system in pangenome.get_system_by_source(system_source):
                gene_families, gf2fam, fam2source = sys2fam_context[system.model.name]
                future = executor.submit(system_projection, system, fam_index, gene_families,
                                         gf2fam, fam2source, association)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
                if canonical:
                    for canonic in system.canonical:
                        gene_families, gf2fam, fam2source = sys2fam_context[canonic.model.name]
                        future = executor.submit(system_projection, system, fam_index, gene_families,
                                                 gf2fam, fam2source, association)
                        future.add_done_callback(lambda p: progress.update())
                        futures.append(future)

            for future in futures:
                result = future.result()
                pangenome_projection = pd.concat([pangenome_projection, result[0]], ignore_index=True)
                organisms_projection = pd.concat([organisms_projection, result[1]], ignore_index=True)
    pan_cols_name = ["system number", "system name", "functional unit name", "organism", "model_GF", "context_GF",
                     "partition", "completeness", "strict", "extended", "split"]
    org_cols_name = ["system number", "system name", "functional unit name", "subsystem number", "organism",
                     "gene family", "partition", "annotation", "secondary_names", "gene.ID", "gene.name", "start",
                     "stop", "strand", "is_fragment", "genomic organization", "product"]
    if 'RGPs' in association:
        pan_cols_name += ['RGPs']
        org_cols_name += ['RGPs']
    if 'spots' in association:
        pan_cols_name += ['spots']
        org_cols_name += ['spots']
    pangenome_projection.columns = pan_cols_name
    pangenome_projection.sort_values(
        by=["system number", "system name", "functional unit name", "organism", "completeness"],
        ascending=[True, True, True, True, True],
        inplace=True)  # TODO Try to order system number numerically
    organisms_projection.columns = org_cols_name
    organisms_projection.sort_values(
        by=["system name", "system number", "subsystem number", "functional unit name", "organism", "start", "stop"],
        ascending=[True, True, False, True, True, True, True],
        inplace=True)
    logging.getLogger("PANORAMA").debug('System projection done')
    return pangenome_projection, organisms_projection


def custom_agg(series: pd.Series):
    """
    Aggregate a column

    Args:
        series: series to aggregate

    Returns:
        The aggregated series
    """
    # if all value are identical, first is chosen
    if series.nunique() == 1:
        return series.iloc[0]
    # else, all values are joined
    else:
        return ', '.join(series.dropna())


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
        if ',' in val:
            parts = [float(x) for x in val.replace('"', '').split(',')]
            return min(parts)  # Take the minimum for sorting
        return float('inf')  # If it cannot be converted, place it at the end


def get_org_df(org_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get the reformated projection dataframe for an organism

    Args:
        org_df: Dataframe for the corresponding organism

    Returns:
        pd.DataFrame: Dataframe reformated for an organism
    """
    org_df = org_df.drop(columns=["organism"])
    org_df_cols = org_df.columns.tolist()
    org_df_grouped = org_df.groupby(["gene family", "system name", "functional unit name"],
                                    as_index=False)
    org_df_grouped = org_df_grouped.agg({"system number": ','.join, "subsystem number": custom_agg,
                                         "partition": 'first', "annotation": custom_agg, "gene.ID": custom_agg,
                                         "gene.name": custom_agg, "secondary_names": custom_agg,
                                         "start": custom_agg,
                                         "stop": custom_agg, "strand": custom_agg, "is_fragment": custom_agg,
                                         "genomic organization": custom_agg, "product": custom_agg,
                                         "RGPs": custom_agg, "spots": custom_agg})
    org_df_grouped = org_df_grouped[org_df_cols]

    # Create a temporary column for sorting based on the numeric values extracted
    org_df_grouped["sort_key"] = org_df_grouped["system number"].apply(extract_numeric_for_sorting)

    # Sort the DataFrame using the temporary column, but keep the original values
    org_df_sorted = org_df_grouped.sort_values(by=["sort_key", "system name", "start", "stop"],
                                               ascending=[True, True, True, True]).drop(columns=["sort_key"])
    return org_df_sorted


def write_projection_systems(output: Path, pangenome_projection: pd.DataFrame, organisms_projection: pd.DataFrame,
                             organisms: List[str] = None, force: bool = False, disable_bar: bool = False):
    """
    Write the projected systems to output files.

    Args:
        output (Path): Path to the output directory.
        pangenome_projection (pd.DataFrame): DataFrame containing the pangenome projection.
        organisms_projection (pd.DataFrame): DataFrame containing the organism projections.
        organisms (List[str], optional): List of organisms to project (default is all organisms).
        force (bool, optional): Force write to the output directory (default is False).

    Returns:
        None
    """

    proj_dir = mkdir(output / "projection", force=force)
    if organisms is not None:
        pangenome_projection = pangenome_projection[~pangenome_projection["organism"].isin(organisms)]
        organisms_projection = organisms_projection[~organisms_projection["organism"].isin(organisms)]

    for organism_name in tqdm(pangenome_projection["organism"].unique(), unit='organisms', disable=disable_bar):
        org_df = organisms_projection.loc[organisms_projection["organism"] == organism_name]
        result = get_org_df(org_df)
        result.to_csv(proj_dir / f"{organism_name}.tsv", sep="\t", index=False)

    pangenome_grouped = pangenome_projection.groupby(by=["system number", "system name"], as_index=False)
    pangenome_grouped = pangenome_grouped.agg({"functional unit name": custom_agg, "organism": custom_agg,
                                               "model_GF": custom_agg, "context_GF": custom_agg,
                                               "partition": get_partition, "completeness": np.mean,
                                               "strict": sum, "extended": sum, "split": sum,
                                               "RGPs": custom_agg, "spots": custom_agg})
    pangenome_grouped = pangenome_grouped[pangenome_projection.columns.tolist()]

    # Create a temporary column for sorting based on the numeric values extracted
    pangenome_grouped["sort_key"] = pangenome_grouped["system number"].apply(extract_numeric_for_sorting)

    # Sort the DataFrame using the temporary column, but keep the original values
    pangenome_sorted = pangenome_grouped.sort_values(by=["sort_key", "system name", "organism"],
                                                     ascending=[True, True, True]).drop(columns=["sort_key"])

    pangenome_sorted.to_csv(output / 'systems.tsv', sep="\t", index=False)
