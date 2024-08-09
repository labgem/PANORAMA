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

# installed libraries
from tqdm import tqdm
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism, Gene

# local libraries
from panorama.utils import mkdir, init_lock
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.systems.system import System
from panorama.systems.models import Family
from panorama.systems.detection import get_metadata_to_families, dict_families_context, check_for_needed


def project_system_on_organisms(graph: nx.Graph, system: System, organism: Organism,
                                gene_fam2mod_fam: Dict[str, Set[Family]],
                                association: List[str] = None) -> Tuple[List[List[str]], List[int], str]:
    """
    Project a system onto an organism's pangenome.

    Args:
        graph (nx.Graph): Genomic context graph of the system for the given organism.
        system (System): The system to be projected.
        organism (Organism): The organism on which the system is to be projected.
        gene_fam2mod_fam (Dict[str, Set[Family]]): A dictionary mapping gene families to model families.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[List[List[str]], List[int], str]: A tuple containing:
            - A list of projected system information for the organism.
            - A list with counts of each system organization type (strict, extended, split).
            - The reconciled system partition.
    """
    def has_short_path(node_list, n):
        """
        Checks if there exists at least one path of length less than `n`
        connecting any two nodes in the given list of nodes in the graph.

        Args:
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

    def write_projection_line(gene: Gene) -> List[str]:
        """
        Write the projection information for a single gene.

        Args:
            gene (Gene): Gene to write projection for.

        Returns:
            List[str]: List of elements representing the projection for the gene.
        """
        line_projection = [gene.family.name, gene.family.named_partition, fam_annot, gene.ID, gene.local_identifier,
                           gene.start, gene.stop, gene.strand, gene.is_fragment, sys_state_in_org, gene.product]
        if 'RGPs' in association:
            rgp = gene.RGP
            if rgp is not None:
                system.add_region(rgp)
                line_projection.append(str(rgp))
            else:
                line_projection.append('')
        if 'spots' in association:
            spot = gene.spot
            if spot is not None:
                system.add_spot(gene.spot)
                line_projection.append(str(spot))
            else:
                line_projection.append('')
        return list(map(str, [system.ID, sub_id, system.name, organism.name] + line_projection))

    def conciliate_system_partition(system_partition: Set[str]) -> str:
        """
        Conciliate the partition of the system.

        Args:
            system_partition (Set[str]): All found partitions for genes coding the system.

        Returns:
            str: The reconciled system partition.

        Raises:
            Exception: If no partition is found. This may happen if partitions are not loaded or computed.
        """
        if len(system_partition) == 1:
            return system_partition.pop()
        else:
            if "persistent" in system_partition:
                return "persistent|accessory"
            else:
                return 'accessory'

    projection = []
    partitions = set()
    counter = [0, 0, 0]  # count strict, extended, and split CC

    model_genes = {gene for gene in graph.nodes if gene.family in system.models_families}
    sub_id = 1
    for cc in nx.connected_components(graph):
        model_cc = cc.intersection(model_genes)
        if len(model_cc) > 0:
            if model_cc == model_genes:
                func_unit = list(system.model.func_units)[0]
                if len(model_cc) == 1 or has_short_path(list(model_cc), func_unit.transitivity):
                    counter[0] += 1
                    sys_state_in_org = "strict"
                else:
                    counter[1] += 1
                    sys_state_in_org = "extended"
            else:
                counter[2] += 1
                sys_state_in_org = "split"
            for cc_gene in cc:
                if cc_gene.family.name in gene_fam2mod_fam:
                    metasource, metaid = system.get_metainfo(cc_gene.family)
                    if metaid != 0:
                        for mod_family in gene_fam2mod_fam[cc_gene.family.name]:
                            avail_name = {mod_family.name}.union(mod_family.exchangeable)
                            metadata = cc_gene.family.get_metadata(metasource, metaid)
                            if metadata.protein_name in avail_name:
                                fam_annot = metadata.protein_name
                                break
                            elif "secondary_name" in metadata.fields:
                                found = False
                                fam_annot = []
                                for name in metadata.secondary_name.split(","):
                                    if name in avail_name:
                                        found = True
                                        fam_annot.append(name)
                                fam_annot = ",".join(fam_annot)
                                if found:
                                    break
                        partitions.add(cc_gene.family.named_partition)
                    else:
                        fam_annot = ""
                else:
                    fam_annot = ""
                projection.append(write_projection_line(cc_gene))
            sub_id += 1
    return projection, counter, conciliate_system_partition(partitions)


def compute_genes_graph(families: Set[GeneFamily], organism: Organism, t: int = 0, w: int = 1) -> nx.Graph:
    """
    Compute the genes graph for a given genomic context in an organism.

    Args:
        families (Set[GeneFamily]): Set of gene families.
        organism (Organism): The organism of interest.
        t (int, optional): The transitive value (default is 0).
        w (int, optional): The window size for gene connection (default is 1).

    Returns:
        nx.Graph: A genomic context graph for the given organism.
    """
    genes_graph = nx.Graph()
    for family in families:
        genes_graph.add_nodes_from({gene for gene in family.genes if gene.organism == organism})
    for gene in genes_graph.nodes:
        if gene.position < gene.contig.number_of_genes:
            right_genes = gene.contig.get_genes(begin=gene.position, end=gene.position + w + 1, outrange_ok=True)
        else:
            right_genes = [gene]

        left_genes = gene.contig.get_genes(begin=gene.position - w, end=gene.position + 1, outrange_ok=True)
        for l_idx, l_gene in enumerate(left_genes, start=1):
            if l_gene in genes_graph.nodes:
                for t_gene in left_genes[l_idx:t + 1]:
                    if t_gene in genes_graph.nodes:
                        genes_graph.add_edge(t_gene, l_gene, transitivity=l_idx)
        for r_idx, r_gene in enumerate(right_genes, start=1):
            if r_gene in genes_graph.nodes:
                for t_gene in right_genes[r_idx:t + 1]:
                    if t_gene in genes_graph.nodes:
                        genes_graph.add_edge(t_gene, r_gene, transitivity=r_idx)
    return genes_graph


def system_projection(system: System, annot2fam: Dict[str, Dict[str, Set[GeneFamily]]],
                      fam_index: Dict[GeneFamily, int], association: List[str] = None
                      ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project a system onto all organisms in a pangenome.

    Args:
        system (System): The system to project.
        annot2fam (Dict[str, Dict[str, Set[GeneFamily]]]): Dictionary mapping annotations to gene families.
        fam_index (Dict[GeneFamily, int]): Index mapping gene families to their positions.
        association (List[str], optional): List of associations to include (e.g., 'RGPs', 'spots').

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two DataFrames containing the projected system for the pangenome and organisms.
    """
    pangenome_projection, organisms_projection = [], []
    func_unit = list(system.model.func_units)[0]
    gene_families, gf2fam, fam2source = dict_families_context(system.model, annot2fam)
    gene_families &= set(system.families)
    gene_families_name = {gf.name for gf in gene_families}
    gf2fam = {gf: fam for gf, fam in gf2fam.items() if gf in gene_families_name}

    for organism in system.models_organisms:
        org_fam = {fam for fam in system.families if organism.bitarray[fam_index[fam]] == 1}
        org_mod_fam = org_fam & set(system.models_families)
        check_needed, _ = check_for_needed(org_mod_fam, gf2fam, fam2source, func_unit)
        if check_needed:
            pan_proj = [system.ID, system.name, organism.name]
            genes_graph = compute_genes_graph(org_fam, organism, func_unit.transitivity, func_unit.window)
            org_proj, counter, partition = project_system_on_organisms(genes_graph, system, organism,
                                                                       gf2fam, association)
            pangenome_projection.append(pan_proj + [partition, len(org_fam) / len(system)] + counter)
            if 'RGPs' in association:
                rgps = {rgp.name for rgp in system.regions if rgp.organism == organism}
                if len(rgps) == 1:
                    pangenome_projection[-1].extend(rgps)
                elif len(rgps) > 1:
                    join_rgps = [','.join(rgps)]
                    pangenome_projection[-1].extend(join_rgps)
            if 'spots' in association:
                spots = {str(spot) for spot in system.spots if organism in spot.organisms}
                if len(spots) == 1:
                    pangenome_projection[-1].extend(spots)
                elif len(spots) > 1:
                    join_spots = [','.join(spots)]
                    pangenome_projection[-1].extend(join_spots)
            organisms_projection += org_proj
    logging.getLogger("PANORAMA").debug(f"System projection done for systems: {system.name}")
    return pd.DataFrame(pangenome_projection).drop_duplicates(), pd.DataFrame(organisms_projection).drop_duplicates()


def project_pangenome_systems(pangenome: Pangenome, system_source: str, association: List[str] = None, threads: int = 1,
                              lock: Lock = None, disable_bar: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project systems onto all organisms in a pangenome.

    Args:
        pangenome (Pangenome): The pangenome to project.
        system_source (str): Source of the systems to project.
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
    fam_index = pangenome.compute_org_bitarrays()
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        logging.getLogger("PANORAMA").info(f'Begin system projection for source : {system_source}')
        with tqdm(total=pangenome.number_of_systems(system_source, with_canonical=False), unit='system',
                  disable=disable_bar) as progress:
            futures = []
            for system in pangenome.get_system_by_source(system_source):
                future = executor.submit(system_projection, system, meta2fam, fam_index, association)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                result = future.result()
                pangenome_projection = pd.concat([pangenome_projection, result[0]], ignore_index=True)
                organisms_projection = pd.concat([organisms_projection, result[1]], ignore_index=True)
    pan_cols_name = ["system number", "system name", "organism", "partition",
                     "completeness", "strict", "extended", "split"]
    org_cols_name = ["system number", "subsystem number", "system name", "organism", "gene family",
                     "partition", "annotation", "gene.ID", "gene.name", "start", "stop", "strand",
                     "is_fragment", "genomic organization", "product"]
    if 'RGPs' in association:
        pan_cols_name += ['RGPs']
        org_cols_name += ['RGPs']
    if 'spots' in association:
        pan_cols_name += ['spots']
        org_cols_name += ['spots']
    pangenome_projection.columns = pan_cols_name
    pangenome_projection.sort_values(by=["system number", "system name", "organism", "completeness"],
                                     ascending=[True, True, True, True],
                                     inplace=True)  # TODO Try to order system number numerically
    organisms_projection.columns = org_cols_name
    organisms_projection.sort_values(by=["system name", "system number", "subsystem number",
                                         "organism", "start", "stop"],
                                     ascending=[True, True, False, True, True, True],
                                     inplace=True)
    logging.getLogger("PANORAMA").debug('System projection done')
    return pangenome_projection, organisms_projection


def write_projection_systems(output: Path, pangenome_projection: pd.DataFrame, organisms_projection: pd.DataFrame,
                             organisms: List[str] = None, force: bool = False):
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

    for organism_name in pangenome_projection["organism"].unique():
        org_df = organisms_projection.loc[organisms_projection["organism"] == organism_name]
        org_df = org_df.drop(columns=["organism"])
        org_df.sort_values(by=["system number", "system name", "start", "stop"],
                           ascending=[True, True, True, True], inplace=True)
        org_df.to_csv(proj_dir/f"{organism_name}.tsv", sep="\t", index=False)
    pangenome_projection.to_csv(output/'systems.tsv', sep="\t", index=False)
