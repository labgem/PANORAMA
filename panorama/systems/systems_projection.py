#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
import logging
from typing import Any, Dict, List, Set, Tuple
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
from panorama.systems.detection import get_annotation_to_families, dict_families_context, check_for_needed


def project_system_on_organisms(graph: nx.Graph, system: System, organism: Organism, fam2annot: Dict[str, Set[Family]],
                                association: List[str] = None) -> Tuple[List[List[str]], List[int], str]:
    """
    Project a system on a pangenome organism

    Args:
        graph: Genomic context graph of the system for the given organism
        system: The system to be projected
        organism: The organism on which the system is to be projected
        fam2annot: Dictionary mapping gene family to annotation

    Returns:
        The projected system on the organism and the number of each system organisation in projected organism
    """

    def write_projection_line(gene: Gene) -> List[str]:
        """
        Write a projection for one gene

        Args:
            gene: Gene to write projection

        Returns:
            List of element to write projection for a gene
        """
        line_projection = [gene.family.name, gene.family.named_partition, annot, gene.ID, gene.local_identifier,
                           gene.start, gene.stop, gene.strand, gene.is_fragment, sys_state_in_org, gene.product]
        if any(asso in association for asso in ['RGPs', 'spots']) and gene.RGP is not None:
            system.add_region(gene.RGP)
            if 'RGPs' in association:
                line_projection.append(gene.RGP.name)
            if 'spots' in association:
                line_projection.append(gene.spot)

        return list(map(str, [system.ID, sub_id, system.name, organism.name] + line_projection))

    def conciliate_system_partition(system_partition: Set[str]) -> str:
        """
        Conciliate the partition of the system

        Args:
            system_partition: All found partitions for gene that code the system

        Returns:
            Reconciled system partition

        Raises:
            Exception if not any partition are found. Could happen if partition are not loaded or computed
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
    counter = [0, 0, 0]  # count strict, conserved and split CC

    # Get all genes for which the family correspond to model
    model_genes = {gene for gene in graph.nodes if gene.family in system.models_families}
    sub_id = 1
    for cc in nx.connected_components(graph):
        model_cc = cc.intersection(model_genes)
        if len(model_cc) > 0:
            if model_cc == model_genes:  # Contain all model families in the system
                if len(cc) == len(model_genes):  # Subgraph with only model families
                    counter[0] += 1
                    sys_state_in_org = "strict"
                else:
                    counter[1] += 1
                    sys_state_in_org = "extended"
            else:  # system is split
                counter[2] += 1
                sys_state_in_org = "split"
            for cc_gene in cc:
                annot = list(fam2annot.get(cc_gene.family.name))[0].name if fam2annot.get(
                    cc_gene.family.name) is not None else None  # TODO look at the chosen annotation in case of multiple annotation per gene
                if annot is not None:
                    partitions.add(cc_gene.family.named_partition)
                projection.append(write_projection_line(cc_gene))
            sub_id += 1
    return projection, counter, conciliate_system_partition(partitions)


def compute_genes_graph(families: Set[GeneFamily], organism: Organism, t: int = 0, w: int = 1) -> nx.Graph:
    """Compute the genes_graph for a given genomic context in an organism

    Args:
        graph: the genomic context graph
        organism: the organism of interest
        t: the transitive value

    Returns:
        A genomic context graph for the given organism
    """

    genes_graph = nx.Graph()
    for family in families:
        # a = set({gene for gene in family.genes if gene.organism == organism})
        genes_graph.add_nodes_from({gene for gene in family.genes if gene.organism == organism})
    for gene in genes_graph.nodes:
        if gene.position < gene.contig.number_of_genes:
            right_genes = gene.contig.get_genes(begin=gene.position, end=gene.position + w + 1)
        else:
            right_genes = [gene]

        left_genes = gene.contig.get_genes(begin=gene.position - w, end=gene.position + 1)
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


def system_projection(system: System, annot2fam: Dict[str, Set[GeneFamily]], fam_index: Dict[GeneFamily, int],
                      association: List[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project system on all pangenome organisms
    Args:
        system: system to project
        annot2fam: Dictionary with for each annotation a set of gene families

    Returns:
        2 Dataframe with projected system, one for the pangenome and another for the organisms
    """
    pangenome_projection, organisms_projection = [], []
    func_unit = list(system.model.func_units)[0]
    t = func_unit.max_separation + 1
    _, fam2annot = dict_families_context(system.model, annot2fam)

    for organism in system.models_organisms:
        org_fam = {fam for fam in system.families if organism.bitarray[fam_index[fam]] == 1}
        if check_for_needed(org_fam, fam2annot, func_unit):
            pan_proj = [system.ID, system.name, organism.name]
            genes_graph = compute_genes_graph(org_fam, organism, t, t + 1)
            org_proj, counter, partition = project_system_on_organisms(genes_graph, system, organism,
                                                                       fam2annot, association)
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


def project_pangenome_systems(pangenome: Pangenome, system_source: str, annotation_sources: List[str],
                              association: List[str] = None, threads: int = 1, lock: Lock = None,
                              disable_bar: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
         Write all systems in pangenomes and project on organisms

        Args:
            pangenome: Pangenome to project
            system_source: source of the systems to project
            annotation_sources: sources of the annotation
            threads: Number of threads available (default: 1)
            lock: Global lock for multiprocessing execution (default: None)
            disable_bar: Allow to disable progress bar (default: False)

        Returns:
            Projection in 2 dataframe, one for each organism and one for the pangenome
        """
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    annot2fam = get_annotation_to_families(pangenome, annotation_sources)
    fam_index = pangenome.compute_org_bitarrays()
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        logging.getLogger("PANORAMA").info(f'Begin system projection for source : {system_source}')
        with tqdm(total=pangenome.number_of_systems(system_source, with_canonical=False), unit='system',
                  disable=disable_bar) as progress:
            futures = []
            for system in pangenome.get_system_by_source(system_source):
                system.families_sources = annotation_sources
                future = executor.submit(system_projection, system, annot2fam, fam_index, association)
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


def write_projection_systems(pangenome_name: str, output: Path, source: str, pangenome_projection: pd.Dataframe,
                             organisms_projection: pd.DataFrame, organisms: List[str] = None, force: bool = False):
    """
     Write all systems in pangenomes and project on organisms

    Args:
        pangenome_name: name of the current pangenome
        output: Path to output directory
        source: Annotation source
        organisms: List of organisms to project (defaults to all organisms)
        pangenome_projection: Dataframe of pangenome projection
        organisms_projection: Dataframe of organisms projection
        force: Force to write into the output directory (default: False)

    Returns:
        Projection in 2 dataframe, one for each organism and one for the pangenome
    """

    mkdir(output / f"{pangenome_name}/projection_{source}", force=force)
    if organisms is not None:
        pangenome_projection.drop(pangenome_projection["organism" in organisms].index)
        organisms_projection.drop(organisms_projection["organism" in organisms].index)

    for organism_name in pangenome_projection["organism"].unique():
        org_df = organisms_projection.loc[organisms_projection["organism"] == organism_name]
        org_df = org_df.drop(columns=["organism"])
        org_df.sort_values(by=["system number", "system name", "start", "stop"],
                           ascending=[True, True, True, True], inplace=True)
        org_df.to_csv(f"{output}/{pangenome_name}/projection_{source}/{organism_name}.tsv", sep="\t", index=False)
    pangenome_projection.to_csv(f"{output}/{pangenome_name}/systems_{source}.tsv", sep="\t", index=False)
