#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import logging
from typing import Dict, List, Set, Tuple, Union
from multiprocessing import Lock
from pathlib import Path

# installed libraries
from tqdm import tqdm
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism, Gene
from ppanggolin.context.searchGeneContext import compute_gene_context_graph

# local libraries
from panorama.utils import mkdir, init_lock
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.systems.system import System
from panorama.systems.models import Family
from panorama.systems.detection import get_annotation_to_families, dict_families_context, check_for_needed


def project_system_on_organisms(graph: nx.Graph, system: System, organism: Organism,
                                fam2annot: Dict[GeneFamily, Set[Family]]) -> Tuple[List[List[str]], List[int], str]:
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
        if gene.name != "":
            gene_name = gene.name
        else:
            gene_name = gene.local_identifier if gene.local_identifier != "" else gene.ID
        line_projection = [gene.family.name, gene.family.named_partition, annot, gene_name,
                           gene.start, gene.stop, gene.strand, gene.is_fragment, sys_state_in_org]
        return list(map(str, [system.ID, system.name, organism.name] + line_projection))

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
    model_family = {node for node in graph.nodes if
                    system.source in node.family.sources}  # Get all gene that have an annotation to code the system
    for cc in nx.connected_components(graph):
        model_cc = {gene for gene in cc if system.source in gene.family.sources}
        if model_cc.issuperset(model_family):  # Contain all model families in the system
            if len(cc) == len(model_family):  # Subgraph with only model families
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
    return projection, counter, conciliate_system_partition(partitions)


def compute_genes_graph(graph: nx.Graph, organism: Organism, t: int = 0) -> nx.Graph:
    """Compute the genes_graph for a given genomic context in an organism

    Args:
        graph: the genomic context graph
        organism: the organism of interest
        t: the transitive value

    Returns:
        A genomic context graph for the given organism
    """

    def compute_for_multigenic(family: GeneFamily, node: Gene):
        """
        Compute the edges in case of multigenic families

        Args:
            family: The multigenic families
            node: the gene of the other family to link with
        """
        for g in family.get_genes_per_org(organism):
            left_genes = g.contig.get_genes(begin=g.position,
                                            end=g.position + t) if g.position < g.contig.number_of_genes else [g]
            right_genes = g.contig.get_genes(begin=g.position - t, end=g.position)
            if node in left_genes or node in right_genes:
                genes_graph.add_edge(node, g)

    genes_graph = nx.Graph()
    for edge in graph.edges:
        if edge[0].is_multigenic_in_org(organism):
            if edge[1].is_multigenic_in_org(organism):
                for gene in edge[0].get_genes_per_org(organism):
                    compute_for_multigenic(edge[1], gene)
            else:
                v = list(edge[1].get_genes_per_org(organism))[0]
                compute_for_multigenic(edge[0], v)
        else:
            u = list(edge[0].get_genes_per_org(organism))[0]
            if edge[1].is_multigenic_in_org(organism):
                compute_for_multigenic(edge[1], u)
            else:
                v = list(edge[1].get_genes_per_org(organism))[0]
                genes_graph.add_edge(u, v)
    return genes_graph


def system_projection(system: System, annot2fam: Dict[str, Set[GeneFamily]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Project system on all pangenome organisms
    Args:
        system: system to project
        annot2fam: Dictionary with for each annotation a set of gene families

    Returns:
        2 Dataframe with projected system, one for the pangenome and another for the organisms
    """
    pangenome_projection, organisms_projection = [], []
    system_orgs = set.union(*list([set(gf.organisms) for gf in system.families]))
    func_unit = list(system.model.func_units)[0]
    t = func_unit.max_separation + 1
    graph = compute_gene_context_graph(families=set(system.families), transitive=t,
                                       window_size=t + 1, disable_bar=True)
    _, fam2annot = dict_families_context(func_unit, annot2fam)
    for organism in system_orgs:
        org_graph = graph.copy()
        org_graph.remove_nodes_from([n for n in graph.nodes if organism not in n.organisms])
        edges_to_remove = [(u, v) for u, v, e in org_graph.edges(data=True) if organism not in e['genomes']]
        org_graph.remove_edges_from(edges_to_remove)
        if check_for_needed(set(org_graph.nodes), fam2annot, func_unit):
            pan_proj = [system.ID, system.name, organism.name]
            genes_graph = compute_genes_graph(org_graph, organism, func_unit.max_separation + 1)
            org_proj, counter, partition = project_system_on_organisms(genes_graph, system, organism, fam2annot)
            pangenome_projection.append(pan_proj + [partition, len(org_graph.nodes) / len(system)] + counter)
            organisms_projection += org_proj
    return pd.DataFrame(pangenome_projection).drop_duplicates(), pd.DataFrame(organisms_projection).drop_duplicates()


def project_pangenome_systems(pangenome: Pangenome, source: str, threads: int = 1, lock: Lock = None,
                              disable_bar: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
         Write all systems in pangenomes and project on organisms

        Args:
            pangenome: Pangenome to project
            output: Path to output directory
            source: Annotation source
            organisms: List of organisms to project (defaults to all organisms)
            threads: Number of threads available (default: 1)
            lock: Global lock for multiprocessing execution (default: None)
            force: Force to write into the output directory (default: False)
            disable_bar: Allow to disable progress bar (default: False)

        Returns:
            Projection in 2 dataframe, one for each organism and one for the pangenome
        """
    pangenome_projection = pd.DataFrame()
    organisms_projection = pd.DataFrame()
    annot2fam = get_annotation_to_families(pangenome, source)
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        logging.getLogger("PANORAMA").info(f'Begin system projection for source : {source}')
        with tqdm(total=pangenome.number_of_systems(source, with_canonical=False), unit='system',
                  disable=disable_bar) as progress:
            futures = []
            for system in pangenome.get_system_by_source(source):
                future = executor.submit(system_projection, system, annot2fam)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                result = future.result()
                pangenome_projection = pd.concat([pangenome_projection, result[0]], ignore_index=True)
                organisms_projection = pd.concat([organisms_projection, result[1]], ignore_index=True)

    pangenome_projection.columns = ["system number", "system name", "organism", "partition",
                                    "completeness", "strict", "conserved", "split"]
    pangenome_projection.sort_values(by=["system number", "system name", "system number", "organism", "completeness"],
                                     ascending=[True, True, True, True, True], inplace=True)  # Try to order system number numericaly
    organisms_projection.columns = ["system number", "system name", "organism", "gene family", "partition",
                                    "annotation",
                                    "gene", "start", "stop", "strand", "is_fragment", "genomic organization"]
    organisms_projection.sort_values(by=["system name", "system number", "organism", "start", "stop"],
                                     ascending=[True, True, True, True, True], inplace=True)
    # pangenome_projection.insert(0, "pangenome name", pangenome_name)
    return pangenome_projection, organisms_projection


def write_projection_systems(pangenome_name: str, output: Path, source: str, pangenome_projection: pd.Dataframe,
                             organisms_projection:pd.DataFrame, organisms: List[str] = None, force: bool = False):
    """
     Write all systems in pangenomes and project on organisms

    Args:
        pangenome: Pangenome to project
        output: Path to output directory
        source: Annotation source
        organisms: List of organisms to project (defaults to all organisms)
        threads: Number of threads available (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        force: Force to write into the output directory (default: False)
        disable_bar: Allow to disable progress bar (default: False)

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
