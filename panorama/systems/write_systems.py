#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import argparse
import logging
from typing import Dict, List, Set, Tuple, Union
from multiprocessing import Manager, Lock
from pathlib import Path

# installed libraries
from tqdm import tqdm
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism, Gene
from ppanggolin.context.searchGeneContext import compute_gene_context_graph

# local libraries
from panorama.utils import mkdir, init_lock
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.systems.system import System
from panorama.systems.models import Family
from panorama.systems.detection import get_annotation_to_families, dict_families_context, check_for_needed


def check_write_systems_args(args: argparse.Namespace) -> None:
    """ Checks the provided arguments to ensure that they are valid.

    Args:
        args: The parsed arguments.

    Raises:
        argparse.ArgumentTypeError: If number of sources is not the same as models
        argparse.ArgumentTypeError: if annotation are given and their number is not the same as systems sources.
    """
    if len(args.sources) != len(args.models):
        raise argparse.ArgumentError(argument=None, message="Number of sources and models are different.")
    if args.annotation_sources is not None:
        if len(args.annotation_sources) != len(args.sources):
            raise argparse.ArgumentError(argument=None, message="Number of annotation sources is different from "
                                                                "number of systems sources.")
    else:
        args.annotation_sources = args.sources


def check_pangenome_write_systems(pangenome: Pangenome, sources: List[str]) -> None:
    """
     Check and load pangenome information before adding annotation

    Args:
        pangenome:  Pangenome object
        sources: Sources used to detect systems

    Raises:
        KeyError: Provided systems source is not in the pangenome
        Exception: Systems have not been detected in pangenome
        AttributeError: If there is no metadata associated to families
    """
    if pangenome.status["systems"] != "inFile":
        raise Exception("Systems have not been detected."
                        "Use 'panorama detect' subcommand to detect systems in pangenomes.")
    else:
        for systems_source in sources:
            if systems_source not in pangenome.status["systems_sources"]:
                logging.getLogger("PANORAMA").error(f"Systems in pangenome {pangenome.name} are: "
                                                    f"{pangenome.status['systems_sources']}")
                raise KeyError(f"There is no systems in pangenome {pangenome.name}, for the source: {systems_source}."
                               f"Look at 'panorama detect' subcommand to detect systems in {pangenome.name}.")


def project_system_on_organisms(graph: nx.Graph, system: System, organism: Organism,
                                fam2annot: Dict[GeneFamily, Set[Family]])\
        -> Tuple[List[List[Union[int, str]]], List[int]]:
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
    projection = []
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
                sys_state_in_org = "conserved"
        else:  # system is split
            counter[2] += 1
            sys_state_in_org = "split"
        for gene in cc:
            annot = list(fam2annot.get(gene.family.name))[0].name if fam2annot.get(
                gene.family.name) is not None else None
            line_projection = [gene.family.name, gene.family.named_partition, annot,
                               gene.name if gene.name != "" else gene.local_identifier if gene.local_identifier != "" else gene.ID,
                               gene.start, gene.stop, gene.strand, gene.is_fragment, sys_state_in_org]
            projection.append([system.ID, system.name, organism.name] + line_projection)
    return projection, counter


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


def write_system_projection(system: System, annot2fam: Dict[str, Set[GeneFamily]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
            pan_proj = [system.ID, system.name, organism.name, len(org_graph.nodes) / len(system)]
            genes_graph = compute_genes_graph(org_graph, organism, func_unit.max_separation + 1)
            org_proj, counter = project_system_on_organisms(genes_graph, system, organism, fam2annot)
            pangenome_projection.append(pan_proj + counter)
            organisms_projection += org_proj
    return pd.DataFrame(pangenome_projection).drop_duplicates(), pd.DataFrame(organisms_projection).drop_duplicates()


def write_projection_systems(pangenome: Pangenome, output: Path, source: str, organisms: List[str] = None,
                             threads: int = 1, lock: Lock = None, force: bool = False,
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
        logging.getLogger("PANORAMA").info(f'Write system projection for source : {source}')
        with tqdm(total=pangenome.number_of_systems(source, with_canonical=False), unit='system',
                  disable=disable_bar) as progress:
            futures = []
            for system in pangenome.get_system_by_source(source):
                future = executor.submit(write_system_projection, system, annot2fam)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                result = future.result()
                pangenome_projection = pd.concat([pangenome_projection, result[0]], ignore_index=True)
                organisms_projection = pd.concat([organisms_projection, result[1]], ignore_index=True)

    pangenome_projection.columns = ["system number", "system name", "organism",
                                    "completeness", "strict", "conserved", "split"]
    pangenome_projection.sort_values(by=["system name", "system number", "organism", "completeness"],
                                     ascending=[True, True, True, True], inplace=True)
    organisms_projection.columns = ["system number", "system name", "organism", "gene family", "partition",
                                    "annotation",
                                    "gene", "start", "stop", "strand", "is_fragment", "genomic organization"]
    organisms_projection.sort_values(by=["system name", "system number", "organism", "start", "stop"],
                                     ascending=[True, True, True, True, True], inplace=True)
    mkdir(output / f"{pangenome.name}/projection_{source}", force=force)
    if organisms is not None:
        pangenome_projection.drop(pangenome_projection["organism" in organisms].index)
        organisms_projection.drop(organisms_projection["organism" in organisms].index)
    for organism_name in pangenome_projection["organism"].unique():
        org_df = organisms_projection.loc[organisms_projection["organism"] == organism_name]
        org_df = org_df.drop(columns=["organism"])
        org_df.sort_values(by=["system number", "system name", "start", "stop"],
                           ascending=[True, True, True, True], inplace=True)
        org_df.to_csv(f"{output}/{pangenome.name}/projection_{source}/{organism_name}.tsv", sep="\t", index=False)
    pangenome_projection.to_csv(f"{output}/{pangenome.name}/systems_{source}.tsv", sep="\t", index=False)
    pangenome_projection.insert(0, "pangenome name", pangenome.name)
    return pangenome_projection, organisms_projection


def write_pangenomes_systems(pangenomes: Pangenomes, output: Path, projection: bool = False, association: str = None,
                             proksee: str = None, organisms: List[str] = None, threads: int = 1, lock: Lock = None,
                             force: bool = False, disable_bar: bool = False) -> None:
    """
    Write flat files about systems for all pangenomes

    Args:
        pangenomes: Pangenome objects with all pangenome
        output: Path to write flat files about systems
        projection: Flag to enable/disable pangenome projection (default: False)
        association: Flag to enable/disable pangenome association (default: False)
        proksee: Flag to enable/disable pangenome proksee generation (default: False)
        organisms: List of organism names to write (default: all organisms)
        threads: Number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        force: Flag to allow overwriting files (default: False)
        disable_bar: Flag to disable the progress bar (default: False)
    """
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Begin write systems for {pangenome.name}")
        for source in pangenome.systems_sources:
            logging.getLogger("PANORAMA").debug(f"Begin project systems for {pangenome.name} and source {source}.")
            if projection:
                projection = write_projection_systems(pangenome, output, source, organisms=organisms, threads=threads,
                                                      lock=lock, force=force, disable_bar=disable_bar)
            if proksee:
                raise NotImplementedError("Proksee not implemented")
            if association:
                raise NotImplementedError("Association not implemented")
            #     write_systems_projection()


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    from panorama.format.read_binaries import load_pangenomes
    from panorama.utility.utility import check_models

    check_write_systems_args(args)
    models_list = []
    for models in args.models:
        models_list.append(check_models(models, disable_bar=args.disable_prog_bar))
    outdir = mkdir(args.output, force=args.force)
    manager = Manager()
    lock = manager.Lock()
    need_info = {"need_annotations": True, "need_families": True, "need_graph": True,
                 "need_metadata": True, "metatypes": ["families"], "sources": args.annotation_sources,
                 "need_systems": True, "systems_sources": args.sources, "models": models_list}
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, check_function=check_pangenome_write_systems,
                                 need_info=need_info, sources=args.sources, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar)

    write_pangenomes_systems(pangenomes, outdir, projection=args.projection, proksee=args.proksee,
                             association=args.association, organisms=args.organisms, threads=args.threads,
                             lock=lock, force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write_systems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    required.add_argument('-m', '--models', required=True, type=Path, nargs="+",
                          help="Path to model list file. You can specify multiple models from different source. "
                               "For that separate the model list files by a space and "
                               "make sure you give them in the same order as the sources.")
    required.add_argument("-s", "--sources", required=True, type=str, nargs="+",
                          help="Name of the systems sources. You can specify multiple sources. "
                               "For that separate names by a space and "
                               "make sure you give them in the same order as the sources.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--projection", required=False, action="store_true",
                          help="Project the systems on organisms. If organisms are specified, "
                               "projection will be done only for them.")
    optional.add_argument("--association", required=False, type=str, default=None,
                          choices=["all", "rgp-modules", "rgp-spots", "modules-spots", "modules", "rgp"],
                          help="Write association between systems and others pangenomes elements")
    optional.add_argument("--proksee", required=False, type=str, default=None, nargs='+',
                          choices=["all", "base", "modules", "rgp", "spots", "annotations"],
                          help="Write a proksee file with systems. "
                               "If you want only the systems with genes, gene families and partition, use base value."
                               "Write rgps, spots or modules if you want them.")
    required.add_argument("--annotation_sources", required=False, type=str, nargs="+", default=None,
                          help="Name of the annotation sources if different from systems."
                               " You can specify multiple sources. For that separate names by a space and "
                               "make sure you give them in the same order as the sources.")
    optional.add_argument("--organisms", required=False, type=str, default=None, nargs='+')
    optional.add_argument("--threads", required=False, type=int, default=1)
