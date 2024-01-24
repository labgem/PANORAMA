#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import time
import argparse
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from typing import Dict, Iterable, List, Set, Union
from multiprocessing import Manager, Lock

# installed libraries
import networkx as nx
from ppanggolin.genome import Gene
from ppanggolin.context.searchGeneContext import compute_gene_context_graph, compute_edge_metrics

# local libraries
from panorama.utils import init_lock
from panorama.utility.utility import check_models
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.format.read_binaries import load_pangenomes
from panorama.geneFamily import GeneFamily
from panorama.models import Models, Model, FuncUnit
from panorama.system import System
from panorama.pangenomes import Pangenome, Pangenomes


def check_pangenome_detection(pangenome: Pangenome, source: str, force: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: Source used to detect system
    :param force: Force to erase pangenome systems from source

    :raise KeyError: Provided source is not in the pangenome
    """
    if pangenome.status["systems"] == "inFile" and source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=source)
        else:
            raise Exception(f"Systems are already detected based on the source : {source}."
                            f" Use the --force option to erase the already computed systems.")
    if pangenome.status["metadata"]["families"] == "inFile":
        if source not in pangenome.status["metasources"]["families"]:
            raise KeyError(f"There is no metadata associate to familie "
                           f"from source {source} in pangenome {pangenome.name}.")
    elif pangenome.status["metadata"]["families"] not in ["Computed", "Loaded"]:
        raise Exception(f"There is no metadata associate to families in your pangenome {pangenome.name}. "
                        "Please see the command annotation before to detect systems")


def get_annotation_to_families(pangenome: Pangenome, source: str) -> Dict[str, Set[GeneFamily]]:
    """ Get for eache annotation a set of families with this annotation

    :param pangenome: Pangenome with gene families
    :param source: name of the annotation source

    :return: Dictionnary with for each annotation a set of gene families
    """
    annot2fam = {}
    for gf in pangenome.gene_families:
        metadata = gf.get_metadata_by_source(source)
        if metadata is not None:
            for meta in metadata:
                if meta.protein_name in annot2fam:
                    annot2fam[meta.protein_name].add(gf)
                else:
                    annot2fam[meta.protein_name] = {gf}
    return annot2fam


def dict_families_context(func_unit: FuncUnit, annot2fam: dict) -> (dict, dict):
    """Recover all families in the function unit

    :param func_unit: function unit object of model
    :param annot2fam: dictionary of annotated families

    :return families: dictionary of families interesting
    """
    families = dict()
    fam2annot = dict()
    for fam_model in func_unit.families:
        if fam_model.name in annot2fam:
            for gf in annot2fam[fam_model.name]:
                families[gf.name] = gf
                if gf.name in fam2annot:
                    fam2annot[gf.name].add(fam_model)
                else:
                    fam2annot[gf.name] = {fam_model}
        for exchangeable in fam_model.exchangeable:
            if exchangeable in annot2fam:
                for gf in annot2fam[exchangeable]:
                    families[gf.name] = gf
                    if gf.name in fam2annot:
                        fam2annot[gf.name].add(fam_model)
                    else:
                        fam2annot[gf.name] = {fam_model}
    return families, fam2annot


def compute_gc_graph(families: Iterable[GeneFamily], transitive: int = 4, window_size: int = 0,
                     jaccard: float = 0.85, disable_bar: bool = False) -> nx.Graph:
    """
    Compute the graph with transitive closure size and window provided as parameter and filter edges based on jaccard index

    Args:
        families:
        transitive:
        window_size:
        jaccard:
        disable_bar:

    Returns:
        nx.Graph: The gene context with edges filtered based on jaccard index

    Todo:
        This function is copy/paste from PPanGGOLiN make this as function in PPanGGOLiN to call it directly
    """
    start_time = time.time()

    logging.getLogger("PANORAMA").debug("Building the graph...")

    gene_context_graph = compute_gene_context_graph(families=families, transitive=transitive,
                                                    window_size=window_size, disable_bar=disable_bar)

    logging.getLogger("PANORAMA").debug(
        f"Took {round(time.time() - start_time, 2)} seconds to build the graph to find common gene contexts")

    logging.getLogger("PANORAMA").debug(
        f"Context graph made of {nx.number_of_nodes(gene_context_graph)} nodes and "
        f"{nx.number_of_edges(gene_context_graph)} edges")

    compute_edge_metrics(gene_context_graph, jaccard)

    filter_flag = f'is_jaccard_gene_>_{jaccard}'

    edges_to_remove = [(n, v) for n, v, d in gene_context_graph.edges(data=True) if not d[filter_flag]]
    gene_context_graph.remove_edges_from(edges_to_remove)

    return gene_context_graph


def search_fu_with_one_fam(func_unit: FuncUnit, annot2fam: dict, source: str, graph: nx.Graph) -> List[System]:
    """Search defense system with only one family necessary

    :param func_unit: Functional unit
    :param annot2fam: Dictionnary which link annotation with gene family
    :param source: name of the annotation source
    :param graph: connected component graph

    :return: Detected system in pangenome
    """
    detected_systems = []
    for mandatory_fam in func_unit.mandatory:
        if mandatory_fam.name in annot2fam:
            for pan_fam in annot2fam[mandatory_fam.name]:
                if pan_fam not in graph.nodes:
                    detected_systems.append(
                        System(system_id=0, model=func_unit.model, source=source, gene_families={pan_fam}))
        for exchangeable in mandatory_fam.exchangeable:
            if exchangeable in annot2fam:
                for pan_fam in annot2fam[exchangeable]:
                    if pan_fam not in graph.nodes:
                        detected_systems.append(System(system_id=0, model=func_unit.model,
                                                       source=source, gene_families={pan_fam}))
    return detected_systems


def extract_cc(node: GeneFamily, graph: nx.Graph, seen: set) -> Set[Union[Gene, GeneFamily]]:
    """ Get connected component of the gene family

    :param node: node corresponding to gene family
    :param graph: graph of connected component
    :param seen: set of already check gene families

    :return: set of connected component
    """
    nextlevel = {node}
    cc = set()
    while len(nextlevel) > 0:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in seen:
                cc.add(v)
                seen.add(v)
                nextlevel |= set(graph.neighbors(v))
    return cc


def check_cc(cc: set, families: dict, fam2annot: dict, func_unit: FuncUnit) -> bool:
    """Check parameters

    :param cc: set of connected component
    :param families: families find in context
    :param fam2annot: link between gene family and annotation
    :param func_unit:functional unit

    :return: Boolean true if parameter respected
    """
    count_forbidden, count_mandatory, count_total = (0, 0, 0)
    forbidden_list, mandatory_list, accessory_list = (list(map(lambda x: x.name, func_unit.forbidden)),
                                                      list(map(lambda x: x.name, func_unit.mandatory)),
                                                      list(map(lambda x: x.name, func_unit.accessory)))
    for fam, annot_set in fam2annot.items():
        for annot in annot_set:
            if annot.max_separation == -1:
                cc.add(families[fam])
    for node in sorted(cc, key=lambda n: len(fam2annot.get(n.name)) if fam2annot.get(n.name) is not None else 0):
        annotations = fam2annot.get(node.name)
        if annotations is not None:
            for annot in annotations:
                if annot.presence == 'forbidden' and annot.name in forbidden_list:  # if node is forbidden
                    count_forbidden += 1
                    forbidden_list.remove(annot.name)
                    if count_forbidden > func_unit.max_forbidden:
                        return False
                    break
                elif annot.presence == 'mandatory' and annot.name in mandatory_list:  # if node is mandatory
                    count_mandatory += 1
                    count_total += 1
                    mandatory_list.remove(annot.name)
                    break
                elif annot.presence == 'accessory' and annot.name in accessory_list:  # if node is accessory
                    count_total += 1
                    accessory_list.remove(annot.name)
                    break
    if (count_mandatory >= func_unit.min_mandatory or func_unit.min_mandatory == -1) and \
            (count_total >= func_unit.min_total or func_unit.min_total == -1):
        if (func_unit.max_mandatory >= count_mandatory or func_unit.max_mandatory == -1) and \
                (func_unit.max_total >= count_total or func_unit.max_total == -1):
            return True
        else:
            return False
    else:
        return False


def verify_param(g: nx.Graph(), families: dict, fam2annot: dict, model: Model, func_unit: FuncUnit, source: str,
                 detected_systems: list):
    """
    Check if the models parameters are respected

    Args:
        g: Connected component graph
        families: families find in context
        fam2annot: dictionary of families interesting
        model: Defined model
        func_unit: Functional unit
        source: annotation source
        detected_systems: detected system list
    """
    seen = set()

    remove_node = set()
    for node in g.nodes():
        if node not in seen:
            cc = extract_cc(node, g, seen)  # extract connect component
            if check_cc(cc, families, fam2annot, func_unit):
                detected_systems.append(System(system_id=0, model=model, gene_families=cc, source=source))
            else:
                remove_node |= cc
    g.remove_nodes_from(remove_node)


def search_system(model: Model, annot2fam: dict, source: str) -> List[System]:
    """
    Search if model system is in pangenome

    Args:
        model: Model to search in pangenome
        annot2fam: Dictionary with for each annotation a set of gene families
        source: Name of the annotation source

    Returns:
        List[System]: List of systems detected in pangenome for the given model
    """
    for func_unit in model.func_units:
        detected_systems = []
        families, fam2annot = dict_families_context(func_unit, annot2fam)
        g = compute_gc_graph(families=families.values(),
                             transitive=func_unit.max_separation + 1,
                             window_size=func_unit.max_separation + 2,
                             jaccard=0.25,
                             disable_bar=True)
        if func_unit.min_total == 1:
            detected_systems += search_fu_with_one_fam(func_unit, annot2fam, source, g)
        verify_param(g, families, fam2annot, model, func_unit, source, detected_systems)
        return detected_systems


def search_systems(models: Models, pangenome: Pangenome, source: str, threads: int = 1, disable_bar: bool = False):
    """
    Search systems present in the pangenome for all models

    Args:
        models: Models to search in pangenomes
        pangenome: Pangenome with gene families
        source: name of the annotation source
        threads: number of available threads
        disable_bar: Flag to disable progress bar
    """
    annot2fam = get_annotation_to_families(pangenome=pangenome, source=source)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=models.size, unit='model', disable=disable_bar) as progress:
            futures = []
            for model in models:
                future = executor.submit(search_system, model, annot2fam, source)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            detected_systems = []
            for future in futures:
                result = future.result()
                detected_systems += result
    for system in detected_systems:
        pangenome.add_system(system)
    logging.getLogger("PANORAMA").info(f"System detection done in pangenome {pangenome.name}")
    if len(detected_systems) > 0:
        pangenome.status["systems"] = "Computed"
        logging.getLogger("PANORAMA").info(f"Write systems in pangenome {pangenome.name}")
        write_pangenome(pangenome, pangenome.file, source=source, disable_bar=disable_bar)
        logging.getLogger("PANORAMA").info(f"Systems written in pangenome {pangenome.name}")
    else:
        logging.getLogger("PANORAMA").info("No system detected")


def search_systems_in_pangenomes(models: Models, pangenomes: Pangenomes, source: str, threads: int = 1,
                                 task: int = 1, lock: Lock = None, disable_bar: bool = False):
    """
    Search systems in pangenomes by multithreading on pangenomes

    Args:
        models: Models to search in pangenomes
        pangenomes: Getter object with Pangenome
        source: name of the annotation source
        threads: number of available threads
        task: number of parallel workers
        lock: Global lock for multiprocessing execution
        disable_bar: Disable progress bar

    Todo:
        - Replace the ThreadPoolExecutor with a for loop and shift systems progress bar.
    """
    with ThreadPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []

            for pangenome in pangenomes:
                logging.getLogger("PANORAMA").debug(
                    f"Align gene families to HMM for {pangenome.name} with {threads // task} threads...")
                future = executor.submit(search_systems, models, pangenome, source, threads // task, disable_bar)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def launch(args):
    """
    Launch functions to detect systems in pangenomes

    Args:
        args: argument given in CLI
    """
    models = check_models(args.models, disable_bar=args.disable_prog_bar)
    manager = Manager()
    lock = manager.Lock()
    need_info = {"need_annotations": True, "need_families": True, "need_metadata": True,
                 "metatypes": ["families"], "sources": [args.source]}
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_detection, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, source=args.source, force=args.force)
    search_systems_in_pangenomes(models=models, pangenomes=pangenomes, source=args.source, threads=args.threads,
                                 task=args.task, lock=lock, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for detection command

    Returns:
        argparse.ArgumentParser: parser arguments for align command
    """
    parser = sub_parser.add_parser("detection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_detection(parser)
    return parser


def parser_detection(parser):
    """
    Add argument to parser for detection command

    Args:
        parser: parser for detection argument

    TODO:
        - add an option to write projection
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument('-m', '--models', required=True, type=Path, nargs='?',
                          help="Path to model list file."
                               "Note: Use panorama utils --models to create the models list file")
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source where panorama as to select in pangenomes')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    optional.add_argument("--task", required=False, nargs='?', type=int, default=1,
                          help="Number of simultaneous task.")
