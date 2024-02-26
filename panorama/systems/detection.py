#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import argparse
from itertools import combinations
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from typing import Dict, List, Set, Tuple
from multiprocessing import Manager, Lock
from collections import defaultdict

# installed libraries
import networkx as nx
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import compute_gene_context_graph

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Models, Model, FuncUnit, Family
from panorama.systems.system import System
from panorama.pangenomes import Pangenome, Pangenomes


def check_detection_parameters(args: argparse.Namespace) -> None:
    """
    Checks the provided arguments to ensure that they are valid.

    Args:
        args: The parsed arguments.

    Raises:
        argparse.ArgumentTypeError: If jaccard is not with the correct type
    """
    args.jaccard = restricted_float(args.jaccard)
    args.annotation_source = args.source if args.annotation_source is None else args.annotation_source


def check_pangenome_detection(pangenome: Pangenome, annotation_source: str, systems_source: str,
                              force: bool = False) -> None:
    """
     Check and load pangenome information before adding annotation

    Args:
        pangenome:  Pangenome object
        annotation_source: Source used to annotate gene famillies
        systems_source: Source used to detect system
        force: Force to erase pangenome systems from source

    Raises:
        KeyError: Provided annotation source is not in the pangenome
        Exception: If System already exists in the Pangenome
        AttributeError: If there is no metadata associated to families
    """
    from panorama.format.write_binaries import erase_pangenome
    if pangenome.status["systems"] == "inFile" and systems_source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=systems_source)
        else:
            raise Exception(f"Systems are already detected based on the source : {systems_source}."
                            f" Use the --force option to erase the already computed systems.")
    if pangenome.status["metadata"]["families"] == "inFile":
        if annotation_source not in pangenome.status["metasources"]["families"]:
            raise KeyError(f"There is no metadata associate to families "
                           f"from source {annotation_source} in pangenome {pangenome.name}.")
    elif pangenome.status["metadata"]["families"] not in ["Computed", "Loaded"]:
        raise AttributeError(f"There is no metadata associate to families in your pangenome {pangenome.name}. "
                             "Please see the command annotation before to detect systems")


def get_annotation_to_families(pangenome: Pangenome, source: str) -> Dict[str, Set[GeneFamily]]:
    """
    Get for each annotation a set of families with this annotation

    Args:
        pangenome: Pangenome with gene families
        source: name of the annotation source

    Returns:
        Dictionary with for each annotation a set of gene families
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


def dict_families_context(func_unit: FuncUnit, annot2fam: Dict[str, Set[GeneFamily]]) \
        -> Tuple[Set[GeneFamily], Dict[GeneFamily, Set[Family]]]:
    """
    Recover all families in the function unit

    Args:
        func_unit: function unit object of model
        annot2fam: dictionary of annotated families

    Returns:
        The set of families of interest in the functional unit and
        a dictionary which link families to their annotation for the functional unit
    """
    families = set()
    fam2annot = defaultdict(set)
    for fam_model in func_unit.families:
        if fam_model.name in annot2fam:
            for gf in annot2fam[fam_model.name]:
                families.add(gf)
                fam2annot[gf.name].add(fam_model)
        for exchangeable in fam_model.exchangeable:
            if exchangeable in annot2fam:
                for gf in annot2fam[exchangeable]:
                    families.add(gf)
                    fam2annot[gf.name].add(fam_model)
    return families, fam2annot


def filter_local_context(graph: nx.Graph, families: Set[GeneFamily], jaccard_threshold: float = 0.8) -> None:
    """
    Filtered a graph based on a local jaccard index

    Args:
        graph: A sub pangenome graph
        families: list of families that code for the system
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
    """
    # Compute local jaccard
    graph_orgs = {gene.organism for u, v, data in graph.edges(data=True) for genes in data['genes'].values()
                  for gene in genes if u in set(families) and v in set(families)}
    for f1, f2, data in graph.edges(data=True):
        f1_genes_in_orgs = {gene for gene in f1.genes if gene.organism in graph_orgs}
        f2_genes_in_orgs = {gene for gene in f2.genes if gene.organism in graph_orgs}
        f1_gene_proportion = len(data['genes'][f1]) / len(f1_genes_in_orgs) if len(f1_genes_in_orgs) > 0 else 0
        f2_gene_proportion = len(data['genes'][f2]) / len(f2_genes_in_orgs) if len(f2_genes_in_orgs) > 0 else 0

        data['f1'] = f1.name
        data['f2'] = f2.name
        data['f1_jaccard_gene'] = f1_gene_proportion
        data['f2_jaccard_gene'] = f2_gene_proportion

        data[f'is_jaccard_gene_>_{jaccard_threshold}'] = ((f1_gene_proportion >= jaccard_threshold) and
                                                          (f2_gene_proportion >= jaccard_threshold))

    # Filter cc
    filter_flag = f'is_jaccard_gene_>_{jaccard_threshold}'

    edges_to_remove = [(n, v) for n, v, d in graph.edges(data=True) if not d[filter_flag]]
    graph.remove_edges_from(edges_to_remove)


def check_for_forbidden(families: Set[GeneFamily], fam2annot: Dict[GeneFamily, Set[Family]],
                        func_unit: FuncUnit) -> bool:
    """
    Check if there is forbidden in families.

    Args:
        families: set of families
        fam2annot: link between gene family and annotation
        func_unit: functional unit

    Returns:
        Boolean true if forbidden condition are encountered in families
    """
    count_forbidden = 0
    forbidden_list = list(map(lambda x: x.name, func_unit.forbidden))
    for node in sorted(families, key=lambda n: len(fam2annot.get(n.name)) if fam2annot.get(n.name) is not None else 0):
        annots = fam2annot.get(node.name)
        if annots is not None:
            for annot in annots:
                if annot.presence == 'forbidden' and annot.name in forbidden_list:  # if node is forbidden
                    count_forbidden += 1
                    forbidden_list.remove(annot.name)
                    if count_forbidden > func_unit.max_forbidden:
                        return True
                    break
    return False


def check_for_needed(families: Set[GeneFamily], fam2annot: Dict[GeneFamily, Set[Family]],
                     func_unit: FuncUnit) -> bool:
    """
    Check if presence/absence rules for needed families are respected

    Args:
        families: set of families
        fam2annot: link between gene family and annotation
        func_unit: functional unit

    Returns:
        Boolean true if parameter respected
    """
    count_mandatory, count_total = (0, 0)
    mandatory_list, accessory_list = (list(map(lambda x: x.name, func_unit.mandatory)),
                                      list(map(lambda x: x.name, func_unit.accessory)))
    for node in sorted(families, key=lambda n: len(fam2annot.get(n.name)) if fam2annot.get(n.name) is not None else 0):
        annots = fam2annot.get(node.name)
        if annots is not None:
            for annot in annots:
                if annot.presence == 'mandatory' and annot.name in mandatory_list:  # if node is mandatory
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


def search_system_in_graph(graph: nx.Graph, families: Set[GeneFamily], fam2annot: Dict[GeneFamily, Set[Family]],
                           func_unit: FuncUnit, source: str, jaccard_threshold: float = 0.8,
                           seen: Set[Tuple[GeneFamily]] = None, depth: int = 0, max_depth: int = 1) -> Set[System]:
    """Search systems corresponding to model in a graph

    Args:
        graph: A graph with families in one connected component
        families: A set of families that code for the searched model
        fam2annot: A dictionary to link gene families to their annotation
        func_unit: One functional unit corresponding to the model
        source: Name of the source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        seen: Set of combination of families already pass
        depth: Depth of the recursion search
        max_depth: Maximum depth of the recursion

    Returns:
        A set of all detected systems in the graph
    """
    detected_systems = set()
    seen = set() if seen is None else seen
    depth += 1
    if check_for_needed(graph.nodes(), fam2annot, func_unit):
        involved_families = set()
        local_graph = graph.copy()
        filter_local_context(local_graph, families, jaccard_threshold)
        for cc in nx.connected_components(local_graph):
            if check_for_needed(cc, fam2annot, func_unit) and not check_for_forbidden(cc, fam2annot, func_unit):
                detected_systems.add(System(system_id=0, model=func_unit.model, gene_families=cc, source=source))
                involved_families |= cc.intersection(families)
                seen.add(tuple(cc.intersection(families)))
        if (len(involved_families) < len(families) and
                len(families) - 1 >= func_unit.min_total and
                depth <= max_depth and not
                check_for_forbidden(graph.nodes(), fam2annot,
                                    func_unit)):  # Enhance by using combination only with the forbidden
            for families_combination in combinations(families, len(families) - 1):
                if families_combination not in seen:
                    seen.add(families_combination)
                    families_combination = set(families_combination)
                    if check_for_needed(families_combination, fam2annot, func_unit):
                        dissociate_family = families.difference(families_combination)
                        filtered_graph = graph.copy()
                        filtered_graph.remove_nodes_from(dissociate_family)
                        for cc in nx.connected_components(filtered_graph):
                            cc_graph = filtered_graph.subgraph(cc).copy()
                            families_in_cc_graph = families_combination & cc
                            detected_systems.update(
                                search_system_in_graph(cc_graph, families_in_cc_graph, fam2annot, func_unit, source,
                                                       jaccard_threshold, seen, depth, max_depth))
    return detected_systems


def search_system(model: Model, annot2fam: Dict[str, Set[GeneFamily]], source: str,
                  jaccard_threshold: float = 0.8, max_depth: int = 0) -> Set[System]:
    """
    Search if model system is in pangenome

    Args:
        model: Model to search in pangenome
        annot2fam: Dictionary with for each annotation a set of gene families
        source: Name of the annotation source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        max_depth: maximum depth of recursion to search for systems (default: 0)

    Returns:
        List[System]: List of systems detected in pangenome for the given model
    """
    detected_systems = set()
    for func_unit in model.func_units:
        families, fam2annot = dict_families_context(func_unit, annot2fam)
        if check_for_needed(families, fam2annot, func_unit):
            t = func_unit.max_separation + 1
            context = compute_gene_context_graph(families=families, transitive=t,
                                                 window_size=t + 1, disable_bar=True)
            for cc_graph in [context.subgraph(c).copy() for c in nx.connected_components(context)]:
                families_in_cc = ({fam for fam in families
                                   if fam2annot[fam.name] not in list(map(lambda x: x.name, func_unit.neutral))}
                                  & cc_graph.nodes)
                detected_systems.update(
                    search_system_in_graph(cc_graph, families_in_cc, fam2annot, func_unit, source, jaccard_threshold,
                                           max_depth=max_depth))
        return detected_systems


def search_systems(models: Models, pangenome: Pangenome, source: str, jaccard_threshold: float = 0.8,
                   max_depth: int = 0, threads: int = 1, lock: Lock = None, disable_bar: bool = False):
    """
    Search systems present in the pangenome for all models

    Args:
        models: Models to search in pangenomes
        pangenome: Pangenome with gene families
        source: name of the annotation source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        max_depth: maximum depth of recursion to search for systems (default: 0)
        threads: number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        disable_bar: Flag to disable progress bar
    """
    from panorama.utils import init_lock
    from panorama.format.write_binaries import write_pangenome

    annot2fam = get_annotation_to_families(pangenome=pangenome, source=source)
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=models.size, unit='model', disable=disable_bar) as progress:
            futures = []
            for model in models:
                future = executor.submit(search_system, model, annot2fam, source, jaccard_threshold, max_depth)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            detected_systems = []
            for future in futures:
                result = future.result()
                detected_systems += result

    for system in sorted(detected_systems, key=lambda x: len(x), reverse=True):
        pangenome.add_system(system)

    logging.getLogger("PANORAMA").info(f"System systems done in pangenome {pangenome.name}")

    if len(detected_systems) > 0:
        pangenome.status["systems"] = "Computed"
        logging.getLogger("PANORAMA").info(f"Write systems in pangenome {pangenome.name}")
        write_pangenome(pangenome, pangenome.file, source=source, disable_bar=disable_bar)
        logging.getLogger("PANORAMA").info(f"Systems written in pangenome {pangenome.name}")
    else:
        logging.getLogger("PANORAMA").info("No system detected")


def search_systems_in_pangenomes(models: Models, pangenomes: Pangenomes, source: str, jaccard_threshold: float = 0.8,
                                 max_depth: int = 0, threads: int = 1, lock: Lock = None, disable_bar: bool = False):
    """
    Search systems in pangenomes by multithreading on pangenomes

    Args:
        models: Models to search in pangenomes
        pangenomes: Getter object with Pangenome
        source: name of the annotation source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        max_depth: maximum depth of recursion to search for systems (default: 0)
        threads: number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        disable_bar: Disable progress bar (default: False)
    """
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Begin system systems for {pangenome.name}")
        search_systems(models, pangenome, source, jaccard_threshold, max_depth, threads, lock, disable_bar)


def launch(args):
    """
    Launch functions to detect systems in pangenomes

    Args:
        args: argument given in CLI
    """
    from panorama.utility.utility import check_models
    from panorama.format.read_binaries import load_pangenomes

    check_detection_parameters(args)
    models = check_models(args.models, disable_bar=args.disable_prog_bar)
    manager = Manager()
    lock = manager.Lock()
    need_info = {"need_annotations": True, "need_families": True, "need_metadata": True,
                 "metatypes": ["families"], "sources": [args.annotation_source]}
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_detection, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, annotation_source=args.annotation_source,
                                 systems_source=args.source, force=args.force)
    search_systems_in_pangenomes(models=models, pangenomes=pangenomes, source=args.source,
                                 jaccard_threshold=args.jaccard, max_depth=args.max_depth,
                                 threads=args.threads, lock=lock, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for systems command

    Returns:
        argparse.ArgumentParser: parser arguments for align command
    """
    parser = sub_parser.add_parser("systems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_detection(parser)
    return parser


def parser_detection(parser):
    """
    Add argument to parser for systems command

    Args:
        parser: parser for systems argument

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
    optional.add_argument('--annotation_source', required=False, type=str, default=None,
                          help="Name of the annotation source to load if different from system source")
    optional.add_argument('--jaccard', required=False, type=float, default=0.8,
                          help="minimum jaccard similarity used to filter edges between gene families. "
                               "Increasing it will improve precision but lower sensitivity a lot.")
    optional.add_argument('--max_depth', required=False, type=int, nargs='?', default=0,
                          help="Set the maximum recursion depth. Greater is the maximum longer it will take.")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
