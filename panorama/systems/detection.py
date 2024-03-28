#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import argparse
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from typing import Dict, List, Set, Tuple, FrozenSet
from multiprocessing import Manager, Lock
from collections import defaultdict

# installed libraries
import networkx as nx
from ppanggolin.genome import Gene
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import compute_gene_context_graph

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Models, Model, FuncUnit, Family
from panorama.systems.system import System
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.utils import init_lock
from panorama.format.write_binaries import write_pangenome


def check_detection_parameters(args: argparse.Namespace) -> None:
    """
    Checks the provided arguments to ensure that they are valid.

    Args:
        args: The parsed arguments.

    Raises:
        argparse.ArgumentTypeError: If jaccard is not with the correct type
    """
    args.jaccard = restricted_float(args.jaccard)
    args.annotation_sources = [args.source] if args.annotation_sources is None else args.annotation_sources


def check_pangenome_detection(pangenome: Pangenome, annotation_sources: List[str], systems_source: str,
                              force: bool = False) -> None:
    """
     Check and load pangenome information before adding annotation

    Args:
        pangenome:  Pangenome object
        annotation_sources: Source used to annotate gene famillies
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
            raise ValueError(f"Systems are already detected based on the source : {systems_source}. "
                             f"Use the --force option to erase the already computed systems.")
    if pangenome.status["metadata"]["families"] == "inFile":
        for annotation_source in annotation_sources:
            if annotation_source not in pangenome.status["metasources"]["families"]:
                raise KeyError(f"There is no metadata associate to families "
                               f"from source {annotation_source} in pangenome {pangenome.name}.")
    elif pangenome.status["metadata"]["families"] not in ["Computed", "Loaded"]:
        raise AttributeError(f"There is no metadata associate to families in your pangenome {pangenome.name}. "
                             "Please see the command annotation before to detect systems")


def get_annotation_to_families(pangenome: Pangenome, sources: List[str]) -> Dict[str, Set[GeneFamily]]:
    """
    Get for each annotation a set of families with this annotation

    Args:
        pangenome: Pangenome with gene families
        sources: list of annotation source name

    Returns:
        Dictionary with for each annotation a set of gene families
    """
    annot2fam = defaultdict(set)
    for source in sources:
        for gf in pangenome.gene_families:
            metadata = gf.get_metadata_by_source(source)
            if metadata is not None:
                for meta in metadata:
                    annot2fam[meta.protein_name].add(gf)
    return annot2fam


def dict_families_context(model: Model, annot2fam: Dict[str, Set[GeneFamily]]) \
        -> Tuple[Set[GeneFamily], Dict[str, Set[Family]]]:
    """
    Recover all families in the function unit

    Args:
        model: model containing the families
        annot2fam: dictionary of annotated families

    Returns:
        The set of families of interest in the functional unit and
        a dictionary which link families to their annotation for the functional unit
    """
    families = set()
    fam2annot = defaultdict(set)
    for fam_model in model.families:
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
    def get_family_genes_in_organisms(family: GeneFamily) -> Set[Gene]:
        if family.name not in fam2genes_in_orgs:
            family_genes_in_orgs = {gene for gene in family.genes if gene.organism in organisms_of_interest}
            fam2genes_in_orgs[family.name] = family_genes_in_orgs
        else:
            family_genes_in_orgs = fam2genes_in_orgs[family.name]
        return family_genes_in_orgs

    # Compute local jaccard
    organisms_of_interest = set()
    if len(families) > 1:
        subgraph = graph.subgraph(families)
        for _, _, data in subgraph.edges(data=True):
            organisms_of_interest |= data['genomes']

    elif len(families) == 1:
        organisms_of_interest = set(list(families)[0].organisms)
    else:
        raise AssertionError("No families of interest")

    fam2genes_in_orgs = {}
    edges2remove = set()
    for f1, f2, data in graph.edges(data=True):
        f1_genes_in_orgs = get_family_genes_in_organisms(f1)
        f2_genes_in_orgs = get_family_genes_in_organisms(f2)
        f1_gene_proportion = len(data['genes'][f1]) / len(f1_genes_in_orgs) if len(f1_genes_in_orgs) > 0 else 0
        f2_gene_proportion = len(data['genes'][f2]) / len(f2_genes_in_orgs) if len(f2_genes_in_orgs) > 0 else 0

        data['f1'] = f1.name
        data['f2'] = f2.name
        data['f1_jaccard_gene'] = f1_gene_proportion
        data['f2_jaccard_gene'] = f2_gene_proportion

        if not ((f1_gene_proportion >= jaccard_threshold) and (f2_gene_proportion >= jaccard_threshold)):
            edges2remove.add((f1, f2))

    # Filter cc
    # filter_flag = f'is_jaccard_gene_>_{jaccard_threshold}'
    #
    # edges_to_remove = [(n, v) for n, v, d in graph.edges(data=True) if not d[filter_flag]]
    graph.remove_edges_from(edges2remove)


def check_for_forbidden(families: Set[GeneFamily], fam2annot: Dict[str, Set[Family]],
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
        annots = fam2annot.get(node.name, [])
        for annot in annots:
            if annot.presence == 'forbidden' and annot.name in forbidden_list:  # if node is forbidden
                count_forbidden += 1
                forbidden_list.remove(annot.name)
                if count_forbidden > func_unit.max_forbidden:
                    return True
                break
    return False


def check_for_needed(families: Set[GeneFamily], fam2annot: Dict[str, Set[Family]],
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
        annots = fam2annot.get(node.name, [])
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


def get_subcombinations(combi: Set[GeneFamily],
                        combinations: List[FrozenSet[GeneFamily]]) -> List[FrozenSet[GeneFamily]]:
    """Remove combinations from combinations and all subcombination
    """
    index = 0
    remove_combinations = []
    while index < len(combinations):
        if combinations[index].issubset(combi):
            remove_combinations.append(combinations.pop(index))
        else:
            index += 1
    return remove_combinations


def search_system_in_cc(graph: nx.Graph, families: Set[GeneFamily], fam2annot: Dict[str, Set[Family]],
                        func_unit: FuncUnit, source: str, jaccard_threshold: float = 0.8,
                        combinations: List[FrozenSet[GeneFamily]] = None) -> List[System]:
    """Search systems corresponding to model in a graph

    Args:
        graph: A graph with families in one connected component
        families: A set of families that code for the searched model
        fam2annot: A dictionary to link gene families to their annotation
        func_unit: One functional unit corresponding to the model
        source: Name of the source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)

    Returns:
        A set of all detected systems in the graph
    """
    detected_systems = []
    while len(combinations) > 0:
        families_combination = set(combinations.pop(0))
        if check_for_needed(families_combination, fam2annot, func_unit):
            if not check_for_forbidden(families_combination, fam2annot, func_unit):
                local_graph = graph.copy()
                local_graph.remove_nodes_from(families.difference(families_combination))
                filter_local_context(local_graph, families_combination, jaccard_threshold)
                for cc in sorted(nx.connected_components(local_graph), key=len, reverse=True):
                    cc: Set[GeneFamily]
                    if check_for_needed(cc, fam2annot, func_unit) and not check_for_forbidden(cc, fam2annot, func_unit):
                        detected_systems.append(System(system_id=0, model=func_unit.model,
                                                       gene_families=cc, source=source))
                        get_subcombinations(cc.intersection(families_combination), combinations)
        else:
            # We can remove all sub combinations because the bigger one does not have the needed
            get_subcombinations(families_combination, combinations)

    return detected_systems


def search_system_in_context(graph: nx.Graph, families: Set[GeneFamily], fam2annot: Dict[str, Set[Family]],
                             func_unit: FuncUnit, source: str, jaccard_threshold: float = 0.8,
                             combinations: List[FrozenSet[GeneFamily]] = None) -> List[System]:
    detected_systems = []
    for cc_graph in [graph.subgraph(c).copy() for c in sorted(nx.connected_components(graph),
                                                              key=len, reverse=True)]:
        families_in_cc = ({fam for fam in families
                           if fam2annot[fam.name] not in list(map(lambda x: x.name, func_unit.neutral))}
                          & cc_graph.nodes)
        if len(families_in_cc) > 0:
            combinations_in_cc = get_subcombinations(families_in_cc, combinations)
            new_detected_systems = search_system_in_cc(cc_graph, families_in_cc, fam2annot,
                                                       func_unit, source, jaccard_threshold,
                                                       combinations=combinations_in_cc)
            detected_systems += new_detected_systems
    return detected_systems


def search_system(model: Model, annot2fam: Dict[str, Set[GeneFamily]], source: str,
                  jaccard_threshold: float = 0.8) -> List[System]:
    """
    Search if model system is in pangenome

    Args:
        model: Model to search in pangenome
        annot2fam: Dictionary with for each annotation a set of gene families
        source: Name of the annotation source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)

    Returns:
        List[System]: List of systems detected in pangenome for the given model
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for model {model.name}")

    detected_systems = []
    families, fam2annot = dict_families_context(model, annot2fam)
    for func_unit in model.func_units:
        fu_families = {fam for fam in families if any(annot in func_unit.families for annot in fam2annot[fam.name])}
        if check_for_needed(fu_families, fam2annot, func_unit):
            t = func_unit.max_separation + 1
            context, combinations2orgs = compute_gene_context_graph(families=families, transitive=t,
                                                                    window_size=t + 1, disable_bar=True)
            combinations = sorted(set(combinations2orgs.keys()), key=len, reverse=True)
            detected_systems += search_system_in_context(context, families, fam2annot, func_unit,
                                                         source, jaccard_threshold, combinations)
    return detected_systems


def search_systems(models: Models, pangenome: Pangenome, source: str, annotation_sources: List[str],
                   jaccard_threshold: float = 0.8, threads: int = 1, lock: Lock = None,
                   disable_bar: bool = False):
    """
    Search systems present in the pangenome for all models

    Args:
        models: Models to search in pangenomes
        pangenome: Pangenome with gene families
        source: name of the source for the system
        annotation_sources: list of the annotation source for the families
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        threads: number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        disable_bar: Flag to disable progress bar
    """

    annot2fam = get_annotation_to_families(pangenome=pangenome, sources=annotation_sources)
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=models.size, unit='model', disable=disable_bar) as progress:
            futures = []
            for model in models:
                future = executor.submit(search_system, model, annot2fam, source, jaccard_threshold)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

                detected_systems = []
                for future in futures:
                    result = future.result()
                    detected_systems += result

    for system in sorted(detected_systems, key=lambda x: len(x), reverse=True):
        pangenome.add_system(system)

    logging.getLogger("PANORAMA").info(f"Systems prediction done in pangenome {pangenome.name}")

    if len(detected_systems) > 0:
        # print(len(detected_systems))
        pangenome.status["systems"] = "Computed"
        logging.getLogger("PANORAMA").info(f"Write systems in pangenome {pangenome.name}")
        write_pangenome(pangenome, pangenome.file, source=source, disable_bar=disable_bar)
        logging.getLogger("PANORAMA").info(f"Systems written in pangenome {pangenome.name}")
    else:
        logging.getLogger("PANORAMA").info("No system detected")


def search_systems_in_pangenomes(models: Models, pangenomes: Pangenomes, source: str, annotation_sources: List[str],
                                 jaccard_threshold: float = 0.8, threads: int = 1,
                                 lock: Lock = None, disable_bar: bool = False):
    """
    Search systems in pangenomes by multithreading on pangenomes

    Args:
        models: Models to search in pangenomes
        pangenomes: Getter object with Pangenome
        source: name of the source for the system
        annotation_sources: list of the annotation source for the families
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        threads: number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        disable_bar: Disable progress bar (default: False)
    """
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Begin systems searching for {pangenome.name}")
        search_systems(models, pangenome, source, annotation_sources, jaccard_threshold,
                       threads, lock, disable_bar)


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
                 "metatypes": ["families"], "sources": args.annotation_sources}
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_detection, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, annotation_sources=args.annotation_sources,
                                 systems_source=args.source, force=args.force)
    search_systems_in_pangenomes(models=models, pangenomes=pangenomes, source=args.source,
                                 annotation_sources=args.annotation_sources, jaccard_threshold=args.jaccard,
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
    optional.add_argument('--annotation_sources', required=False, type=str, default=None, nargs='+',
                          help="Name of the annotation sources to load if different from system source. "
                               "Could be more than one, separated by space.")
    optional.add_argument('--jaccard', required=False, type=float, default=0.8,
                          help="minimum jaccard similarity used to filter edges between gene families. "
                               "Increasing it will improve precision but lower sensitivity a lot.")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
