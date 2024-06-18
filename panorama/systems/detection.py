#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to detect biological systems in pangenomes
"""

# default libraries
from __future__ import annotations
import argparse
import time
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from typing import Dict, Iterable, List, Set, Tuple, FrozenSet
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
from panorama.format.write_binaries import write_pangenome, erase_pangenome


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


def check_pangenome_detection(pangenome: Pangenome, metadata_sources: List[str], systems_source: str,
                              force: bool = False) -> None:
    """
     Check and load pangenome information before adding annotation

    Args:
        pangenome:  Pangenome object
        metadata_sources: Sources used to associate annotation to gene famillies
        systems_source: Source used to detect system
        force: Force to erase pangenome systems from source

    Raises:
        KeyError: Provided annotation source is not in the pangenome
        Exception: If System already exists in the Pangenome
        AttributeError: If there is no metadata associated to families
    """
    if pangenome.status["systems"] == "inFile" and systems_source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=systems_source)
        else:
            raise ValueError(f"Systems are already detected based on the source : {systems_source}. "
                             f"Use the --force option to erase the already computed systems.")
    if pangenome.status["metadata"]["families"] == "inFile":
        for annotation_source in metadata_sources:
            if annotation_source not in pangenome.status["metasources"]["families"]:
                raise KeyError(f"There is no metadata associate to families "
                               f"from source {annotation_source} in pangenome {pangenome.name}.")
    elif pangenome.status["metadata"]["families"] not in ["Computed", "Loaded"]:
        raise AttributeError(f"There is no metadata associate to families in your pangenome {pangenome.name}. "
                             "Please see the command annotation before to detect systems")


def get_metadata_to_families(pangenome: Pangenome, sources: Iterable[str]) -> Dict[str, Dict[str, Set[GeneFamily]]]:
    """
    Get for each metadata a set of families with the same protein_name value in the metadata

    Args:
        pangenome: Pangenome with gene families
        sources: list of metadata source name

    Returns:
        Dictionary with for each metadata source a dictionary with for each metadata a set of gene families
    """
    meta2fam = {source: defaultdict(set) for source in sources}
    for source in sources:
        for gf in pangenome.gene_families:
            metadata = gf.get_metadata_by_source(source)
            if metadata is not None:
                for meta in metadata.values():
                    meta2fam[source][meta.protein_name].add(gf)
                    if "secondary_name" in meta.fields:
                        for secondary_name in meta.secondary_name.split(','):
                            meta2fam[source][secondary_name].add(gf)
    return meta2fam


def dict_families_context(model: Model, annot2fam: Dict[str, Dict[str, Set[GeneFamily]]]) \
        -> Tuple[Set[GeneFamily], Dict[str, Set[Family]], Dict[str, str]]:
    """
    Retrieve all families that code for a family in the model

    Args:
        model: model containing the families
        annot2fam: dictionary of annotated families

    Returns:
        The set of gene families of interest in the functional unit and
        a dictionary which link families to their annotation for the functional unit
    """
    gene_families = set()
    gf2fam = defaultdict(set)
    fam2source = {}
    for fam_model in model.families:
        for source, annotation2families in annot2fam.items():
            if fam_model.name in annotation2families:
                for gf in annotation2families[fam_model.name]:
                    gene_families.add(gf)
                    gf2fam[gf.name].add(fam_model)
                    if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                        logging.getLogger("PANORAMA").warning("2 annotation have the same protein name for different "
                                                              "sources. First source encountered will be used.")
                    else:
                        fam2source[fam_model.name] = source

        for exchangeable in fam_model.exchangeable:
            for source, annotation2families in annot2fam.items():
                if exchangeable in annotation2families:
                    for gf in annotation2families[exchangeable]:
                        gene_families.add(gf)
                        gf2fam[gf.name].add(fam_model)
                        if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                            logging.getLogger("PANORAMA").warning(
                                "2 annotation have the same protein name for different "
                                "sources. First source encountered will be used.")
                        else:
                            fam2source[fam_model.name] = source
    return gene_families, gf2fam, fam2source


def filter_local_context(graph: nx.Graph, families: Set[GeneFamily], jaccard_threshold: float = 0.8) -> None:
    """
    Filtered a graph based on a local jaccard index

    Args:
        graph: A sub pangenome graph
        families: list of families that code for the system
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
    """

    def get_family_genes_in_organisms(family: GeneFamily) -> Set[Gene]:
        """
        Get the set of genes of the gene family that are in organisms of interest

        Args:
            family: The gene family to search genes

        Returns:
            The set of genes of the gene family that are in organisms of interest
        """
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


def check_for_forbidden(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[str, Set[Family]],
                        func_unit: FuncUnit) -> bool:
    """
    Check if there is forbidden in families.

    Args:
        gene_families: set of families
        gene_fam2mod_fam: link between gene family and annotation
        func_unit: functional unit

    Returns:
        Boolean true if forbidden condition are encountered in families
    """

    count_forbidden = 0
    forbidden_list = list(map(lambda x: x.name, func_unit.forbidden))

    def get_number_of_mod_fam(gene_family: GeneFamily) -> int:
        """Get the number of model family associated to the gene family.

        Args:
            gene_family: The gene family of interest

        Returns:
            The number of model family associated to the gene family.
        """
        model_families = gene_fam2mod_fam.get(gene_family.name)
        if model_families is not None:
            return len(model_families)
        else:
            return 0

    for node in sorted(gene_families, key=lambda n: get_number_of_mod_fam(n)):
        for family in gene_fam2mod_fam[node.name]:
            if family.presence == 'forbidden' and family.name in forbidden_list:  # if node is forbidden
                count_forbidden += 1
                forbidden_list.remove(family.name)
                if count_forbidden > func_unit.max_forbidden:
                    return True
                break
    return False


def check_for_needed(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[str, Set[Family]],
                     mod_fam2meta_source: Dict[str, str], func_unit: FuncUnit
                     ) -> Tuple[bool, Dict[GeneFamily, Tuple[str, int]]]:
    """
    Check if presence/absence rules for needed families are respected

    Args:
        gene_families: set of families
        gene_fam2mod_fam: link between gene family and model families
        mod_fam2meta_source: link between model families and metadata source
        func_unit: functional unit

    Returns:
        Boolean true if parameter respected
    """
    count_mandatory, count_total = (0, 0)
    mandatory_list, accessory_list = (list(map(lambda x: x.name, func_unit.mandatory)),
                                      list(map(lambda x: x.name, func_unit.accessory)))

    def get_number_of_mod_fam(gene_family: GeneFamily) -> int:
        """Get the number of model family associated to the gene family.

        Args:
            gene_family: The gene family of interest

        Returns:
            The number of model family associated to the gene family.
        """
        model_families = gene_fam2mod_fam.get(gene_family.name)
        if model_families is not None:
            return len(model_families)
        else:
            return 0

    families2meta_info = {}
    for node in sorted(gene_families, key=lambda n: get_number_of_mod_fam(n)):
        for family in gene_fam2mod_fam[node.name]:
            if family.presence == 'mandatory' and family.name in mandatory_list:  # if node is mandatory
                for meta_id, metadata in node.get_metadata_by_source(mod_fam2meta_source[family.name]).items():
                    if metadata.protein_name == family.name:
                        families2meta_info[node] = (mod_fam2meta_source[family.name], meta_id)
                    elif "secondary_name" in metadata.fields:
                        if family.name in metadata.secondary_name.split(","):
                            families2meta_info[node] = (mod_fam2meta_source[family.name], meta_id)
                count_mandatory += 1
                count_total += 1
                mandatory_list.remove(family.name)
                break
            elif family.presence == 'accessory' and family.name in accessory_list:  # if node is accessory
                for meta_id, metadata in node.get_metadata_by_source(mod_fam2meta_source[family.name]).items():
                    if metadata.protein_name == family.name:
                        families2meta_info[node] = (mod_fam2meta_source[family.name], meta_id)
                    elif "secondary_name" in metadata.fields:
                        if family.name in metadata.secondary_name.split(","):
                            families2meta_info[node] = (mod_fam2meta_source[family.name], meta_id)
                count_total += 1
                accessory_list.remove(family.name)
                break

    if (count_mandatory >= func_unit.min_mandatory or func_unit.min_mandatory == -1) and \
            (count_total >= func_unit.min_total or func_unit.min_total == -1):
        if (func_unit.max_mandatory >= count_mandatory or func_unit.max_mandatory == -1) and \
                (func_unit.max_total >= count_total or func_unit.max_total == -1):
            return True, families2meta_info
        else:
            return False, families2meta_info
    else:
        return False, families2meta_info


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


def search_system_in_cc(graph: nx.Graph, families: Set[GeneFamily], gene_fam2mod_fam: Dict[str, Set[Family]],
                        mod_fam2meta_source: Dict[str, str], func_unit: FuncUnit, source: str,
                        jaccard_threshold: float = 0.8, combinations: List[FrozenSet[GeneFamily]] = None
                        ) -> List[System]:
    """Search systems corresponding to model in a graph

    Args:
        graph: A graph with families in one connected component
        families: A set of families that code for the searched model
        gene_fam2mod_fam: A dictionary to link gene families to their annotation
        mod_fam2meta_source: link between model families and metadata source
        func_unit: One functional unit corresponding to the model
        source: Name of the source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        combinations: List of families combination known to exist in genomes.

    Returns:
        A set of all detected systems in the graph
    """
    detected_systems = []
    while len(combinations) > 0:
        families_combination = set(combinations.pop(0))
        check_needed, _ = check_for_needed(families_combination, gene_fam2mod_fam, mod_fam2meta_source, func_unit)
        if check_needed:
            if not check_for_forbidden(families_combination, gene_fam2mod_fam, func_unit):
                local_graph = graph.copy()
                local_graph.remove_nodes_from(families.difference(families_combination))
                filter_local_context(local_graph, families_combination, jaccard_threshold)
                for cc in sorted(nx.connected_components(local_graph), key=len, reverse=True):
                    cc: Set[GeneFamily]
                    check_needed, fam2metainfo = check_for_needed(cc, gene_fam2mod_fam, mod_fam2meta_source, func_unit)
                    if check_needed and not check_for_forbidden(cc, gene_fam2mod_fam, func_unit):
                        detected_systems.append(System(system_id=0, model=func_unit.model, source=source,
                                                       gene_families=cc, families_to_metainfo=fam2metainfo)
                                                )
                        get_subcombinations(cc.intersection(families_combination), combinations)
        else:
            # We can remove all sub combinations because the bigger one does not have the needed
            get_subcombinations(families_combination, combinations)

    return detected_systems


def search_system_in_context(graph: nx.Graph, families: Set[GeneFamily], gene_fam2mod_fam: Dict[str, Set[Family]],
                             mod_fam2meta_source: Dict[str, str], func_unit: FuncUnit, source: str,
                             jaccard_threshold: float = 0.8, combinations: List[FrozenSet[GeneFamily]] = None
                             ) -> List[System]:
    """
    Search systems corresponding to model in a pangenomic context

    Args:
        graph: A graph with families in a pangenomic context
        families: A set of families that code for the searched model
        gene_fam2mod_fam: A dictionary to link gene families to model families
        mod_fam2meta_source: link between model families and metadata source
        func_unit: One functional unit corresponding to the model
        source: Name of the source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        combinations: Known existing combination of families that code for the searched model

    Returns:
        List of detected systems in the pangenomic context
    """
    detected_systems = []
    for cc_graph in [graph.subgraph(c).copy() for c in sorted(nx.connected_components(graph),
                                                              key=len, reverse=True)]:
        families_in_cc, neutral_families = get_functional_unit_gene_families(func_unit, families, gene_fam2mod_fam)
        families_in_cc &= cc_graph.nodes

        if len(families_in_cc) > 0:
            combinations_in_cc = get_subcombinations(families_in_cc, combinations)
            families_in_cc |= neutral_families
            new_detected_systems = search_system_in_cc(cc_graph, families_in_cc, gene_fam2mod_fam, mod_fam2meta_source,
                                                       func_unit, source, jaccard_threshold, combinations_in_cc)
            detected_systems += new_detected_systems
    return detected_systems


def get_functional_unit_gene_families(func_unit: FuncUnit, gene_families: Set[GeneFamily],
                                      gene_fam2mod_fam:  Dict[str, Set[Family]]
                                      ) -> Tuple[Set[GeneFamily], Set[GeneFamily]]:
    """
    Get the gene families that might be in the functional unit

    Args:
        func_unit: The functional unit to consider
        gene_families: Set of gene families that might be in the system
        gene_fam2mod_fam: link between gene family and model family

    Returns:
        gene families that code for the functional unit, the neutral families.
    """
    fu_families = set()
    neutral_families = set()
    for gf in gene_families:
        for family in gene_fam2mod_fam[gf.name]:
            if family in func_unit.neutral:
                neutral_families.add(gf)
                break
            elif family in func_unit.families:
                fu_families.add(gf)
                break
    return fu_families, neutral_families


def search_system(model: Model, meta2fam: Dict[str, Dict[str, Set[GeneFamily]]], source: str,
                  jaccard_threshold: float = 0.8) -> List[System]:
    """
    Search if model system is in pangenome

    Args:
        model: Model to search in pangenome
        meta2fam: Dictionary with for each metadata source a dictionary with a metadate linked to a set of gene_families
        source: Name of the annotation source
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene gene_families (default: 0.8)

    Returns:
        List[System]: List of systems detected in pangenome for the given model
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for model {model.name}")
    begin = time.time()
    detected_systems = []
    gene_families, gf2fam, fam2source = dict_families_context(model, meta2fam)

    for func_unit in model.func_units:
        fu_families = set()
        fu_families.update(*get_functional_unit_gene_families(func_unit, gene_families, gf2fam))
        check_needed, _ = check_for_needed(fu_families, gf2fam, fam2source, func_unit)
        if check_needed:
            t = func_unit.max_separation + 1
            context, combinations2orgs = compute_gene_context_graph(families=fu_families, transitive=t,
                                                                    window_size=t + 1, disable_bar=True)
            combinations = sorted(set(combinations2orgs.keys()), key=len, reverse=True)
            detected_systems += search_system_in_context(context, fu_families, gf2fam, fam2source, func_unit, source,
                                                         jaccard_threshold, combinations)
    logging.getLogger("PANORAMA").debug(f"Done search for model {model.name} in {time.time() - begin} seconds")
    return detected_systems


def search_systems(models: Models, pangenome: Pangenome, source: str, metadata_sources: List[str],
                   jaccard_threshold: float = 0.8, threads: int = 1, lock: Lock = None,
                   disable_bar: bool = False):
    """
    Search systems present in the pangenome for all models

    Args:
        models: Models to search in pangenomes
        pangenome: Pangenome with gene families
        source: name of the source for the system
        metadata_sources: list of the metadata source for the families
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        threads: number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        disable_bar: Flag to disable progress bar
    """

    meta2fam = get_metadata_to_families(pangenome=pangenome, sources=metadata_sources)
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=models.size, unit='model', disable=disable_bar) as progress:
            futures = []
            for model in models:
                future = executor.submit(search_system, model, meta2fam, source, jaccard_threshold)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

                detected_systems = []
                for future in futures:
                    result = future.result()
                    detected_systems += result

    for system in sorted(detected_systems, key=lambda x: (len(x.model.canonical), -len(x))):
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


def search_systems_in_pangenomes(models: Models, pangenomes: Pangenomes, source: str, metadata_sources: List[str],
                                 jaccard_threshold: float = 0.8, threads: int = 1,
                                 lock: Lock = None, disable_bar: bool = False):
    """
    Search systems in pangenomes by multithreading on pangenomes

    Args:
        models: Models to search in pangenomes
        pangenomes: Getter object with Pangenome
        source: name of the source for the system
        metadata_sources: list of the metadata source for the families
        jaccard_threshold: minimum jaccard similarity used to filter edges between gene families (default: 0.8)
        threads: number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        disable_bar: Disable progress bar (default: False)
    """
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Begin systems searching for {pangenome.name}")
        search_systems(models, pangenome, source, metadata_sources, jaccard_threshold,
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
                                 disable_bar=args.disable_prog_bar, metadata_sources=args.annotation_sources,
                                 systems_source=args.source, force=args.force)
    search_systems_in_pangenomes(models=models, pangenomes=pangenomes, source=args.source,
                                 metadata_sources=args.annotation_sources, jaccard_threshold=args.jaccard,
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
