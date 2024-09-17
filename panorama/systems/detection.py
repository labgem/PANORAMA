#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to detect biological systems in pangenomes.
"""

# default libraries
from __future__ import annotations
import argparse
import itertools
import time
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from typing import Dict, List, Set, Tuple, FrozenSet, Union
from multiprocessing import Manager, Lock, get_context

# installed libraries
import pandas as pd
import networkx as nx
from tqdm import tqdm
from ppanggolin.genome import Organism
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import compute_gene_context_graph
from ppanggolin.metadata import Metadata

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.utils import (filter_local_context, check_for_forbidden_families, get_metadata_to_families,
                                    get_gfs_matrix_combination, dict_families_context, find_combinations,
                                    check_needed_families)
from panorama.systems.models import Models, Model, FuncUnit, Family
from panorama.systems.system import System, SystemUnit
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.utils import init_lock
from panorama.format.write_binaries import write_pangenome, erase_pangenome

import sys
import pickle
sys.setrecursionlimit(100000)


def check_detection_args(args: argparse.Namespace) -> Dict[str, Union[bool, str, List[str]]]:
    """
    Checks and processes the provided arguments to ensure they are valid.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Raises:
        argparse.ArgumentTypeError: If 'jaccard' is not a restricted float.

    Returns:
        dict: A dictionary indicating the required families and sources.
    """
    args.jaccard = restricted_float(args.jaccard)
    args.annotation_sources = [args.source] if args.annotation_sources is None else args.annotation_sources
    return {"need_annotations": True, "need_families": True, "need_metadata": True,
            "metatypes": ["families"], "sources": args.annotation_sources}


def check_pangenome_detection(pangenome: Pangenome, metadata_sources: List[str], systems_source: str,
                              force: bool = False) -> None:
    """
    Checks and loads pangenome information before adding families.

    Args:
        pangenome (Pangenome): Pangenome object to be checked.
        metadata_sources (list of str): Sources used to associate families to gene families.
        systems_source (str): Source used to detect systems.
        force (bool, optional): If True, forces the erasure of pangenome systems from the source.

    Raises:
        KeyError: If the provided annotation source is not in the pangenome.
        ValueError: If systems are already detected based on the source and 'force' is not used.
        AttributeError: If there is no metadata associated with families.
    """
    if pangenome.status["systems"] == "inFile" and systems_source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=systems_source)
        else:
            raise ValueError(f"Systems are already detected based on the source: {systems_source}. "
                             f"Use the --force option to erase the already computed systems.")
    if pangenome.status["metadata"]["families"] == "inFile":
        for annotation_source in metadata_sources:
            if annotation_source not in pangenome.status["metasources"]["families"]:
                raise KeyError(f"There is no metadata associated with families "
                               f"from source {annotation_source} in pangenome {pangenome.name}.")
    elif pangenome.status["metadata"]["families"] not in ["Computed", "Loaded"]:
        raise AttributeError(f"There is no metadata associated with families in your pangenome {pangenome.name}. "
                             "Please see the command annotation before detecting systems")


def get_subcombinations(combi: Set[GeneFamily], combinations: List[FrozenSet[GeneFamily]]
                        ) -> List[FrozenSet[GeneFamily]]:
    """
    Removes a combination and all its sub-combinations from a list of combinations.

    Args:
        combi (set of GeneFamily): Combination to be removed.
        combinations (list of FrozenSet[GeneFamily]): List of combinations to be filtered.

    Returns:
        List[FrozenSet[GeneFamily]]: List of removed combinations.
    """
    index = 0
    remove_combinations = []
    while index < len(combinations):
        if combinations[index].issubset(combi):
            remove_combinations.append(combinations.pop(index))
        else:
            index += 1
    return remove_combinations


def get_gf2metainfo(fam2gf: Dict[str, str], mod_fam2meta_source: Dict[str, str],
                    mandatory_gfs2metadata: Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                    accessory_gfs2metadata: Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                    name2gf: Dict[str, GeneFamily]) -> Dict[GeneFamily, Tuple[str, int]]:
    """
    Build a dictionary of gene families link to the metadata coding for the unit.

    Returns:
        Dict[GeneFamily, Tuple[str, int]]: Dictionary with gene families as key and a tuple of metadata source and id as value
    """
    gf2meta_info = {}
    for gf_name, fam in fam2gf.items():
        if fam in mandatory_gfs2metadata:
            gf2metadata = mandatory_gfs2metadata[fam]
        else:
            gf2metadata = accessory_gfs2metadata[fam]
        gf = name2gf[gf_name]
        meta_id, _ = gf2metadata[gf]
        gf2meta_info[gf] = (mod_fam2meta_source[fam], meta_id)
    return gf2meta_info


def search_fu_in_cc(graph: nx.Graph, gfs: Set[GeneFamily], gf_name2fam_name: Dict[str, Set[str]],
                    mod_fam2meta_source: Dict[str, str], func_unit: FuncUnit, source: str, matrix: pd.DataFrame,
                    mandatory_gfs2metadata: Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                    accessory_gfs2metadata: Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                    combinations: List[FrozenSet[GeneFamily]],
                    combinations2orgs: Dict[FrozenSet[GeneFamily], Set[Organism]],
                    jaccard_threshold: float = 0.8, i: int = 0) -> Set[SystemUnit]:
    """
    Searches for functional unit corresponding to a model in a graph.

    Args:
        graph (nx.Graph): A graph with families in one connected component.
        gfs (Set[GeneFamily]): A set of families that code for the searched model.
        gf_name2fam_name (Dict[str, Set[Family]]): Dictionary linking gene families to their families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking model families to metadata sources.
        func_unit (FuncUnit): One functional unit corresponding to the model.
        source (str): Name of the source.
        matrix (pd.DataFrame): Dataframe containing association between gene families and unit families.
        mandatory_gfs2metadata (Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]): Dictionary linking gene families to metadata for mandatory families
        accessory_gfs2metadata (Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]): Dictionary linking gene families to metadata for accessory families
        combinations: Existing combination of gene families in organisms.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
        combinations2orgs (Dict[FrozenSet[GeneFamily], Set[Organism]]): the combination of gene families corresponding to the context that exist in at least one genome

    Returns:
        Set[SystemUnit]: Set of all detected functional unit in the graph.
    """

    detected_su = set()
    mandatory_gfs = {gf for gf2metadata in mandatory_gfs2metadata.values() for gf in gf2metadata.keys() if
                     gf in gfs}
    accessory_gfs = {gf for gf2metadata in accessory_gfs2metadata.values() for gf in gf2metadata.keys() if
                     gf in gfs}
    # with VizTracer(tracer_entries=5000000, output_file=f"detect{i}.json") as tracer:
    name2gf = {gf.name: gf for gf in gfs}
    while len(combinations) > 0:
        # Continue if there is combination
        context_combination = combinations.pop(0)  # combination found when searching contex
        organisms_of_interest = combinations2orgs[context_combination]
        context_combination = set(context_combination)
        context_gfs_name = {gf.name for gf in context_combination} & {gf.name for gf in mandatory_gfs | accessory_gfs}
        if len(context_gfs_name) >= func_unit.min_total:
            filtered_matrix = matrix[list(context_gfs_name)]
            solutions = check_needed_families(filtered_matrix, func_unit)
            if solutions:
                # Keep only mandatory and accessory of current combination and eventual forgiven and neutral families.
                node2remove = gfs.difference(context_combination)
                filtered_graph = filter_local_context(graph, organisms_of_interest, node2remove,
                                                      jaccard_threshold=jaccard_threshold)
                for cc in sorted(nx.connected_components(filtered_graph), key=len, reverse=True):
                    cc: Set[GeneFamily]
                    if len(cc.intersection(context_combination)) >= func_unit.min_total:
                        families_in_cc = context_combination.intersection(cc)
                        final_matrix = filtered_matrix[list({gf.name for gf in families_in_cc} & context_gfs_name)]
                        if not check_for_forbidden_families(families_in_cc, gf_name2fam_name, func_unit):
                            working_combs = find_combinations(final_matrix, func_unit)
                            for working_comb, fam2gf in working_combs:
                                fam2meta_info = get_gf2metainfo(fam2gf, mod_fam2meta_source, mandatory_gfs2metadata,
                                                                accessory_gfs2metadata, name2gf)
                                detected_su.add(SystemUnit(functional_unit=func_unit, source=source, gene_families=cc,
                                                           gene_families_to_metainfo=fam2meta_info))
                                get_subcombinations(cc.intersection(context_combination), combinations)
            else:
                get_subcombinations(context_combination, combinations)
        else:
            # We can remove all sub-combinations because the bigger one does not have the needed
            get_subcombinations(context_combination, combinations)

    return detected_su


def search_unit_in_context(graph: nx.Graph, families: Set[GeneFamily], gf2fam: Dict[str, Set[str]],
                           mod_fam2meta_source: Dict[str, str], matrix: pd.DataFrame,
                           mandatory_gfs2metadata: Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                           accessory_gfs2metadata: Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]],
                           func_unit: FuncUnit, source: str,
                           combinations2orgs: Dict[FrozenSet[GeneFamily], Set[Organism]] = None,
                           jaccard_threshold: float = 0.8) -> Set[SystemUnit]:
    """
    Searches for system unit corresponding to a model in a pangenomic context.

    Args:
        graph (nx.Graph): A graph with families in a pangenomic context.
        families (Set[GeneFamily]): A set of families that code for the searched model.
        gf2fam (Dict[str, Set[Family]]): Dictionary linking gene families to model families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking model families to metadata sources.
        func_unit (FuncUnit): One functional unit corresponding to the model.
        source (str): Name of the source.
        matrix (pd.DataFrame): Dataframe containing association between gene families and unit families.
        mandatory_gfs2metadata (Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]): Dictionary linking gene families to metadata for mandatory families
        accessory_gfs2metadata (Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]): Dictionary linking gene families to metadata for accessory families
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
        combinations2orgs (Dict[FrozenSet[GeneFamily], Set[Organism]]): the combination of gene families corresponding to the context that exist in at least one genome

    Returns:
        Set[SystemUnit]: Set of detected systems in the pangenomic context.
    """
    detected_fu = set()
    families_in_cc, neutral_families = get_functional_unit_gene_families(func_unit, families, gf2fam)
    combinations = sorted(combinations2orgs.keys(), key=lambda fs: (len(fs), sorted(fs)), reverse=True)
    for i, cc in enumerate(sorted(nx.connected_components(graph), key=len, reverse=True)):
        cc: Set[GeneFamily]
        fam_in_cc = families_in_cc & cc

        if len(families_in_cc) > 0:
            combinations_in_cc = get_subcombinations(fam_in_cc, combinations)
            n_fam = neutral_families & cc
            fam_in_cc |= n_fam
            cc_graph = graph.subgraph(cc).copy()
            to = time.time()
            new_detected_fu = search_fu_in_cc(cc_graph, fam_in_cc, gf2fam, mod_fam2meta_source, func_unit, source,
                                              matrix, mandatory_gfs2metadata, accessory_gfs2metadata,
                                              combinations_in_cc, combinations2orgs, jaccard_threshold, i)
            if time.time() - to > 1:
                print(f"search: {time.time() - to} seconds")
            detected_fu |= new_detected_fu
    return detected_fu


def get_functional_unit_gene_families(func_unit: FuncUnit, gene_families: Set[GeneFamily],
                                      gf2fam: Dict[str, Set[str]]) -> Tuple[Set[GeneFamily], Set[GeneFamily]]:
    """
    Retrieves the gene families that might be in the functional unit.

    Args:
        func_unit (FuncUnit): The functional unit to consider.
        gene_families (Set[GeneFamily]): Set of gene families that might be in the system.
        gf2fam (Dict[str, Set[Family]]): Dictionary linking gene families to model families.

    Returns:
        tuple: A tuple containing:
            - Set[GeneFamily]: Gene families that code for the functional unit.
            - Set[GeneFamily]: The neutral families.
    """
    fu_families = set()
    neutral_families = set()
    neutral_names = {fam.name for fam in func_unit.neutral}
    families_names = {fam.name for fam in func_unit.families}
    for gf in gene_families:
        is_neutral = False
        is_functional = False
        for family in gf2fam[gf.name]:
            if family in neutral_names:
                is_neutral = True
            elif family in families_names:
                is_functional = True

        if is_functional:
            fu_families.add(gf)
        elif is_neutral:
            neutral_families.add(gf)
    return fu_families, neutral_families


def search_system_units(model: Model, gene_families: Set[GeneFamily], gf2fam: Dict[str, Set[str]],
                        fam2source: Dict[str, str], source: str, jaccard_threshold: float = 0.8
                        ) -> Dict[str, Set[SystemUnit]]:
    """

    Args:
        model (Model): Model corresponding to the system searched.
        gene_families (Set[GeneFamily]): Set of gene families that might be in the system.
        gf2fam (Dict[str, Set[Family]]): Dictionary linking gene families to their families.
        fam2source (Dict[str, str]): Dictionary linking families to their sources.
        source (str): Name of the annotation source
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.

    Returns:
        Dict[str, Set[SystemUnit]]: System unit found with their name as key and the units as value.
    """
    su_found = {}
    for func_unit in model.func_units:
        fu_families = set()
        fu_families.update(*get_functional_unit_gene_families(func_unit, gene_families, gf2fam))
        if len(fu_families) >= func_unit.min_total:
            matrix, md_gfs2meta, acc_gfs2meta = get_gfs_matrix_combination(fu_families, gf2fam, func_unit, fam2source)
            if check_needed_families(matrix, func_unit):
                context, combinations2orgs = compute_gene_context_graph(families=fu_families,
                                                                        transitive=func_unit.transitivity,
                                                                        window_size=func_unit.window,
                                                                        disable_bar=True)
                combinations2orgs = {combs: org for combs, org in combinations2orgs.items() if
                                     len(combs) >= func_unit.min_total}
                t1 = time.time()
                new_detected_fu = search_unit_in_context(context, fu_families, gf2fam, fam2source, matrix, md_gfs2meta,
                                                         acc_gfs2meta, func_unit, source, combinations2orgs,
                                                         jaccard_threshold)
                if time.time() - t1 > 1:
                    print(f"search unit {func_unit.name} : {time.time() - t1}")
                if len(new_detected_fu) > 0:
                    su_found[func_unit.name] = new_detected_fu
    return su_found


def check_for_needed_units(su_found: Dict[str, Set[SystemUnit]], model: Model) -> bool:
    """
    Checks if the presence/absence rules for needed functional units are respected.

    Args:
        su_found (Dict[str, Set[SystemUnit]]): Dictionary with all system unit found sorted by name.
        model (Model): Model corresponding to the system checked.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    count_mandatory, count_total = (0, 0)
    mandatory_list, accessory_list = (list(map(lambda x: x.name, model.mandatory)),
                                      list(map(lambda x: x.name, model.accessory)))
    for name, fu_set in su_found.items():
        if name in mandatory_list:
            if len(fu_set) > 0:
                count_mandatory += 1
                count_total += 1
        elif name in accessory_list and len(fu_set) > 0:
            count_total += 1

    if (count_mandatory >= model.min_mandatory or model.min_mandatory == -1) and \
            (count_total >= model.min_total or model.min_total == -1):
        return True
    else:
        return False


def check_for_forbidden_unit(su_found: Dict[str, Set[SystemUnit]], model: Model) -> bool:
    """
    Checks if there are forbidden system unit.

    Args:
        su_found (Dict[str, Set[SystemUnit]]): Dictionary with all system unit found sorted by name.
        model (Model): Model corresponding to the system checked.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    forbidden_list = list(map(lambda x: x.name, model.forbidden))

    for name, fu_set in su_found.items():
        if name in forbidden_list and len(fu_set) > 0:
            return True
    return False


def get_system_unit_combinations(su_found: Dict[str, Set[SystemUnit]], model: Model):
    """
       Generate combinations of system units that could code for a system.

       The function generates combinations from mandatory, optional, and neutral categories,
       ensuring that model parameters are respected. Combinations with and without elements
       from neutral categories are generated.

       Args:
           su_found (Dict[str, Set[SystemUnit]]): A dictionary where keys are functional unit names and
                        values are sets of elements belonging to each functional unit model.
           model (Model): Model corresponding to the functional unit.

       Returns:
           list: A list of all possible valid combinations based on the model.
       """
    mandatory_list = list(map(lambda x: x.name, model.mandatory))
    accessory_list = list(map(lambda x: x.name, model.accessory))
    neutral_list = list(map(lambda x: x.name, model.neutral))
    # Extract elements by category
    mandatory_unit = {name: su_found[name] for name in mandatory_list if name in su_found}
    accessory_unit = {name: su_found[name] for name in accessory_list if name in su_found}
    neutral_unit = {name: su_found[name] for name in neutral_list if name in su_found}

    # Initialize a list to store valid combinations
    valid_combinations = []

    # Generate combinations with at least min_mandatory mandatory categories
    for num_mandatory in range(model.min_mandatory, len(mandatory_unit) + 1):
        mandatory_combos = itertools.combinations(mandatory_unit, num_mandatory)

        for mandatory_combo in mandatory_combos:
            # Combine with accessory categories to meet the min_total constraint
            for num_accessory in range(max(0, model.min_total - num_mandatory), len(accessory_unit) + 1):
                accessory_combos = itertools.combinations(accessory_unit, num_accessory)

                for accessory_combo in accessory_combos:
                    # Form a combination of selected mandatory and accessory categories
                    combo_min_acc = list(mandatory_combo) + list(accessory_combo)

                    # Select one element from each category
                    combo_units = []
                    for cat in combo_min_acc:
                        combo_units.append(list(su_found[cat]))

                    # Generate cartesian product of elements (one element per category)
                    product = list(itertools.product(*combo_units))

                    for p in product:
                        # Store the combination without neutral elements
                        final_combo = list(p)
                        valid_combinations.append(final_combo)

                        # Generate and store the combination with neutral elements
                        final_combo_with_neutral = list(final_combo)
                        for neutral_cat in neutral_unit:
                            final_combo_with_neutral.append(
                                list(neutral_unit[neutral_cat])[0])  # Add one neutral element
                        valid_combinations.append(final_combo_with_neutral)

    return valid_combinations


def search_for_system(model: Model, su_found: Dict[str, Set[SystemUnit]], source: str,
                      jaccard_threshold: float = 0.8) -> Set[System]:
    """
    Search a system corresponding to the model based on the unit found.

    Args:
        model (Model): Model corresponding to the system searched.
        su_found (Dict[str, Set[SystemUnit]]): the system unit found for the model.
        source (str): Name of the annotation source
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.


    Returns:
        Set[System]: Systems detected
    """
    gene_families = {fam for su_set in su_found.values() for su in su_set for fam in su.gene_families}
    context, _ = compute_gene_context_graph(families=gene_families, transitive=model.transitivity,
                                            window_size=model.window, disable_bar=True)
    detect_system = set()
    for combination in get_system_unit_combinations(su_found, model):
        fam2su = {}
        contracted_graph = context.copy()
        organisms = set()
        for su in combination:
            fam_list = list(su.gene_families)
            u = fam_list.pop(0)
            for v in fam_list:
                nx.contracted_nodes(contracted_graph, u, v)
            fam2su[u] = su
            organisms |= set(su.organisms)
        filter_local_context(contracted_graph, organisms, jaccard_threshold)
        for cc in sorted(nx.connected_components(contracted_graph), key=len, reverse=True):
            cc: Set[GeneFamily]
            su_in_cc = {fam2su[fam].name: fam2su[fam] for fam in cc if fam in fam2su}
            if check_for_needed_units(su_in_cc, model):
                detect_system.add(System(model, source, units=set(su_in_cc.values())))
    return detect_system


def search_system(model: Model, gene_families: Set[GeneFamily], gf2fam: Dict[str, Set[str]],
                  fam2source: Dict[str, str], source: str, jaccard_threshold: float = 0.8) -> Set[System]:
    """
    Searches for a model system in a pangenome.

    Args:
        model (Model): Model to search in the pangenome.
        meta2fam (dict): Dictionary with each metadata source mapping to a dictionary of metadata linked to sets of gene families.
        source (str): Name of the annotation source.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.

    Returns:
        set: Set of systems detected in the pangenome for the given model.
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for model {model.name}")
    begin = time.time()
    # gene_families, gf2fam, fam2source = dict_families_context(model, meta2fam)
    su_found = search_system_units(model, gene_families, gf2fam, fam2source, source, jaccard_threshold)
    detected_systems = set()
    if check_for_needed_units(su_found, model) and not check_for_forbidden_unit(su_found, model):
        if len(su_found) == 1:  # only one functional unit found so each one correspond to a system
            for fu in next(iter(su_found.values())):
                new_sys = System(model=model, source=source)
                new_sys.add_unit(fu)
                detected_systems.add(new_sys)
        else:
            detected_systems = search_for_system(model, su_found, source, jaccard_threshold)
    logging.getLogger("PANORAMA").debug(f"Done search for model {model.name} in {time.time() - begin} seconds")
    for sys in detected_systems:
        state = pickle.dumps(sys)
        new_sys = pickle.loads(state)
        # for gf in sys.gene_families:
        #     gf._genes_getter = {}
        #     gf._genePerOrg = {}
    gf = list(gene_families)[0]
    return gf


def mp_search_system(args):
    return search_system(*args)


def search_systems(models: Models, pangenome: Pangenome, source: str, metadata_sources: List[str],
                   jaccard_threshold: float = 0.8, threads: int = 1, lock: Lock = None,
                   disable_bar: bool = False):
    """
    Searches for systems present in the pangenome for all models.

    Args:
        models (Models): Models to search in pangenomes.
        pangenome (Pangenome): Pangenome object containing gene families.
        source (str): Name of the source for the system.
        metadata_sources (list of str): List of the metadata sources for the families.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
        threads (int, optional): Number of available threads. Default is 1.
        lock (Lock, optional): Global lock for multiprocessing execution. Default is None.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.
    """
    logging.getLogger("PANORAMA").debug(f"Begin systems detection with {threads} threads in {pangenome.name}")
    begin = time.time()

    families2gfs = get_metadata_to_families(pangenome=pangenome, sources=metadata_sources)
    detected_systems = set()

    with ProcessPoolExecutor(max_workers=threads, mp_context=get_context("fork"),
                             initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=models.size, unit='model', disable=disable_bar) as progress:
            futures = []
            for model in models:
                gene_families, gf2fam, fam2source = dict_families_context(model, families2gfs)
                future = executor.submit(search_system, model, gene_families, gf2fam, fam2source, source, jaccard_threshold)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

                for future in futures:
                    result = future.result()
                    detected_systems |= result

    # args_search = []
    # for model in models:
    #     gene_families, gf2fam, fam2source = dict_families_context(model, families2gfs)
    #     args_search.append((model, gene_families, gf2fam, fam2source, source, jaccard_threshold))
    # with get_context("fork").Pool(processes=threads) as pool:
    #     for res in pool.imap(mp_search_system, args_search):
    #         detected_systems |= res

    for idx, system in enumerate(sorted(detected_systems,
                                        key=lambda x: (len(x.model.canonical), -len(x), -x.number_of_families)),
                                 start=1):
        system.ID = str(idx)
        pangenome.add_system(system)

    logging.getLogger("PANORAMA").debug(f"{pangenome.number_of_systems(source)} systems detected in {pangenome.name}")
    logging.getLogger("PANORAMA").info(f"Systems prediction done for {pangenome.name} in {time.time() - begin} seconds")


def search_systems_in_pangenomes(models: Models, pangenomes: Pangenomes, source: str, metadata_sources: List[str],
                                 jaccard_threshold: float = 0.8, threads: int = 1,
                                 lock: Lock = None, disable_bar: bool = False):
    """
    Searches for systems in pangenomes by multithreading on pangenomes.

    Args:
        models (Models): Models to search in pangenomes.
        pangenomes (Pangenomes): Getter object with Pangenome.
        source (str): Name of the source for the system.
        metadata_sources (list of str): List of the metadata sources for the families.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
        threads (int, optional): Number of available threads. Default is 1.
        lock (Lock, optional): Global lock for multiprocessing execution. Default is None.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.
    """
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Begin systems searching for {pangenome.name}")
        search_systems(models, pangenome, source, metadata_sources, jaccard_threshold,
                       threads, lock, disable_bar)


def write_systems_to_pangenome(pangenome: Pangenome, source: str, disable_bar: bool = False):
    """
    Writes detected systems to the pangenome.

    Args:
        pangenome (Pangenome): Pangenome object containing detected systems.
        source (str): Name of the annotation source.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.
    """
    if pangenome.number_of_systems(source) > 0:
        pangenome.status["systems"] = "Computed"
        logging.getLogger("PANORAMA").info(f"Write systems in pangenome {pangenome.name}")
        write_pangenome(pangenome, pangenome.file, source=source, disable_bar=disable_bar)
        logging.getLogger("PANORAMA").info(f"Systems written in pangenome {pangenome.name}")
    else:
        logging.getLogger("PANORAMA").info("No system detected")


def write_systems_to_pangenomes(pangenomes: Pangenomes, source: str, threads: int = 1, lock: Lock = None,
                                disable_bar: bool = False):
    """
    Writes detected systems into pangenomes.

    Args:
        pangenomes (Pangenomes): Pangenomes object containing all the pangenomes with systems.
        source (str): Metadata source.
        threads (int, optional): Number of available threads. Default is 1.
        lock (Lock, optional): Lock for multiprocessing execution. Default is None.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.
    """
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome in tqdm(pangenomes, unit='pangenome', disable=disable_bar):
                logging.getLogger("PANORAMA").debug(f"Write systems for pangenome {pangenome.name}")
                future = executor.submit(write_systems_to_pangenome, pangenome, source, disable_bar)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def launch(args):
    """
    Launches functions to detect systems in pangenomes.

    Args:
        args (Namespace): Argument given in CLI.
    """
    from panorama.utility.utility import check_models
    from panorama.format.read_binaries import load_pangenomes

    need_info = check_detection_args(args)
    models = check_models(args.models, disable_bar=args.disable_prog_bar)
    manager = Manager()
    lock = manager.Lock()
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_detection, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, metadata_sources=args.annotation_sources,
                                 systems_source=args.source, force=args.force)
    search_systems_in_pangenomes(models=models, pangenomes=pangenomes, source=args.source,
                                 metadata_sources=args.annotation_sources, jaccard_threshold=args.jaccard,
                                 threads=args.threads, lock=lock, disable_bar=args.disable_prog_bar)
    write_systems_to_pangenomes(pangenomes, args.source, args.threads, lock, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Creates a subparser to launch PANORAMA from the command line.

    Args:
        sub_parser (argparse.ArgumentParser): Sub-parser for the 'systems' command.

    Returns:
        argparse.ArgumentParser: Parser with arguments for the 'systems' command.
    """
    parser = sub_parser.add_parser("systems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_detection(parser)
    return parser


def parser_detection(parser):
    """
    Adds arguments to the parser for the 'systems' command.

    Args:
        parser (argparse.ArgumentParser): Parser for the 'systems' command.

    TODO:
        - Add an option to write projection.
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required:")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in a .tsv file.')
    required.add_argument('-m', '--models', required=True, type=Path, nargs='?',
                          help="Path to model list file. Note: Use 'panorama utils --models' to create the models list file.")
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source to select in pangenomes.')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--annotation_sources', required=False, type=str, default=None, nargs='+',
                          help="Name of the annotation sources to load if different from the system source. "
                               "Can specify more than one, separated by space.")
    optional.add_argument('--jaccard', required=False, type=float, default=0.8,
                          help="Minimum Jaccard similarity used to filter edges between gene families. "
                               "Increasing this value improves precision but significantly lowers sensitivity.")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads.")
