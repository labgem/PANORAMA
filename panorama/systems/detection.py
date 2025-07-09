#!/usr/bin/env python3
# coding:utf-8

"""This module provides functions to detect biological systems in pangenomes."""

# Default libraries
from __future__ import annotations
import argparse
import itertools
import time
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Set, Tuple, FrozenSet, Union
from multiprocessing import Manager, Lock

# Installed libraries
import pandas as pd
import networkx as nx
from tqdm import tqdm
from ppanggolin.genome import Organism
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import compute_gene_context_graph
from ppanggolin.metadata import Metadata

# Local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.utils import (
    filter_local_context,
    filter_global_context,
    check_for_families,
    get_metadata_to_families,
    get_gfs_matrix_combination,
    dict_families_context,
    check_needed_families,
)
from panorama.systems.models import Models, Model, FuncUnit, Family
from panorama.systems.system import System, SystemUnit
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.utils import init_lock
from panorama.format.write_binaries import write_pangenome, erase_pangenome


def check_detection_args(
    args: argparse.Namespace,
) -> Dict[str, Union[bool, str, List[str]]]:
    """Checks and processes the provided arguments to ensure they are valid.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
         Dict[str, Union[bool, str, List[str]]]: A dictionary indicating necessary information to load the pangenome.

    Raises:
        argparse.ArgumentTypeError: If 'jaccard' is not a restricted float.
    """
    args.jaccard = restricted_float(args.jaccard)
    args.annotation_sources = (
        [args.source] if args.annotation_sources is None else args.annotation_sources
    )
    return {
        "need_annotations": True,
        "need_families": True,
        "need_metadata": True,
        "metatypes": ["families"],
        "sources": args.annotation_sources,
    }


def check_pangenome_detection(
    pangenome: Pangenome,
    metadata_sources: List[str],
    systems_source: str,
    force: bool = False,
):
    """Checks and loads pangenome information before adding families.

    Args:
        pangenome (Pangenome): Pangenome object to be checked.
        metadata_sources (List[str]): Sources used to associate families with gene families.
        systems_source (str): Source used to detect systems.
        force (bool, optional): If True, forces the erasure of pangenome systems from the source. Defaults to False.

    Raises:
        KeyError: If the provided annotation source is not in the pangenome.
        ValueError: If systems are already detected based on the source and 'force' is not used.
        AttributeError: If there is no metadata associated with families.
    """
    if (
        pangenome.status["systems"] == "inFile"
        and systems_source in pangenome.status["systems_sources"]
    ):
        if force:
            erase_pangenome(pangenome, systems=True, source=systems_source)
        else:
            raise ValueError(
                f"Systems are already detected based on the source: {systems_source}. "
                f"Use the --force option to erase the already computed systems."
            )
    if pangenome.status["metadata"]["families"] == "inFile":
        for annotation_source in metadata_sources:
            if annotation_source not in pangenome.status["metasources"]["families"]:
                raise KeyError(
                    f"There is no metadata associated with families "
                    f"from source {annotation_source} in pangenome {pangenome.name}."
                )
    elif pangenome.status["metadata"]["families"] not in ["Computed", "Loaded"]:
        raise AttributeError(
            f"There is no metadata associated with families in your pangenome {pangenome.name}. "
            "Please see the command annotation before detecting systems"
        )


def get_subcombinations(
    combi: Set[GeneFamily], combinations: List[FrozenSet[GeneFamily]]
) -> List[FrozenSet[GeneFamily]]:
    """Removes a combination and all its sub-combinations from a list of combinations.

    Args:
        combi (Set[GeneFamily]): Combination to be removed.
        combinations (List[FrozenSet[GeneFamily]]): List of combinations to be filtered.

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


def search_unit_in_cc(
    graph: nx.Graph,
    system_gf: Set[GeneFamily],
    func_unit: FuncUnit,
    source: str,
    gf2fam: Dict[GeneFamily, Set[Family]],
    fam2source: Dict[str, str],
    detected_su: Set[SystemUnit],
) -> Set[FrozenSet[GeneFamily]]:
    """Searches for functional units in connected components of a graph.

    Args:
        graph (nx.Graph): The graph to search within.
        system_gf (Set[GeneFamily]): Set of gene families in the system.
        func_unit (FuncUnit): Functional unit to search for.
        source (str): Source of the data.
        gf2fam (Dict[GeneFamily, Set[Family]]): Dictionary linking gene families to their families.
        fam2source (Dict[str, str]): Dictionary linking families to their sources.
        detected_su (Set[SystemUnit]): Set of already detected system units.

    Returns:
        Set[FrozenSet[GeneFamily]]: Set of found combinations of gene families.
    """
    found_combinations = set()
    for cc in sorted(
        [cc for cc in nx.connected_components(graph) if len(cc) >= func_unit.min_total],
        key=lambda c: (len(c), sorted(c)),
        reverse=True,
    ):
        cc: Set[GeneFamily]
        gfs_in_cc = system_gf.intersection(cc)
        check, gf2meta_info = check_for_families(
            gfs_in_cc, gf2fam, fam2source, func_unit
        )
        if check:
            new_unit = SystemUnit(
                functional_unit=func_unit,
                source=source,
                gene_families=cc,
                families_to_metainfo=gf2meta_info,
            )
            is_subset = False
            for unit in sorted(detected_su, key=len, reverse=True):
                if new_unit.is_subset(unit):
                    unit.merge(new_unit)
                    is_subset = True
            if not is_subset:
                detected_su.add(new_unit)
            found_combinations.add(frozenset(gfs_in_cc))
    return found_combinations


def search_unit_in_combination(
    graph: nx.Graph,
    families: Set[GeneFamily],
    gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
    mod_fam2meta_source: Dict[str, str],
    func_unit: FuncUnit,
    source: str,
    matrix: pd.DataFrame,
    combinations: List[FrozenSet[GeneFamily]],
    combinations2orgs: Dict[FrozenSet[GeneFamily], Set[Organism]],
    jaccard_threshold: float = 0.8,
    local: bool = False,
) -> Set[SystemUnit]:
    """Searches for functional unit corresponding to a model in a graph.

    Args:
        graph (nx.Graph): A graph with families in one connected component.
        families (Set[GeneFamily]): A set of families that code for the searched model.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]]): Dictionary linking gene families to their families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking model families to metadata sources.
        func_unit (FuncUnit): One functional unit corresponding to the model.
        source (str): Name of the source.
        matrix (pd.DataFrame): Dataframe containing association between gene families and unit families.
        combinations (List[FrozenSet[GeneFamily]]): Existing combination of gene families in organisms.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        combinations2orgs (Dict[FrozenSet[GeneFamily], Set[Organism]]): The combination of gene families corresponding to the context that exists in at least one genome.
        local: (bool, optional): Whether to filter the context with a local Jaccard index or not. Defaults to False.

    Returns:
        Set[SystemUnit]: Set of all detected functional units in the graph.
    """
    detected_su: Set[SystemUnit] = set()
    mandatory_gfs = {
        gf
        for gf in families
        for fam in gene_fam2mod_fam[gf]
        if fam.presence == "mandatory"
    }
    accessory_gfs = {
        gf
        for gf in families
        for fam in gene_fam2mod_fam[gf]
        if fam.presence == "accessory"
    }

    while len(combinations) > 0:
        logging.getLogger("PANORAMA").debug(f"Search in {len(combinations)}")
        context_combination = combinations.pop(0)
        organisms_of_interest = combinations2orgs[context_combination]
        context_combination = set(context_combination)
        context_gfs_name = {gf.name for gf in context_combination} & {
            gf.name for gf in mandatory_gfs | accessory_gfs
        }

        if len(context_gfs_name) >= func_unit.min_total:
            filtered_matrix = matrix[list(context_gfs_name)]
            logging.getLogger("PANORAMA").debug("Check needed")
            if check_needed_families(filtered_matrix, func_unit):
                local_graph = graph.copy()
                node2remove = families.difference(context_combination)
                local_graph.remove_nodes_from(node2remove)
                logging.getLogger("PANORAMA").debug(
                    f"filter context with jaccard {jaccard_threshold}"
                )
                if local:
                    local_graph = filter_local_context(
                        local_graph,
                        organisms_of_interest,
                        jaccard_threshold=jaccard_threshold,
                    )
                logging.getLogger("PANORAMA").debug("search unit in filtered context")
                founds_combs = search_unit_in_cc(
                    local_graph,
                    context_combination,
                    func_unit,
                    source,
                    gene_fam2mod_fam,
                    mod_fam2meta_source,
                    detected_su,
                )
                for comb in founds_combs:
                    get_subcombinations(set(comb), combinations)
            else:
                get_subcombinations(context_combination, combinations)
        else:
            get_subcombinations(context_combination, combinations)
    return detected_su


def search_unit_in_context(
    graph: nx.Graph,
    families: Set[GeneFamily],
    gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
    mod_fam2meta_source: Dict[str, str],
    matrix: pd.DataFrame,
    func_unit: FuncUnit,
    source: str,
    combinations2orgs: Dict[FrozenSet[GeneFamily], Set[Organism]] = None,
    jaccard_threshold: float = 0.8,
    local: bool = False,
) -> Set[SystemUnit]:
    """Searches for system unit corresponding to a model in a pangenomic context.

    Args:
        graph (nx.Graph): A graph with families in a pangenomic context.
        families (Set[GeneFamily]): A set of families that code for the searched model.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]]): Dictionary linking gene families to model families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking model families to metadata sources.
        func_unit (FuncUnit): One functional unit corresponding to the model.
        source (str): Name of the source.
        matrix (pd.DataFrame): Dataframe containing association between gene families and unit families.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        combinations2orgs (Dict[FrozenSet[GeneFamily], Set[Organism]], optional): The combination of gene families corresponding to the context that exists in at least one genome. Defaults to None.
        local: (bool, optional): Whether to filter the context with a local Jaccard index or not. Defaults to False.

    Returns:
       Set[SystemUnit]: Set of detected systems in the pangenomic context.
    """
    detected_su: Set[SystemUnit] = set()
    families_in_context, neutral_families_in_context = (
        get_functional_unit_gene_families(func_unit, families, gene_fam2mod_fam)
    )
    combinations = sorted(
        combinations2orgs.keys(), key=lambda fs: (len(fs), sorted(fs)), reverse=True
    )
    if not local:
        graph = filter_global_context(graph, jaccard_threshold)

    for cc_graph in [
        graph.subgraph(cc).copy()
        for cc in sorted(nx.connected_components(graph), key=len, reverse=True)
    ]:
        families_in_cc = families_in_context.copy()
        neutral_families = neutral_families_in_context.copy()
        families_in_cc &= cc_graph.nodes
        neutral_families &= cc_graph.nodes
        if len(families_in_cc) > 0:
            families_in_cc |= neutral_families
            combinations_in_cc = get_subcombinations(families_in_cc, combinations)
            if len(combinations_in_cc) > 0:
                logging.getLogger("PANORAMA").debug("Search unit in cc")
                t0 = time.time()
                # TODO filter is never done
                new_detected_su = search_unit_in_combination(
                    cc_graph,
                    families_in_cc,
                    gene_fam2mod_fam,
                    mod_fam2meta_source,
                    func_unit,
                    source,
                    matrix,
                    combinations_in_cc,
                    combinations2orgs,
                    jaccard_threshold,
                )
                logging.getLogger("PANORAMA").debug(
                    f"{len(detected_su)} unit found in cc in {time.time() - t0} seconds."
                )
                logging.getLogger("PANORAMA").debug("filter unit in cc")
                t1 = time.time()
                for new_su in new_detected_su:
                    is_subset = False
                    for unit in sorted(detected_su, key=len, reverse=True):
                        if new_su.is_subset(unit):
                            is_subset = True
                            unit.merge(new_su)
                    if not is_subset:
                        detected_su.add(new_su)
                logging.getLogger("PANORAMA").debug(
                    f"{len(detected_su)} unit found after filtering in {time.time() - t1} seconds."
                )
    return detected_su


def get_functional_unit_gene_families(
    func_unit: FuncUnit,
    gene_families: Set[GeneFamily],
    gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
) -> Tuple[Set[GeneFamily], Set[GeneFamily]]:
    """Retrieves the gene families that might be in the functional unit.

    Args:
        func_unit (FuncUnit): The functional unit to consider.
        gene_families (Set[GeneFamily]): Set of gene families that might be in the system.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]]): Dictionary linking gene families to model families.

    Returns:
        Tuple[Set[GeneFamily], Set[GeneFamily]]:
            - gene families that code for the functional unit
            - neutral gene families of the function unit.
    """
    fu_families = set()
    neutral_families = set()
    for gf in gene_families:
        is_neutral = False
        is_functional = False
        for family in gene_fam2mod_fam[gf]:
            if family in func_unit.neutral:
                is_neutral = True
            elif family in func_unit.families:
                is_functional = True
        if is_functional:
            fu_families.add(gf)
        elif is_neutral:
            neutral_families.add(gf)
    return fu_families, neutral_families


def search_system_units(
    model: Model,
    gene_families: Set[GeneFamily],
    gf2fam: Dict[GeneFamily, Set[Family]],
    fam2source: Dict[str, str],
    source: str,
    jaccard_threshold: float = 0.8,
    sensitivity: int = 1,
) -> Dict[str, Set[SystemUnit]]:
    """Searches for system units corresponding to a model.

    Args:
        model (Model): Model corresponding to the system searched.
        gene_families (Set[GeneFamily]): Set of gene families that might be in the system.
        gf2fam (Dict[GeneFamily, Set[Family]]): Dictionary linking gene families to their families.
        fam2source (Dict[str, str]): Dictionary linking families to their sources.
        source (str): Name of the annotation source.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        sensitivity (int, optional): Sensitivity level for detection.
            1 corresponds to a global Jaccard filtering on the context without looking at all the combinations.
            2 corresponds to a global Jaccard filtering on the specific context of each combination.
            3 corresponds to a local Jaccard filtering on the specific context of each combination.
            Defaults to 1.

    Returns:
        Dict[str, Set[SystemUnit]]: System units found with their name as the key and units as value.

    Raises:
        ValueError: If sensitivity is not 1, 2, or 3.
    """
    su_found = {}
    for func_unit in model.func_units:
        fu_families = set()
        fu_families.update(
            *get_functional_unit_gene_families(func_unit, gene_families, gf2fam)
        )
        if len(fu_families) >= func_unit.min_total:
            matrix = get_gfs_matrix_combination(
                fu_families, gf2fam
            )
            if check_needed_families(matrix, func_unit):
                logging.getLogger("PANORAMA").debug("Extract Genomic context")
                t0 = time.time()
                context, combinations2orgs = compute_gene_context_graph(
                    families=fu_families,
                    transitive=func_unit.transitivity,
                    window_size=func_unit.window,
                    disable_bar=True,
                )
                detected_su: Set[SystemUnit] = set()
                if sensitivity == 1:
                    filtered_context = filter_global_context(context, jaccard_threshold)
                    search_unit_in_cc(
                        filtered_context,
                        fu_families,
                        func_unit,
                        source,
                        gf2fam,
                        fam2source,
                        detected_su,
                    )
                elif sensitivity == 2 or sensitivity == 3:
                    combinations2orgs = {
                        combs: org
                        for combs, org in combinations2orgs.items()
                        if len(combs) >= func_unit.min_total
                    }
                    logging.getLogger("PANORAMA").debug(
                        f"Genomic context extracted in {time.time() - t0} seconds. "
                        f"{len(combinations2orgs)} combinations found"
                    )
                    t1 = time.time()
                    logging.getLogger("PANORAMA").debug("Search unit in context")
                    detected_su = search_unit_in_context(
                        context,
                        fu_families,
                        gf2fam,
                        fam2source,
                        matrix,
                        func_unit,
                        source,
                        combinations2orgs,
                        jaccard_threshold,
                        local=True if sensitivity == 3 else False,
                    )
                    logging.getLogger("PANORAMA").debug(
                        f"{len(detected_su)} unit found in {time.time() - t1} seconds."
                    )
                else:
                    raise ValueError(
                        "Sensitivity is expected to be either 1 or 2 or 3."
                    )
                if len(detected_su) > 0:
                    su_found[func_unit.name] = detected_su
    return su_found


def check_for_needed_units(su_found: Dict[str, Set[SystemUnit]], model: Model) -> bool:
    """Checks if the presence/absence rules for necessary functional units are respected.

    Args:
        su_found (Dict[str, Set[SystemUnit]]): Dictionary with all system units found sorted by name.
        model (Model): Model corresponding to the system checked.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    count_mandatory, count_total = (0, 0)
    mandatory_list, accessory_list = (
        list(map(lambda x: x.name, model.mandatory)),
        list(map(lambda x: x.name, model.accessory)),
    )
    for name, fu_set in su_found.items():
        if name in mandatory_list:
            if len(fu_set) > 0:
                count_mandatory += 1
                count_total += 1
        elif name in accessory_list and len(fu_set) > 0:
            count_total += 1
    return (count_mandatory >= model.min_mandatory or model.min_mandatory == -1) and (
        count_total >= model.min_total or model.min_total == -1
    )


def check_for_forbidden_unit(
    su_found: Dict[str, Set[SystemUnit]], model: Model
) -> bool:
    """Checks if there are forbidden system units.

    Args:
        su_found (Dict[str, Set[SystemUnit]]): Dictionary with all system units found sorted by name.
        model (Model): Model corresponding to the system checked.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    forbidden_list = list(map(lambda x: x.name, model.forbidden))
    for name, fu_set in su_found.items():
        if name in forbidden_list and len(fu_set) > 0:
            return True
    return False


def get_system_unit_combinations(
    su_found: Dict[str, Set[SystemUnit]], model: Model
) -> List[List[SystemUnit]]:
    """Generates combinations of system units that could code for a system.

    The function generates combinations from mandatory, optional, and neutral categories,
    ensuring that model parameters are respected. Combinations with and without elements
    from neutral categories are generated.

    Args:
        su_found (Dict[str, Set[SystemUnit]]): A dictionary where keys are functional unit names and values are sets of
                 elements belonging to each functional unit model.
        model (Model): Model corresponding to the functional unit.

    Returns:
        List[List[SystemUnit]]: A list of all possible valid combinations based on the model.
    """
    mandatory_list = list(map(lambda x: x.name, model.mandatory))
    accessory_list = list(map(lambda x: x.name, model.accessory))
    neutral_list = list(map(lambda x: x.name, model.neutral))

    mandatory_unit = {
        name: su_found[name] for name in mandatory_list if name in su_found
    }
    accessory_unit = {
        name: su_found[name] for name in accessory_list if name in su_found
    }
    neutral_unit = {name: su_found[name] for name in neutral_list if name in su_found}

    valid_combinations = []

    for num_mandatory in range(model.min_mandatory, len(mandatory_unit) + 1):
        mandatory_combos = itertools.combinations(mandatory_unit, num_mandatory)
        for mandatory_combo in mandatory_combos:
            for num_accessory in range(
                max(0, model.min_total - num_mandatory), len(accessory_unit) + 1
            ):
                accessory_combos = itertools.combinations(accessory_unit, num_accessory)
                for accessory_combo in accessory_combos:
                    combo_min_acc = list(mandatory_combo) + list(accessory_combo)
                    combo_units = []
                    for cat in combo_min_acc:
                        combo_units.append(list(su_found[cat]))
                    product = list(itertools.product(*combo_units))
                    for p in product:
                        final_combo = list(p)
                        valid_combinations.append(final_combo)
                        final_combo_with_neutral = list(final_combo)
                        for neutral_cat in neutral_unit:
                            final_combo_with_neutral.append(
                                list(neutral_unit[neutral_cat])[0]
                            )
                        valid_combinations.append(final_combo_with_neutral)
    return valid_combinations


def search_for_system(
    model: Model,
    su_found: Dict[str, Set[SystemUnit]],
    source: str,
    jaccard_threshold: float = 0.8,
) -> Set[System]:
    """Searches for a system corresponding to the model based on the units found.

    Args:
        model (Model): Model corresponding to the system searched.
        su_found (Dict[str, Set[SystemUnit]]): The system units found for the model.
        source (str): Name of the annotation source.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.

    Returns:
        set of System: Systems detected.
    """
    gene_families = {
        fam for su_set in su_found.values() for su in su_set for fam in su.families
    }
    context, _ = compute_gene_context_graph(
        families=gene_families,
        transitive=model.transitivity,
        window_size=model.window,
        disable_bar=True,
    )
    detect_system = set()
    for combination in get_system_unit_combinations(su_found, model):
        fam2su = {}
        contracted_graph = context.copy()
        organisms = set()
        for su in combination:
            fam_list = list(su.families)
            u = fam_list.pop(0)
            for v in fam_list:
                nx.contracted_nodes(contracted_graph, u, v)
            fam2su[u] = su
            organisms |= set(su.organisms)
        filter_local_context(contracted_graph, organisms, jaccard_threshold)
        for cc in sorted(
            nx.connected_components(contracted_graph), key=len, reverse=True
        ):
            cc: Set[GeneFamily]
            su_in_cc = {fam2su[fam].name: fam2su[fam] for fam in cc if fam in fam2su}
            if check_for_needed_units(su_in_cc, model):
                detect_system.add(System(model, source, units=set(su_in_cc.values())))
    return detect_system


def search_system(
    model: Model,
    gene_families: Set[GeneFamily],
    gf2fam: Dict[GeneFamily, Set[Family]],
    fam2source: Dict[str, str],
    source: str,
    jaccard_threshold: float = 0.8,
    sensitivity: int = 1,
) -> Set[System]:
    """Searches for a model system in a pangenome.

    Args:
        model (Model): Model to search in the pangenome.
        gene_families (Set[GeneFamily]): Set of gene families that might be in the system.
        gf2fam (Dict[GeneFamily, Set[Family]]): Dictionary linking gene families to their families.
        fam2source (Dict[str, str]): Dictionary linking families to their sources.
        source (str): Name of the annotation source.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        sensitivity (int, optional): Sensitivity level for detection.
            1 corresponds to a global Jaccard filtering on the context without looking at all the combinations.
            2 corresponds to a global Jaccard filtering on the specific context of each combination.
            3 corresponds to a local Jaccard filtering on the specific context of each combination.
            Defaults to 1.

    Returns:
        set of System: Set of systems detected in the pangenome for the given model.
    """
    logging.getLogger("PANORAMA").debug(f"Begin search for model {model.name}")
    begin = time.time()
    su_found = search_system_units(
        model, gene_families, gf2fam, fam2source, source, jaccard_threshold, sensitivity
    )
    detected_systems = set()
    if check_for_needed_units(su_found, model) and not check_for_forbidden_unit(
        su_found, model
    ):
        if len(su_found) == 1:
            for fu in next(iter(su_found.values())):
                new_sys = System(model=model, source=source)
                new_sys.add_unit(fu)
                detected_systems.add(new_sys)
        else:
            detected_systems = search_for_system(
                model, su_found, source, jaccard_threshold
            )
    logging.getLogger("PANORAMA").debug(
        f"Done search for model {model.name} in {time.time() - begin} seconds"
    )
    return detected_systems


def search_systems(
    models: Models,
    pangenome: Pangenome,
    source: str,
    metadata_sources: List[str],
    jaccard_threshold: float = 0.8,
    sensitivity: int = 1,
    threads: int = 1,
    lock: Lock = None,
    disable_bar: bool = False,
):
    """Searches for systems present in the pangenome for all models.

    Args:
        models (Models): Models to search in pangenomes.
        pangenome (Pangenome): Pangenome object containing gene families.
        source (str): Name of the source for the system.
        metadata_sources (List[str]): List of the metadata sources for the families.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        threads (int, optional): Number of available threads. Defaults to 1.
        lock (Lock, optional): Global lock for multiprocessing execution. Defaults to None.
        disable_bar (bool, optional): If True, disables the progress bar. Defaults to False.
        sensitivity (int, optional): Sensitivity level for detection.
            1 corresponds to a global Jaccard filtering on the context without looking at all the combinations.
            2 corresponds to a global Jaccard filtering on the specific context of each combination.
            3 corresponds to a local Jaccard filtering on the specific context of each combination.
            Defaults to 1.
    """
    logging.getLogger("PANORAMA").debug(
        f"Begin systems detection with {threads} threads in {pangenome.name}"
    )
    begin = time.time()
    meta2fam = get_metadata_to_families(pangenome=pangenome, sources=metadata_sources)
    with ThreadPoolExecutor(
        max_workers=threads, initializer=init_lock, initargs=(lock,)
    ) as executor:
        with tqdm(total=models.size, unit="model", disable=disable_bar) as progress:
            futures = []
            for model in models:
                gene_families, gf2fam, fam2source = dict_families_context(
                    model, meta2fam
                )
                future = executor.submit(
                    search_system,
                    model,
                    gene_families,
                    gf2fam,
                    fam2source,
                    source,
                    jaccard_threshold,
                    sensitivity,
                )
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            detected_systems = set()
            for future in futures:
                result = future.result()
                detected_systems |= result
    for idx, system in enumerate(
        sorted(
            detected_systems,
            key=lambda x: (
                len(x.model.canonical),
                -len(x),
                -x.number_of_model_gene_families,
            ),
        ),
        start=1,
    ):
        system.ID = str(idx)
        pangenome.add_system(system)
    logging.getLogger("PANORAMA").debug(
        f"{pangenome.number_of_systems(source)} systems detected in {pangenome.name}"
    )
    logging.getLogger("PANORAMA").info(
        f"Systems prediction done for {pangenome.name} in {time.time() - begin:.2f} seconds"
    )


def search_systems_in_pangenomes(
    models: Models,
    pangenomes: Pangenomes,
    source: str,
    metadata_sources: List[str],
    jaccard_threshold: float = 0.8,
    sensitivity: int = 1,
    threads: int = 1,
    lock: Lock = None,
    disable_bar: bool = False,
):
    """Searches for systems in pangenomes by multithreading on pangenomes.

    Args:
        models (Models): Models to search in pangenomes.
        pangenomes (Pangenomes): Getter object with Pangenome.
        source (str): Name of the source for the system.
        metadata_sources (List[str]): List of the metadata sources for the families.
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        threads (int, optional): Number of available threads. Defaults to 1.
        lock (Lock, optional): Global lock for multiprocessing execution. Defaults to None.
        disable_bar (bool, optional): If True, disables the progress bar. Defaults to False.
        sensitivity (int, optional): Sensitivity level for detection.
            1 corresponds to a global Jaccard filtering on the context without looking at all the combinations.
            2 corresponds to a global Jaccard filtering on the specific context of each combination.
            3 corresponds to a local Jaccard filtering on the specific context of each combination.
            Defaults to 1.
    """
    t0 = time.time()
    for pangenome in tqdm(
        pangenomes, total=len(pangenomes), unit="pangenome", disable=disable_bar
    ):
        logging.getLogger("PANORAMA").debug(
            f"Begin systems searching for {pangenome.name}"
        )
        search_systems(
            models,
            pangenome,
            source,
            metadata_sources,
            jaccard_threshold,
            sensitivity,
            threads,
            lock,
            disable_bar,
        )
    logging.getLogger("PANORAMA").info(
        f"Systems prediction in all pangenomes done in {time.time() - t0:.2f} seconds"
    )


def write_systems_to_pangenome(
    pangenome: Pangenome, source: str, disable_bar: bool = False
):
    """Writes detected systems to the pangenome.

    Args:
        pangenome (Pangenome): Pangenome object containing detected systems.
        source (str): Name of the annotation source.
        disable_bar (bool, optional): If True, disables the progress bar. Defaults to False.
    """
    if pangenome.number_of_systems(source) > 0:
        pangenome.status["systems"] = "Computed"
        logging.getLogger("PANORAMA").info(
            f"Write systems in pangenome {pangenome.name}"
        )
        write_pangenome(
            pangenome, pangenome.file, source=source, disable_bar=disable_bar
        )
        logging.getLogger("PANORAMA").info(
            f"Systems written in pangenome {pangenome.name}"
        )
    else:
        logging.getLogger("PANORAMA").info("No system detected")


def write_systems_to_pangenomes(
    pangenomes: Pangenomes,
    source: str,
    threads: int = 1,
    lock: Lock = None,
    disable_bar: bool = False,
):
    """Writes detected systems into pangenomes.

    Args:
        pangenomes (Pangenomes): Pangenomes object containing all the pangenomes with systems.
        source (str): Metadata source.
        threads (int, optional): Number of available threads. Defaults to 1.
        lock (Lock, optional): Lock for multiprocessing execution. Defaults to None.
        disable_bar (bool, optional): If True, disables the progress bar. Defaults to False.
    """
    with ThreadPoolExecutor(
        max_workers=threads, initializer=init_lock, initargs=(lock,)
    ) as executor:
        with tqdm(
            total=len(pangenomes), unit="pangenome", disable=disable_bar
        ) as progress:
            futures = []
            for pangenome in tqdm(pangenomes, unit="pangenome", disable=disable_bar):
                logging.getLogger("PANORAMA").debug(
                    f"Write systems for pangenome {pangenome.name}"
                )
                future = executor.submit(
                    write_systems_to_pangenome, pangenome, source, disable_bar
                )
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            for future in futures:
                future.result()


def launch(args):
    """Launches functions to detect systems in pangenomes.

    Args:
        args (argparse.Namespace): Argument given in CLI.
    """
    from panorama.utility.utility import check_models
    from panorama.format.read_binaries import load_pangenomes

    need_info = check_detection_args(args)
    models = check_models(args.models, disable_bar=args.disable_prog_bar)
    manager = Manager()
    lock = manager.Lock()
    pangenomes = load_pangenomes(
        pangenome_list=args.pangenomes,
        need_info=need_info,
        check_function=check_pangenome_detection,
        max_workers=args.threads,
        lock=lock,
        disable_bar=args.disable_prog_bar,
        metadata_sources=args.annotation_sources,
        systems_source=args.source,
        force=args.force,
    )
    search_systems_in_pangenomes(
        models=models,
        pangenomes=pangenomes,
        source=args.source,
        metadata_sources=args.annotation_sources,
        jaccard_threshold=args.jaccard,
        sensitivity=args.sensitivity,
        threads=args.threads,
        lock=lock,
        disable_bar=args.disable_prog_bar,
    )
    write_systems_to_pangenomes(
        pangenomes, args.source, args.threads, lock, disable_bar=args.disable_prog_bar
    )


def subparser(sub_parser) -> argparse.ArgumentParser:
    """Creates a subparser to launch PANORAMA from the command line.

    Args:
        sub_parser (argparse.ArgumentParser): Sub-parser for the 'systems' command.

    Returns:
        argparse.ArgumentParser: Parser with arguments for the 'systems' command.
    """
    parser = sub_parser.add_parser("systems")
    parser_detection(parser)
    return parser


def parser_detection(parser):
    """Adds arguments to the parser for the 'systems' command.

    Args:
        parser (argparse.ArgumentParser): Parser for the 'systems' command.
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="All of the following arguments are required:",
    )
    required.add_argument(
        "-p",
        "--pangenomes",
        required=True,
        type=Path,
        nargs="?",
        help="A list of pangenome .h5 files in a .tsv file.",
    )
    required.add_argument(
        "-m",
        "--models",
        required=True,
        type=Path,
        nargs="?",
        help="Path to model list file. Note: Use 'panorama utils --models' to create the models list file.",
    )
    required.add_argument(
        "-s",
        "--source",
        required=True,
        type=str,
        nargs="?",
        help="Name of the annotation source to select in pangenomes.",
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "--annotation_sources",
        required=False,
        type=str,
        default=None,
        nargs="+",
        help="Name of the annotation sources to load if different from the system source. "
        "Can specify more than one, separated by space.",
    )
    optional.add_argument(
        "--jaccard",
        required=False,
        type=float,
        default=0.8,
        help="Minimum Jaccard similarity used to filter edges between gene families. "
        "Increasing this value improves precision but significantly lowers sensitivity.",
    )
    optional.add_argument(
        "--sensitivity",
        required=False,
        type=int,
        default=3,
        choices=[1, 2, 3],
        help=argparse.SUPPRESS,
    )
    optional.add_argument(
        "--threads",
        required=False,
        nargs="?",
        type=int,
        default=1,
        help="Number of available threads.",
    )
