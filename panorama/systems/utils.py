#!/usr/bin/env python3

"""
This module provides utility functions to detect and write biological systems in pangenomes.
"""

# default libraries
from __future__ import annotations
import logging
from typing import Dict, Iterable, List, Set, Tuple
from collections import defaultdict

# installed libraries
from itertools import combinations
import numpy as np
import pandas as pd
import networkx as nx
from ppanggolin.genome import Organism
from ppanggolin.metadata import Metadata

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family
from panorama.pangenomes import Pangenome

pd.options.mode.copy_on_write = True  # Remove when pandas3.0 available. See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy


def filter_local_context(graph: nx.Graph, organisms: Set[Organism], excluded: set[GeneFamily],
                         jaccard_threshold: float = 0.8) -> nx.Graph:
    """
    Filters a graph based on a local Jaccard index.

    Args:
        graph (nx.Graph): A sub-pangenome graph.
        organisms (Set[Organism]): Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """

    # Compute local Jaccard
    filtered_graph = nx.Graph()
    fam2genome_proportion = {}
    for f1, f2, data in graph.edges(data=True):
        if not (f1 in excluded or f2 in excluded):
            filtered_graph.add_node(f1)
            filtered_graph.add_node(f2)
            if f1 not in fam2genome_proportion:
                f1_orgs = set(f1.organisms).intersection(organisms)
                fam2genome_proportion[f1] = len(data["genomes"]) / len(f1_orgs) if len(f1_orgs) > 0 else 0
            if f2 not in fam2genome_proportion:
                f2_orgs = set(f2.organisms).intersection(organisms)
                fam2genome_proportion[f2] = len(data["genomes"]) / len(f2_orgs) if len(f2_orgs) > 0 else 0

            data['f1'] = f1.name
            data['f2'] = f2.name
            data['f1_jaccard_gene'] = fam2genome_proportion[f1]
            data['f2_jaccard_gene'] = fam2genome_proportion[f2]

            if fam2genome_proportion[f1] >= jaccard_threshold and fam2genome_proportion[f2] >= jaccard_threshold:
                filtered_graph.add_edge(f1, f2)

    return filtered_graph


# def get_number_of_mod_fam(gene_family: GeneFamily, gf2fam) -> int:
#     """Gets the number of model families associated with the gene family.
#
#     Args:
#         gene_family (GeneFamily): The gene family of interest.
#
#     Returns:
#         int: The number of model families associated with the gene family.
#     """
#     model_families = gene_fam2mod_fam.get(gene_family.name)
#     if model_families is not None:
#         return len(model_families)
#     else:
#         return 0


def check_for_forbidden_families(gene_families: Set[GeneFamily], gf2fam: Dict[str, Set[str]],
                                 func_unit: FuncUnit) -> bool:
    """
    Checks if there are forbidden conditions in the families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gf2fam (Dict[str, Set[Family]]): Dictionary linking gene families to model families.
        func_unit (FuncUnit): Functional unit to check against.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    forbidden_list = list(map(lambda x: x.name, func_unit.forbidden))
    name2fam = {fam.name: fam for fam in func_unit.families}
    for node in gene_families:
        if node in gf2fam:
            for name in gf2fam[node.name]:
                family = name2fam[name]
                if family.presence == 'forbidden' and family.name in forbidden_list:  # if node is forbidden
                    return True
    return False


def add_metadata_to_dict(gene_family: GeneFamily, family: Family, presence_gfs2metadata,
                         mod_fam2meta_source: Dict[str, str]):
    """
    Associate gene families to metadata in a dictionary

    Args:
        presence_gfs2metadata: Dictionary to save gene families and their metadata
    """
    avail_name = {family.name}.union(family.exchangeable)
    for meta_id, metadata in gene_family.get_metadata_by_source(mod_fam2meta_source[family.name]).items():
        if metadata.protein_name in avail_name:
            if gene_family in presence_gfs2metadata[family.name]:
                _, current_metadata = presence_gfs2metadata[family.name][gene_family]
                if metadata.score > current_metadata.score:
                    presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)
            else:
                presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)
        elif "secondary_name" in metadata.fields:
            if any(name in avail_name for name in metadata.secondary_name.split(",")):
                if gene_family in presence_gfs2metadata[family.name]:
                    _, current_metadata = presence_gfs2metadata[family.name][gene_family]
                    if metadata.score > current_metadata.score:
                        presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)
                else:
                    presence_gfs2metadata[family.name][gene_family] = (meta_id, metadata)


def get_gfs_matrix_combination(gene_families: Set[GeneFamily], gf2fam: Dict[str, Set[str]],
                               func_unit: FuncUnit, mod_fam2meta_source: Dict[str, str]
                               ) -> Tuple[pd.DataFrame, Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]], Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]]:
    """
    Build a matrix of association between gene families and families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gf2fam (Dict[GeneFamily, Set[Family]): Dictionary linking gene families to model families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking families to metadata sources.

    Returns:
        pd.DataFrame: Matrix of association between gene families and families.
        Dict[GeneFamily, Tuple[int, Metadata]]: Dictionary linking gene families metadata to mandatory families
        Dict[GeneFamily, Tuple[int, Metadata]]: Dictionary linking gene families metadata to accessory families
    """
    mandatory_gfs2metadata = defaultdict(dict)
    accessory_gfs2metadata = defaultdict(dict)
    gfs = set()
    name2fam = {fam.name: fam for fam in func_unit.families}
    for gene_family in gene_families:
        for name in gf2fam[gene_family.name]:
            family = name2fam[name]
            if family.presence == "mandatory":
                add_metadata_to_dict(gene_family, family, mandatory_gfs2metadata, mod_fam2meta_source)
                gfs.add(gene_family)
            elif family.presence == "accessory":
                add_metadata_to_dict(gene_family, family, accessory_gfs2metadata, mod_fam2meta_source)
                gfs.add(gene_family)

    fams = list(mandatory_gfs2metadata.keys()) + list(accessory_gfs2metadata.keys())
    score_matrix = np.zeros((len(fams), len(gfs)))
    gfs = list(gfs)
    for i, fam in enumerate(fams):
        if fam in mandatory_gfs2metadata:
            gfs2metadata = mandatory_gfs2metadata[fam]
        else:
            gfs2metadata = accessory_gfs2metadata[fam]
        for j, gf in enumerate(gfs):
            score_matrix[i, j] = 1 if gf in gfs2metadata else 0
    return (pd.DataFrame(score_matrix, index=fams, columns=[gf.name for gf in gfs]),
            mandatory_gfs2metadata, accessory_gfs2metadata)


def bitset_from_row(row) -> int:
    """
    Converts a row of the binary matrix into an integer bitset.
    Each bit represents the presence of an individual for that characteristic.
    Args:
        row: A row of the binary matrix

    Returns:
        int: integer bitset representing the presence of an individual for that characteristic
    """
    bitset = 0
    for i, val in enumerate(row):
        if val == 1:
            bitset |= (1 << i)  # Met le i-ème bit à 1
    return bitset


def search_comb(comb: Tuple[int, ...], mandatory_bitsets: List[int],
                accessory_bitsets: List[int]) -> Tuple[Dict[int, int], int, int]:
    """
    Search for a working combination in bitsets

    Args:
        comb (Tuple[int, ...]): Combination of gene families index
        mandatory_bitsets (List[int]): Bitsets with index of mandatory families
        accessory_bitsets (List[int]): Bitsets with index of accessory families

    Returns:
        Dict[int, int]: Association between gene families and families
        Integer: Number of mandatory families covered
        Integer: Number of accessory families covered
    """
    covered_families = set()
    gf2fam = {}

    covered_mandatory = 0
    for family_index, bitset in enumerate(mandatory_bitsets):
        covering_gfs = set(ind for ind in comb if (1 << ind) & bitset)
        if len(covering_gfs) >= 1:
            for gf in covering_gfs:
                gf2fam[gf] = family_index
            if family_index not in covered_families:
                covered_mandatory += 1
                covered_families.add(family_index)

    covered_accessory = 0
    for family_index, bitset in enumerate(accessory_bitsets, start=len(mandatory_bitsets)):
        covering_gfs = set(ind for ind in comb if (1 << ind) & bitset)
        if len(covering_gfs) >= 1:
            for gf in covering_gfs:
                gf2fam[gf] = family_index
            if family_index not in covered_families:
                covered_accessory += 1
                covered_families.add(family_index)

    return gf2fam, covered_mandatory, covered_accessory


def find_combinations(matrix: pd.DataFrame, func_unit: FuncUnit) -> List[Tuple[Set[str], Dict[str, str]]]:
    """
    Search working combination of gene families that respect families presence absence model rules

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        List: List of working combination
            Set[str]: Set of selected gene families name that correspond to a working combination
            Dict[str, str]: Association between gene families and families for the combination.
    """

    mandatory = {fam.name for fam in func_unit.mandatory}
    mandatory_indices = [i for i, fam in enumerate(matrix.index.values) if fam in mandatory]

    if len(mandatory_indices) < func_unit.min_mandatory or matrix.shape[0] < func_unit.min_total:
        return []

    bitsets = [bitset_from_row(matrix.iloc[i, :]) for i in range(matrix.shape[0])]
    mandatory_bitsets = [bitsets[i] for i in mandatory_indices]
    accessory_bitsets = [bitsets[i] for i in range(matrix.shape[0]) if i not in mandatory_indices]

    solutions = []
    for size in sorted(range(func_unit.min_total, matrix.shape[1] + 1), reverse=True):
        for comb in combinations(range(matrix.shape[1]), size):
            gf2fam, covered_mandatory, covered_accessory = search_comb(comb, mandatory_bitsets, accessory_bitsets)
            if (covered_mandatory >= func_unit.min_mandatory and
                    covered_mandatory + covered_accessory >= func_unit.min_total):
                solutions.append(({matrix.columns[i] for i in comb},
                                  {matrix.columns[j]: matrix.index[i] for j, i in gf2fam.items()}))

    return solutions


def is_subset(comb1: Tuple[int, ...], comb2: Tuple[int, ...]) -> bool:
    """
    Check if a combination of gene families is subset of another one

    Args:
        comb1 (Tuple[int, ...]): Combination to check if it's a subset
        comb2 (Tuple[int, ...]): Combination supposed bigger

    Returns:
        Boolean: True if first combination is a subset of the second one, False otherwise
    """
    return set(comb1).issubset(comb2)


def check_needed_families(matrix: pd.DataFrame, func_unit: FuncUnit) -> bool:
    """
    Search if it exists a combination of gene families that respect families presence absence unit rules

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        Boolean: True if it exists, False otherwise
    """

    mandatory = {fam.name for fam in func_unit.mandatory}
    mandatory_indices = [i for i, fam in enumerate(matrix.index.values) if fam in mandatory]

    if len(mandatory_indices) < func_unit.min_mandatory or matrix.shape[0] < func_unit.min_total:
        return False

    bitsets = [bitset_from_row(matrix.iloc[i, :]) for i in range(matrix.shape[0])]
    mandatory_bitsets = [bitsets[i] for i in mandatory_indices]
    accessory_bitsets = [bitsets[i] for i in range(matrix.shape[0]) if i not in mandatory_indices]

    not_valid = []
    for size in sorted(range(func_unit.min_total, matrix.shape[1] + 1), reverse=True):
        for comb in combinations(range(matrix.shape[1]), size):
            # If a larger comb can't respect presence absence rules, so a tinier can not too.
            if not any(is_subset(comb, larger_comb) for larger_comb in not_valid):
                _, covered_mandatory, covered_accessory = search_comb(comb, mandatory_bitsets,
                                                                      accessory_bitsets)
                if (covered_mandatory >= func_unit.min_mandatory and
                        covered_mandatory + covered_accessory >= func_unit.min_total):
                    return True
                else:
                    not_valid.append(comb)
    return False


def get_metadata_to_families(pangenome: Pangenome, sources: Iterable[str]) -> Dict[str, Dict[str, Set[GeneFamily]]]:
    """
    Retrieves a mapping of metadata to sets of gene families for each metadata source.

    Args:
        pangenome (Pangenome): Pangenome object containing gene families.
        sources (iterable of str): List of metadata source names.

    Returns:
        dict: A dictionary where each metadata source maps to another dictionary of metadata to sets of gene families.
    """
    meta2gfs = {source: defaultdict(set) for source in sources}
    for source in sources:
        for gf in pangenome.gene_families:
            metadata = gf.get_metadata_by_source(source)
            if metadata is not None:
                for meta in metadata.values():
                    meta2gfs[source][meta.protein_name].add(gf)
                    if "secondary_name" in meta.fields and meta.secondary_name != "":
                        for secondary_name in meta.secondary_name.split(','):
                            meta2gfs[source][secondary_name].add(gf)
    return meta2gfs


def dict_families_context(model: Model, families2gfs: Dict[str, Dict[str, Set[GeneFamily]]]
                          ) -> Tuple[Set[GeneFamily], Dict[str, Set[str]], Dict[str, str]]:
    """
    Retrieves all gene families associated with the families in the model.

    Args:
        model (Model): Model containing the families.
        families2gfs (dict): Dictionary of annotated families.

    Returns:
        tuple: A tuple containing:
            - Set[GeneFamily]: Gene families of interest in the functional unit.
            - dict: Dictionary linking gene families to their families.
            - dict: Dictionary linking families to their sources.
    """
    gene_families = set()
    gf2fam = defaultdict(set)
    fam2source = {}

    for family in model.families:
        for source, fam2gfs in families2gfs.items():
            if family.name in fam2gfs:
                for gf in fam2gfs[family.name]:
                    gene_families.add(gf)
                    gf2fam[gf.name].add(family.name)
                    if family.name in fam2source and fam2source[family.name] != source:
                        logging.getLogger("PANORAMA").warning("Two families have the same protein name for different "
                                                              "sources. First source encountered will be used.")
                    else:
                        fam2source[family.name] = source

        for exchangeable in family.exchangeable:
            for source, fam2gfs in families2gfs.items():
                if exchangeable in fam2gfs:
                    for gf in fam2gfs[exchangeable]:
                        gene_families.add(gf)
                        gf2fam[gf.name].add(family.name)
                        if family.name in fam2source and fam2source[family.name] != source:
                            logging.getLogger("PANORAMA").warning(
                                "Two families have the same protein name for different "
                                "sources. First source encountered will be used.")
                        else:
                            fam2source[family.name] = source
    return gene_families, gf2fam, fam2source