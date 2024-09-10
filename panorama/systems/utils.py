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
# from pulp import LpMaximize, LpProblem, LpVariable, lpSum, LpBinary, PULP_CBC_CMD
from ppanggolin.genome import Organism
from ppanggolin.metadata import Metadata

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family
from panorama.pangenomes import Pangenome

pd.options.mode.copy_on_write = True  # Remove when pandas3.0 available. See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy


def filter_local_context(graph: nx.Graph, organisms: Set[Organism],
                         jaccard_threshold: float = 0.8) -> None:
    """
    Filters a graph based on a local Jaccard index.

    Args:
        graph (nx.Graph): A sub-pangenome graph.
        organisms (Set[Organism]): Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """

    # Compute local Jaccard
    edges2remove = set()
    for f1, f2, data in graph.edges(data=True):
        f1_orgs = set(f1.organisms).intersection(organisms)
        f2_orgs = set(f2.organisms).intersection(organisms)
        f1_gene_proportion = len(data["genomes"]) / len(f1_orgs) if len(f1_orgs) > 0 else 0
        f2_gene_proportion = len(data["genomes"]) / len(f2_orgs) if len(f2_orgs) > 0 else 0

        data['f1'] = f1.name
        data['f2'] = f2.name
        data['f1_jaccard_gene'] = f1_gene_proportion
        data['f2_jaccard_gene'] = f2_gene_proportion

        if not ((f1_gene_proportion >= jaccard_threshold) and (f2_gene_proportion >= jaccard_threshold)):
            edges2remove.add((f1, f2))

    graph.remove_edges_from(edges2remove)


def check_for_forbidden_families(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                                 func_unit: FuncUnit) -> bool:
    """
    Checks if there are forbidden conditions in the families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[str, Set[Family]]): Dictionary linking gene families to model families.
        func_unit (FuncUnit): Functional unit to check against.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
    """
    forbidden_list = list(map(lambda x: x.name, func_unit.forbidden))

    def get_number_of_mod_fam(gene_family: GeneFamily) -> int:
        """Gets the number of model families associated with the gene family.

        Args:
            gene_family (GeneFamily): The gene family of interest.

        Returns:
            int: The number of model families associated with the gene family.
        """
        model_families = gene_fam2mod_fam.get(gene_family.name)
        if model_families is not None:
            return len(model_families)
        else:
            return 0

    for node in sorted(gene_families, key=lambda n: get_number_of_mod_fam(n)):
        if node in gene_fam2mod_fam:
            for family in gene_fam2mod_fam[node]:
                if family.presence == 'forbidden' and family.name in forbidden_list:  # if node is forbidden
                    return True
    return False


def get_gfs_matrix_combination(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                               mod_fam2meta_source: Dict[str, str]
                               ) -> Tuple[
    pd.DataFrame, Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]], Dict[str, Dict[GeneFamily, Tuple[int, Metadata]]]]:
    """
    Build a matrix of association between gene families and families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]): Dictionary linking gene families to model families.
        mod_fam2meta_source (Dict[str, str]): Dictionary linking families to metadata sources.

    Returns:
        pd.DataFrame: Matrix of association between gene families and families.
        Dict[GeneFamily, Tuple[int, Metadata]]: Dictionary linking gene families metadata to mandatory families
        Dict[GeneFamily, Tuple[int, Metadata]]: Dictionary linking gene families metadata to accessory families
    """

    def add_metadata_to_dict(presence_gfs2metadata):
        """
        Associate gene families to metadata in a dictionary

        Args:
            presence_gfs2metadata: Dictionary to save gene families and their metadata
        """
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

    mandatory_gfs2metadata = defaultdict(dict)
    accessory_gfs2metadata = defaultdict(dict)
    gfs = set()
    for gene_family in gene_families:
        for family in gene_fam2mod_fam[gene_family]:
            avail_name = {family.name}.union(family.exchangeable)
            if family.presence == "mandatory":
                add_metadata_to_dict(mandatory_gfs2metadata)
                gfs.add(gene_family)
            elif family.presence == "accessory":
                add_metadata_to_dict(accessory_gfs2metadata)
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


# def greedy_algorithm(matrix: pd.DataFrame, func_unit: FuncUnit) -> Tuple[Set[str], Set[str]]:
#     """
#     Find a working combination of gene families covering the needed mandatory and accessory families
#
#     Args:
#         matrix: The association matrix between gene families and families
#         func_unit: The functional unit to search for.
#
#     Returns:
#         Set[str]: Set of selected gene families name that correspond to a working combination
#         Set[str]: Set of covered families names that correspond to a working combination
#     """
#     covered_families = set()
#     selected_gfs = set()
#
#     mandatory = {fam.name for fam in func_unit.mandatory}
#     min_needed = False
#     while not min_needed:
#         best_gf = None
#         best_coverage = 0
#
#         for gf in matrix.columns:
#             if gf in selected_gfs:
#                 continue
#             coverage = len([family for family in matrix.index
#                             if matrix.loc[family, gf] == 1 and family not in covered_families])
#
#             if coverage > best_coverage:
#                 best_coverage = coverage
#                 best_gf = gf
#
#         if best_gf is None:
#             break
#
#         selected_gfs.add(best_gf)
#         for family in matrix.index:
#             if matrix.loc[family, best_gf] == 1:
#                 covered_families.add(family)
#         if (len(covered_families) >= func_unit.min_total and
#                 len(covered_families.intersection(mandatory)) >= func_unit.min_mandatory):
#             min_needed = True
#     if min_needed:
#         return selected_gfs, covered_families
#     else:
#         return set(), set()


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

    def bitset_from_row(row) -> int:
        """
        Converts a row of the binary matrix into an integer bitset.
        Each bit represents the presence of an individual for that characteristic.
        Args:
            row: A row of the binary matrix

        Returns:
            int: integer bitset representing the presence of an individual for that characteristic
        """
        bset = 0
        for i, val in enumerate(row):
            if val == 1:
                bset |= (1 << i)  # Met le i-ème bit à 1
        return bset

    num_individuals = matrix.shape[1]
    num_features = matrix.shape[0]

    mandatory = {fam.name for fam in func_unit.mandatory}
    mandatory_indices = [i for i, fam in enumerate(matrix.index.values) if fam in mandatory]

    if len(mandatory_indices) < func_unit.min_mandatory or matrix.shape[0] < func_unit.min_total:
        return []

    bitsets = [bitset_from_row(matrix.iloc[i, :]) for i in range(matrix.shape[0])]
    mandatory_bitsets = [bitsets[i] for i in mandatory_indices]
    accessory_bitsets = [bitsets[i] for i in range(num_features) if i not in mandatory_indices]

    solutions = []
    for size in sorted(range(func_unit.min_total, num_individuals + 1), reverse=True):
        for comb in combinations(range(num_individuals), size):
            covered_features = set()
            association = {}

            covered_mandatory = 0
            for feature_index, bitset in enumerate(mandatory_bitsets):
                covering_individuals = set(ind for ind in comb if (1 << ind) & bitset)
                if len(covering_individuals) >= 1:
                    for ind in covering_individuals:
                        association[ind] = feature_index
                    if feature_index not in covered_features:
                        covered_mandatory += 1
                        covered_features.add(feature_index)

            covered_accessory = 0
            for feature_index, bitset in enumerate(accessory_bitsets, start=len(mandatory_bitsets)):
                covering_individuals = set(ind for ind in comb if (1 << ind) & bitset)
                if len(covering_individuals) >= 1:
                    for ind in covering_individuals:
                        association[ind] = feature_index
                    if feature_index not in covered_features:
                        covered_accessory += 1
                        covered_features.add(feature_index)

            if (covered_mandatory >= func_unit.min_mandatory and
                    covered_mandatory + covered_accessory >= func_unit.min_total):
                solutions.append(({matrix.columns[i] for i in comb},
                                  {matrix.columns[j]: matrix.index[i] for j, i in association.items()}))

    return solutions


def check_needed_families(matrix: pd.DataFrame, func_unit: FuncUnit) -> bool:
    """
    Search if it exists a combination of gene families that respect families presence absence unit rules

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        Boolean: True if it exists, False otherwise
    """

    def is_subset(comb1, comb2):
        """
        Vérifie si comb1 est un sous-ensemble de comb2.
        """
        return set(comb1).issubset(comb2)

    def bitset_from_row(row) -> int:
        """
        Converts a row of the binary matrix into an integer bitset.
        Each bit represents the presence of an individual for that characteristic.
        Args:
            row: A row of the binary matrix

        Returns:
            int: integer bitset representing the presence of an individual for that characteristic
        """
        bset = 0
        for i, val in enumerate(row):
            if val == 1:
                bset |= (1 << i)  # Met le i-ème bit à 1
        return bset

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
                covered_features = set()
                association = {}

                covered_mandatory = 0
                for feature_index, bitset in enumerate(mandatory_bitsets):
                    covering_individuals = set(ind for ind in comb if (1 << ind) & bitset)
                    if len(covering_individuals) >= 1:
                        for ind in covering_individuals:
                            association[ind] = feature_index
                        if feature_index not in covered_features:
                            covered_mandatory += 1
                            covered_features.add(feature_index)

                covered_accessory = 0
                for feature_index, bitset in enumerate(accessory_bitsets, start=len(mandatory_bitsets)):
                    covering_individuals = set(ind for ind in comb if (1 << ind) & bitset)
                    if len(covering_individuals) >= 1:
                        for ind in covering_individuals:
                            association[ind] = feature_index
                        if feature_index not in covered_features:
                            covered_accessory += 1
                            covered_features.add(feature_index)

                if (covered_mandatory >= func_unit.min_mandatory and
                        covered_mandatory + covered_accessory >= func_unit.min_total):
                    return True
                else:
                    not_valid.append(comb)
    return False


# def ilp_optimization(matrix: pd.DataFrame, func_unit: FuncUnit):
#     """
#     Integer Linear Programming (ILP) to find an optimal combination of families that maximizes fam coverage
#     while satisfying constraints on mandatory and total annotations.
#
#     Args:
#         matrix (pd.DataFrame): A binary matrix where rows represent annotations and columns represent families.
#         func_unit (int): Minimum total number of annotations (mandatory + accessory) that must be covered.
#
#     Returns:
#         tuple: A tuple containing:
#             - list: Selected families that maximize coverage.
#             - set: Covered annotations based on the selected families.
#             - float: Execution time of the ILP optimization.
#     """
#     # Create a linear programming problem with the objective to maximize
#     model = LpProblem(name="gf-selection", sense=LpMaximize)
#
#     # Decision variables: 1 if the gf is selected, 0 otherwise
#     x = {gf: LpVariable(f"x_{gf}", cat=LpBinary) for gf in matrix.columns}
#
#     mandatory = matrix.index.values.tolist()
#     # Constraint: Each mandatory fam must be covered by at least one gf
#     for family in mandatory:
#         model += (lpSum(matrix.loc[family, gf] * x[gf] for gf in matrix.columns) >= 1, f"mandatory_{family}")
#
#     # Constraint: The total number of covered mandatory annotations must be at least min_mandatory
#     model += (
#         lpSum(matrix.loc[fam, gf] * x[gf] for fam in mandatory for gf in matrix.columns) >= func_unit.min_mandatory,
#         "min_mandatory")
#
#     # Constraint: The total number of covered annotations (mandatory + accessory) must be at least min_total
#     model += (
#         lpSum(matrix.loc[fam, gf] * x[gf] for fam in matrix.index for gf in matrix.columns) >= func_unit.min_total,
#         "min_total")
#
#     # Constraint: Ensure each annotation is covered by at most one family
#     # for fam in matrix.index:
#     #     model += (lpSum(matrix.loc[fam, gf] * x[gf] for gf in matrix.columns) <= 1, f"unique_cover_{fam}")
#
#     # Objective function: Maximize the total number of covered annotations
#     model += lpSum(matrix.loc[fam, gf] * x[gf] for fam in matrix.index for gf in matrix.columns)
#
#     # Solve the ILP problem
#     model.solve(PULP_CBC_CMD(msg=False))
#
#     # Extract the selected families from the solution
#     selected_gfs = [gf for gf in matrix.columns if x[gf].value() == 1]
#
#     # Determine the annotations covered by the selected families
#     coverage = {fam for fam in matrix.index if any(matrix.loc[fam, gf] == 1 for gf in selected_gfs)}
#
#     # Create a mapping of annotations to families that cover them
#     fam_gf_map = {}
#     for family in coverage:
#         covering_families = [gf for gf in selected_gfs if matrix.loc[family, gf] == 1]
#         fam_gf_map[family] = covering_families
#     return selected_gfs, fam_gf_map

def get_metadata_to_families(pangenome: Pangenome, sources: Iterable[str]) -> Dict[str, Dict[str, Set[GeneFamily]]]:
    """
    Retrieves a mapping of metadata to sets of gene families for each metadata source.

    Args:
        pangenome (Pangenome): Pangenome object containing gene families.
        sources (iterable of str): List of metadata source names.

    Returns:
        dict: A dictionary where each metadata source maps to another dictionary of metadata to sets of gene families.
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
        -> Tuple[Set[GeneFamily], Dict[GeneFamily, Set[Family]], Dict[str, str]]:
    """
    Retrieves all gene families associated with the families in the model.

    Args:
        model (Model): Model containing the families.
        annot2fam (dict): Dictionary of annotated families.

    Returns:
        tuple: A tuple containing:
            - Set[GeneFamily]: Gene families of interest in the functional unit.
            - dict: Dictionary linking gene families to their families.
            - dict: Dictionary linking families to their sources.
    """
    gene_families = set()
    gf2fam = defaultdict(set)
    fam2source = {}
    for fam_model in model.families:
        for source, annotation2families in annot2fam.items():
            if fam_model.name in annotation2families:
                for gf in annotation2families[fam_model.name]:
                    gene_families.add(gf)
                    gf2fam[gf].add(fam_model)
                    if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                        logging.getLogger("PANORAMA").warning(
                            "Two families have the same protein name for different "
                            "sources. First source encountered will be used.")
                    else:
                        fam2source[fam_model.name] = source

        for exchangeable in fam_model.exchangeable:
            for source, annotation2families in annot2fam.items():
                if exchangeable in annotation2families:
                    for gf in annotation2families[exchangeable]:
                        gene_families.add(gf)
                        gf2fam[gf].add(fam_model)
                        if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                            logging.getLogger("PANORAMA").warning(
                                "Two families have the same protein name for different "
                                "sources. First source encountered will be used.")
                        else:
                            fam2source[fam_model.name] = source
    return gene_families, gf2fam, fam2source
