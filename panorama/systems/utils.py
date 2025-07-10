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


def filter_global_context(graph: nx.Graph, jaccard_threshold: float = 0.8) -> nx.Graph[GeneFamily]:

    new_graph = nx.Graph()
    # Copy all nodes from the original graph to the new graph
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        # Calculate Jaccard gene proportions for both families
        f1_proportion = len(data["genomes"]) / len(set(f1.organisms))
        f2_proportion = len(data["genomes"]) / len(set(f2.organisms))

        # Update the local copy of the edge data
        data.update({
            'f1': f1.name,
            'f2': f2.name,
            'f1_jaccard_gene': f1_proportion,
            'f2_jaccard_gene': f2_proportion
        })

        # Add the edge to the new graph if it meets the Jaccard threshold
        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            # Add the edge with the updated edge data
            new_graph.add_edge(f1, f2, **data)

    return new_graph


def filter_local_context(graph: nx.Graph, organisms: Set[Organism],
                         jaccard_threshold: float = 0.8) -> nx.Graph[GeneFamily]:
    """
    Filters a graph based on a local Jaccard index.

    Args:
        graph (nx.Graph): A sub-pangenome graph.
        organisms (Set[Organism]): Organisms where edges between families of interest exist. Default is None
        jaccard_threshold (float, optional): Minimum Jaccard similarity used to filter edges between gene families. Default is 0.8.
    """
    def get_gene_proportion(gene_family: GeneFamily) -> float:
        """Returns the Jaccard gene proportion for a given gene family."""
        # Compute the proportion if not cached
        gf_orgs = set(gene_family.organisms).intersection(organisms)
        if len(gf_orgs) > 0:
            return len(data["genomes"].intersection(organisms)) / len(gf_orgs)
        else:
            return 0

    new_graph = nx.Graph()
    # Copy all nodes from the original graph to the new graph
    new_graph.add_nodes_from(graph.nodes(data=True))

    for f1, f2, data in graph.edges(data=True):
        # Calculate Jaccard gene proportions for both families
        f1_proportion = get_gene_proportion(f1)
        f2_proportion = get_gene_proportion(f2)

        # Update the local copy of the edge data
        data.update({
            'f1': f1.name,
            'f2': f2.name,
            'f1_jaccard_gene': f1_proportion,
            'f2_jaccard_gene': f2_proportion
        })

        # Add the edge to the new graph if it meets the Jaccard threshold
        if f1_proportion >= jaccard_threshold and f2_proportion >= jaccard_threshold:
            # Add the edge with the updated edge data
            new_graph.add_edge(f1, f2, **data)

    return new_graph


def check_for_families(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]],
                       mod_fam2meta_source: Dict[str, str], func_unit: FuncUnit
                       ) -> Tuple[bool, Dict[GeneFamily, Tuple[str, int]]]:
    """
    Checks if there are forbidden conditions in the families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[str, Set[Family]]): Dictionary linking gene families to model families.
        func_unit (FuncUnit): Functional unit to check against.

    Returns:
        bool: True if forbidden conditions are encountered, False otherwise.
        dict: Dictionary mapping gene families to metadata information.
    """
    def fam_sort_key(item):
        family, meta_info = item
        score = meta_info[2]
        presence_priority = {"mandatory": 0, "accessory": 1, "forbidden": 2}
        presence_rank = presence_priority.get(family.presence, 3)
        return (-score, presence_rank, family.name) # negative score for descending order

    gf2meta_info = {}
    mandatory_seen = accessory_seen = set()
    
    for gf in sorted(gene_families, key=lambda n: len(gene_fam2mod_fam[n])):
        fam2meta_info = {}
        for family in gene_fam2mod_fam[gf]:
            avail_name = {family.name}.union(family.exchangeable)
            # if family.presence in ("mandatory", "accessory", "forbidden"):
            for meta_id, metadata in gf.get_metadata_by_source(mod_fam2meta_source[family.name]).items():
                if (metadata.protein_name in avail_name or ("secondary_name" in metadata.fields and
                    any(name in avail_name for name in metadata.secondary_name.split(",")))):
                    fam2meta_info[family] = (mod_fam2meta_source[family.name], meta_id, metadata.score)

        sorted_fam2meta_info = sorted(fam2meta_info.items(), key=fam_sort_key)
        family, meta_info = sorted_fam2meta_info[0]
        
        if family.presence == "forbidden":
            return False, {}
        
        if family.presence == "mandatory" and family.name not in mandatory_seen:
            mandatory_seen.add(family.name)
        elif family.presence == "accessory" and family.name not in accessory_seen:
            accessory_seen.add(family.name)
        gf2meta_info[gf] = meta_info[:-1]
        
    if (len(mandatory_seen) >= func_unit.min_mandatory and 
        len(mandatory_seen | accessory_seen) >= func_unit.min_total):
        return True, gf2meta_info
    return False, {}


def get_gfs_matrix_combination(gene_families: Set[GeneFamily], gene_fam2mod_fam: Dict[GeneFamily, Set[Family]]) -> pd.DataFrame:
    """
    Build a matrix of association between gene families and families.

    Args:
        gene_families (Set[GeneFamily]): Set of gene families.
        gene_fam2mod_fam (Dict[GeneFamily, Set[Family]): Dictionary linking gene families to model families.

    Returns:
        pd.DataFrame: Matrix of association between gene families and families.
    """

    gfs = set()
    model_fams = set()
    for gf in gene_families:
        for fam in gene_fam2mod_fam[gf]:
            if fam.presence in ["mandatory", "accessory"]:
                gfs.add(gf)
                model_fams.add(fam.name)
                
    gfs = list(gfs)
    model_fams = list(model_fams)
    score_matrix = np.zeros((len(model_fams), len(gfs)), dtype=int)

    for j, gf in enumerate(gfs):
        annotations = [f.name for f in gene_fam2mod_fam[gf]]
        for i, fam in enumerate(model_fams):
            if fam in annotations:
                score_matrix[i, j] = 1 
                
    return pd.DataFrame(score_matrix, index=model_fams, columns=[gf.name for gf in gfs])


def check_needed_families(matrix: pd.DataFrame, func_unit: FuncUnit) -> bool:
    """
    Check if there are enough mandatory and total families to satisfy the functional unit rules.

    Args:
        matrix: The association matrix between gene families and families
        func_unit: The functional unit to search for.

    Returns:
        Boolean: True if satisfied, False otherwise
    """

    matrix = matrix.loc[matrix.sum(axis=1) > 0] # Remove all-zero rows
    
    mandatory_fams = {fam.name for fam in func_unit.mandatory}
    mandatory_count = sum(fam in mandatory_fams for fam in matrix.index)

    return (mandatory_count >= func_unit.min_mandatory and
            len(matrix) >= func_unit.min_total)


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
                    if "secondary_name" in meta.fields and meta.secondary_name != "":
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
        exchangeables = fam_model.exchangeable | {fam_model.name} 
        for exchangeable in exchangeables: 
            for source, annotation2families in annot2fam.items():
                if exchangeable in annotation2families:
                    
                    if fam_model.name in fam2source and fam2source[fam_model.name] != source:
                        logging.getLogger("PANORAMA").warning(f"Protein annotation {fam_model.name} is encountered in multiple sources." 
                                                               "All sources will be used, but only first one will be associated with " 
                                                               "the model family.")
                    else:
                        fam2source[fam_model.name] = source
                    
                    for gf in annotation2families[exchangeable]:
                        gene_families.add(gf)
                        gf2fam[gf].add(fam_model)     

    return gene_families, gf2fam, fam2source




# The following functions are not used anywhere anymore

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