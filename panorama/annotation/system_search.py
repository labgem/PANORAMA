#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations

# installed libraries
import networkx as nx
import matplotlib.pyplot as plt
import re
from tqdm import tqdm
from ppanggolin.context.searchGeneContext import compute_gene_context_graph
from ppanggolin.genome import Gene

# local libraries
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.annotation.rules import Systems, System, FuncUnit


def min_mandatory_func_unit(system: System, func_unit: FuncUnit):
    if func_unit.parameters["min_mandatory"] is not None:
        return func_unit.parameters["min_mandatory"]
    else:
        return system.parameters["min_mandatory"]


def max_forbidden_func_unit(system: System, func_unit: FuncUnit):
    if func_unit.parameters["max_forbidden"] is not None:
        return func_unit.parameters["max_forbidden"]
    else:
        return system.parameters["max_forbidden"]


def min_total_func_unit(system: System, func_unit: FuncUnit):
    if func_unit.parameters["min_total"] is not None:
        return func_unit.parameters["min_total"]
    else:
        return system.parameters["min_total"]


def list_mandatory_family(func_unit: FuncUnit):
    list_mandatory = []
    for fam in func_unit.mandatory.keys():
        list_mandatory.append(fam)
    return list_mandatory


def dict_families_context(func_unit: FuncUnit, annot2fam: dict) -> (dict, dict):
    """
    Recover all families in the function unit
    :param func_unit: function unit object of system
    :param annot2fam: dictionary of annotated families
    :return families: dictionary of families interesting
    """
    families = dict()
    fam2annot = dict()
    for fam_sys in func_unit.families.keys():
        # families.update({annot: pan_fam for annot, pan_fam in annot2fam.items() if re.match(f"^{fam_sys}", annot)})
        for annot, pan_fam in annot2fam.items():
            if re.match(f"^{fam_sys}", annot):
                for family in pan_fam:
                    families[family.name] = family
                    fam2annot[family.name] = fam_sys
    return families, fam2annot


def bool_condition(system: System, func_unit: FuncUnit, list_mandatory: list, group: set, pred_dict: dict,
                   count_forbidden: int):
    bool_list = [False, False, False, False]
    if len(list_mandatory) <= len(func_unit.families.keys())-min_mandatory_func_unit(system, func_unit):
        bool_list[0] = True
    if group not in pred_dict.values():
        bool_list[1] = True
    if count_forbidden <= max_forbidden_func_unit(system, func_unit):
        bool_list[2] = True
    if len(group)-count_forbidden >= min_mandatory_func_unit(system, func_unit):
        bool_list[3] = True
    return bool_list




def verify_param(g: nx.Graph(), fam2annot: dict, system: System, func_unit: FuncUnit):
    """
    Verify parameters

    """
    pred_dict = dict()  # dictionnaire pour les systèmes prédits
    nb_pred = 0  # compteur de prédiction du system

    for node in g.nodes():  # pour chaque noeud du graphe
        group = {node.name}  # groupe d'un ensemble de noeud du graphe
        count_forbidden = 0
        list_mandatory = list_mandatory_family(func_unit)
        if fam2annot.get(node.name) in func_unit.forbidden.keys():  # if node is forbidden
            count_forbidden += 1
        elif fam2annot.get(node.name) in list_mandatory:
            list_mandatory.remove(fam2annot.get(node.name))
        if len(node.edges) >= min_mandatory_func_unit(system, func_unit) - 1 and \
                count_forbidden <= max_forbidden_func_unit(system, func_unit):  # si nb arête plus petit que min_m et compteur f plus petit que max_f
            for edge in node.edges:     # pr chaque arête
                neighbor = edge.source if edge.source != node else edge.target
                if fam2annot.get(neighbor.name) in func_unit.forbidden.keys(): 
                    count_forbidden += 1
                elif fam2annot.get(neighbor.name) in list_mandatory:
                    list_mandatory.remove(fam2annot.get(neighbor.name))
                if fam2annot.get(neighbor.name) in func_unit.families.keys():
                    group.add(neighbor.name)
        if all(bool_condition(system, func_unit, list_mandatory, group, pred_dict, count_forbidden)):
            nb_pred += 1
            pred_dict[nb_pred] = group



def launch_system_search(system: System, annot2fam: dict):
    for func_unit in system.func_units.values():
        families, fam2annot = dict_families_context(func_unit, annot2fam)
        g = compute_gene_context_graph(families, func_unit.parameters["max_separation"])
        verify_param(g, fam2annot, system, func_unit)
        nx.draw(g)
        plt.show()
