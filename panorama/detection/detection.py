#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import argparse
from pathlib import Path
import logging
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from itertools import combinations

# installed libraries
import networkx as nx
import pandas as pd
from ppanggolin.context.searchGeneContext import compute_gene_context_graph
import matplotlib.pyplot as plt

# local libraries
from panorama.detection.systems import Systems, System, FuncUnit
from panorama.utils import check_tsv_sanity
from panorama.format.read_binaries import check_pangenome_info
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome


def check_pangenome_detection(pangenome: Pangenome, source: str, force: bool = False, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param disable_bar: Disable bar
    """
    if source in pangenome.status["annotation_source"]:
        check_pangenome_info(pangenome, need_annotations=True, need_families=True,
                             need_annotation_fam=True, need_modules=True, disable_bar=disable_bar)
    else:
        raise Exception("Annotation source not in pangenome")


def read_systems(systems_path: Path, disable_bar: bool = False) -> Systems:
    """ Read all json files systems in the directory

    :param systems_path: path of systems directory
    """
    systems = Systems()
    for file in tqdm(list(systems_path.glob("*.json")), unit='system', desc="Read system", disable=disable_bar):
        with open(file.resolve().as_posix()) as json_file:
            data = json.load(json_file)
            system = System()
            try:
                system.read_system(data)
            except Exception:
                raise Exception(f"Problem to read json {file}")
            else:
                systems.add_sys(system)
    return systems


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
    for fam in func_unit.mandatory:
        list_mandatory.append(fam)
    return list_mandatory


def get_annotation_to_families(pangenome: Pangenome, source: str):
    annot2fam = {}
    for gf in pangenome.gene_families:
        annotations = gf.get_source(source)
        if annotations is not None:
            for annotation in annotations:
                if annotation.name in annot2fam:
                    annot2fam[annotation.name].add(gf)
                else:
                    annot2fam[annotation.name] = {gf}
    return annot2fam


def dict_families_context(func_unit: FuncUnit, annot2fam: dict) -> (dict, dict):
    """
    Recover all families in the function unit
    :param func_unit: function unit object of system
    :param annot2fam: dictionary of annotated families
    :return families: dictionary of families interesting
    """
    families = dict()
    fam2annot = dict()
    for fam_sys in func_unit.families:
        # families.update({annot: pan_fam for annot, pan_fam in annot2fam.items() if re.match(f"^{fam_sys}", annot)})
        if fam_sys.name in annot2fam:
            for gf in annot2fam[fam_sys.name]:
                families[gf.name] = gf
                if gf.name in fam2annot:
                    fam2annot[gf.name].add(fam_sys)
                else:
                    fam2annot[gf.name] = {fam_sys}
    return families, fam2annot


# def bool_condition(system: System, func_unit: FuncUnit, list_mandatory: list, group: set, pred_dict: dict,
#                    count_forbidden: int):
#     bool_list = [False, False, False, False]
#     if len(list_mandatory) <= len(func_unit.families.keys()) - min_mandatory_func_unit(system, func_unit):
#         bool_list[0] = True
#     if group not in pred_dict.values():
#         bool_list[1] = True
#     if count_forbidden <= max_forbidden_func_unit(system, func_unit):
#         bool_list[2] = True
#     if len(group) - count_forbidden >= min_mandatory_func_unit(system, func_unit):
#         bool_list[3] = True
#     return bool_list


def search_fu_with_one_fam(func_unit: FuncUnit, annot2fam: dict, pred_dict: dict, nb_pred: int):
    for mandatory_fam in func_unit.mandatory:
        if mandatory_fam.name in annot2fam:
            for pan_fam in annot2fam[mandatory_fam.name]:
                pred_dict[nb_pred] = {pan_fam}
                nb_pred += 1
    return nb_pred


def verify_param(g: nx.Graph(), fam2annot: dict, system: System, func_unit: FuncUnit, pred_dict: dict, nb_pred: int):
    """
    Verify parameters

    """
    seen = set()

    def extract_cc(node: GeneFamily, graph: nx.Graph, seen: set):
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

    def check_cc(cc: set, fam2annot: dict, func_unit: FuncUnit):
        count_forbidden, count_mandatory, count_accesory = (0, 0, 0)
        forbidden_list, mandatory_list, accessory_list = (func_unit.forbidden_name(), func_unit.mandatory_name(),
                                                          func_unit.accessory_name())
        for node in cc:
            annotations = fam2annot.get(node.name)
            if annotations is not None:
                for annot in annotations:
                    if annot.type == 'forbidden':  # if node is forbidden
                        if annot.name in forbidden_list:
                            count_forbidden += 1
                            forbidden_list.remove(annot.name)
                            if count_forbidden > func_unit.parameters["max_forbidden"]:
                                return False
                    elif annot.type == 'mandatory':  # if node is mandatory
                        if annot.name in mandatory_list:
                            count_mandatory += 1
                            mandatory_list.remove(annot.name)
                    if annot.type == 'accessory':  # if node is accessory
                        if annot.name in accessory_list:
                            count_accesory += 1
                            accessory_list.remove(annot.name)
        if count_mandatory >= func_unit.parameters['min_mandatory'] and \
                count_accesory + count_mandatory >= func_unit.parameters['min_total']:
            return True
        else:
            return False

    remove_node = set()
    for node in g.nodes():  # pour chaque noeud du graphe
        if node not in seen:
            cc = extract_cc(node, g, seen)  # extract coonnect component
            if check_cc(cc, fam2annot, func_unit):
                pred_dict[nb_pred] = cc
                nb_pred += 1
            else:
                remove_node |= cc
    g.remove_nodes_from(remove_node)
    return pred_dict


def search_system(system: System, annot2fam: dict):
    if system.name in ['disarm_type_I']:
        print("pika")
    for func_unit in system.func_units:
        pred_dict = {}
        nb_pred = 0
        if func_unit.parameters['min_total'] == 1:
            nb_pred = search_fu_with_one_fam(func_unit, annot2fam, pred_dict, nb_pred)
        families, fam2annot = dict_families_context(func_unit, annot2fam)
        g = compute_gene_context_graph(families, func_unit.parameters["max_separation"] + 1, disable_bar=True)
        # nx.draw(g)
        # plt.show()
        pred_dict = verify_param(g, fam2annot, system, func_unit, pred_dict, nb_pred)
        if len(pred_dict) > 0:
            return pred_dict, system.name


def system_to_module(pangenome: Pangenome, system: System, predictions: set):
    """Associate a system to modules"""
    for module in pangenome.modules:
        for prediction in predictions:
            if prediction.issubset(module.families):
                print("pika")
    pass


def project_system(pangenome: Pangenome, proj_dict: dict, pred_res: dict, system: System, source: str):
    mandatory_family = [fam.name for fam in system.families if fam.type == 'mandatory']
    accessory_family = [fam.name for fam in system.families if fam.type == 'accessory']
    for fu in system.func_units:
        for organism in pangenome.organisms:
            # if organism.name in ['GCF_000472745.1_ASM47274v1_genomic'] and system.name == 'DRT_Other':
            #     print("pika")
            for prediction in pred_res.values():
                for combi_pred in combinations(prediction, fu.parameters["min_total"]):
                    if set(combi_pred).issubset(organism.families):
                        mandatory_dict = {mandatory: False for mandatory in mandatory_family}
                        accessory_dict = {accessory: False for accessory in accessory_family}
                        for gf in prediction:
                            if gf.get_source(source) is not None:
                                for annotation in [annot.name for annot in gf.get_source(source)]:
                                    if annotation in mandatory_family:
                                        mandatory_dict[annotation] = True
                                    elif annotation in accessory_family:
                                        accessory_dict[annotation] = True
                        if sum(mandatory_dict.values()) >= fu.parameters['min_mandatory']:
                            if sum(mandatory_dict.values()) + sum(accessory_dict.values()) >= fu.parameters['min_total']:
                                if organism.name in proj_dict:
                                    proj_dict[organism.name].add(system.name)
                                else:
                                    proj_dict[organism.name] = {system.name}
            # print("pikapi")
    # for keys, value in pred_res.items():
    #     for gf_combi in combinations(value, 2):
    #         org_inter = set()
    #         for fam in gf_combi:
    #             if len(org_inter) == 0:
    #                 org_inter = fam.organisms
    #             else:
    #                 org_inter = org_inter.intersection(fam.organisms)
    #         for org in org_inter:
    #             if org.name in ['GCF_006539645.1_ASM653964v1_genomic', 'GCF_001038185.1_ASM103818v1_genomic',
    #                             'GCF_000807315.1_ASM80731v1_genomic', 'GCF_014194605.1_ASM1419460v1_genomic',
    #                             'GCF_009864815.1_ASM986481v1_genomic', 'GCF_013752735.1_ASM1375273v1_genomic',
    #                             'GCF_000472985.1_ASM47298v1_genomic', 'GCF_000773865.1_ASM77386v1_genomic',
    #                             'GCF_000284375.1_ASM28437v1_genomic']:
    #                 print("pika")
    #             gf_present = set()
    #             for contig in org.contigs:
    #                 for gene in contig.genes:
    #                     if gene.family in gf_combi:
    #                         gf_present.add(gene.family)
    #             mandatory_dict = {mandatory: False for mandatory in mandatory_family}
    #             for gf in gf_present:
    #                 for annotation in [annot.name for annot in gf.get_source(source)]:
    #                     if annotation in mandatory_family:
    #                         mandatory_dict[annotation] = True
    #             if all(x for x in mandatory_dict.values()):
    #                 if org.name in proj_dict:
    #                     proj_dict[org.name].add(system.name)
    #                 else:
    #                     proj_dict[org.name] = {system.name}


def search_systems(systems: Systems, pangenome: Pangenome, source: str, threads: int = 1, disable_bar: bool = False):
    """
    Search present system in the pangenome
    :param systems:
    :param pangenome:
    :param source:
    :param threads:
    :param disable_bar:

    :return:
    """
    proj_dict = {}
    annot2fam = get_annotation_to_families(pangenome=pangenome, source=source)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=systems.size, unit='system', disable=disable_bar) as progress:
            futures = []
            for system in systems:
                future = executor.submit(search_system, system, annot2fam)
                futures.append(future)

            prediction_results = {}
            for future in futures:
                result = future.result()
                future.add_done_callback(lambda p: progress.update())
                if result is not None:
                    prediction_results[result[1]] = result[0]
                    # system_to_module(pangenome, systems.get_sys(result[1]), result[0].values())
                    project_system(pangenome, proj_dict, result[0], systems.get_sys(result[1]), source)
    proj = pd.DataFrame.from_dict(proj_dict, orient='index')
    sys_df = pd.DataFrame(prediction_results, columns=['System',
                                                       'Nb Detection']).sort_values('System').reset_index(drop=True)
    proj.to_csv("projection5.tsv", sep="\t", index=['Organisms'], index_label='Organisms',
                header=[f"System {i}" for i in range(1, proj.shape[1] + 1)])
    sys_df.to_csv("system5.tsv", sep="\t", header=['System', "Nb_detected"])


def launch(args):
    """
    Launch functions to detect systems in pangenomes

    :param args: Argument given
    """
    pan_to_path = check_tsv_sanity(args.pangenomes)
    systems = read_systems(args.systems)
    for pangenome_name, pangenome_info in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
        pangenome.add_file(pangenome_info["path"])
        check_pangenome_detection(pangenome, source=args.source, force=args.force, disable_bar=args.disable_prog_bar)
        res = search_systems(systems, pangenome, args.source, args.threads, args.disable_prog_bar)
        logging.getLogger().info(f"Write Annotation in pangenome {pangenome_name}")
        # write_pangenome(pangenome, pangenome_info["path"], source=args.source, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("annotation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_detection(parser)
    return parser


def parser_detection(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument('-s', '--systems', required=False, type=Path, default=None,
                          help="Path to systems directory")
    required.add_argument("-S", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source where panorama as to select in pangenomes')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_detection(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
