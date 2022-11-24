#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations

# installed libraries
import networkx as nx
import re
from ppanggolin.context.searchGeneContext import compute_gene_context_graph

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.annotation.systems import Systems, System, FuncUnit


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
            for pan_fam in annot2fam[fam_sys.name]:
                families[pan_fam.name] = pan_fam
                fam2annot[pan_fam.name] = fam_sys
    return families, fam2annot


def bool_condition(system: System, func_unit: FuncUnit, list_mandatory: list, group: set, pred_dict: dict,
                   count_forbidden: int):
    bool_list = [False, False, False, False]
    if len(list_mandatory) <= len(func_unit.families.keys()) - min_mandatory_func_unit(system, func_unit):
        bool_list[0] = True
    if group not in pred_dict.values():
        bool_list[1] = True
    if count_forbidden <= max_forbidden_func_unit(system, func_unit):
        bool_list[2] = True
    if len(group) - count_forbidden >= min_mandatory_func_unit(system, func_unit):
        bool_list[3] = True
    return bool_list


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


    def check_cc(cc: set, fam2annot: dict, system: System, func_unit: FuncUnit):
        count_forbidden, count_mandatory, count_accesory = (0, 0, 0)
        forbidden_list, mandatory_list, accessory_list = (func_unit.forbidden_name(), func_unit.mandatory_name(),
                                                          func_unit.accessory_name())
        for node in cc:
            annot = fam2annot.get(node.name)
            if annot is not None:
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
            if check_cc(cc, fam2annot, system, func_unit):
                pred_dict[nb_pred] = cc
                nb_pred += 1
            else:
                remove_node |= cc
    # if len(pred_dict) > 0:
    #     nx.draw(g)
    #     plt.show()
    g.remove_nodes_from(remove_node)
    return pred_dict


def launch_system_search(system: System, annot2fam: dict):
    for func_unit in system.func_units:
        pred_dict = {}
        nb_pred = 0
        if func_unit.parameters['min_total'] == 1:
            nb_pred = search_fu_with_one_fam(func_unit, annot2fam, pred_dict, nb_pred)
        families, fam2annot = dict_families_context(func_unit, annot2fam)
        g = compute_gene_context_graph(families, func_unit.parameters["max_separation"] + 1, disable_bar=True)
        pred_dict = verify_param(g, fam2annot, system, func_unit, pred_dict, nb_pred)
        if len(pred_dict) > 0:
            return pred_dict


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("annotation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_annot(parser)
    return parser


def parser_annot(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--tsv', type=Path, nargs='?',
                                help='Gene families annotation in TSV file. See our github for more detail about format')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?',
                                help="File with all HMM or a directory with one HMM by file")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--prediction_size", required=False, type=int, default=1,
                          help="Number of prediction associate with gene families")
    optional.add_argument('-s', '--systems', required=False, type=Path, default=None,
                          help="Path to systems directory")
    optional.add_argument("--meta", required=False, type=Path, default=None,
                          help="Metadata link to HMM with protein name, description and cutoff")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_annot(main_parser)
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
