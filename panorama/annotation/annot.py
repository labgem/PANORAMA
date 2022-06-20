#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import csv
import json
import re
from pathlib import Path

# installed libraries
import networkx as nx
import matplotlib.pyplot as plt
from ppanggolin.formats.readBinaries import check_pangenome_info

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.annotation.rules import System, Systems
from panorama.pangenomes import Pangenome


def check_pangenome_annotation(pangenome: Pangenome, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param disable_bar: Disable bar
    """
    check_pangenome_info(pangenome, need_families=True, need_graph=True, need_annotations=True, disable_bar=disable_bar)


def read_systems(file: Path, systems=Systems()):
    """ Read all json files in the directory

    :param systems_path: path of systems directory
    :param systems: class Systems with all systems
    """
    # for file in systems_path.glob("*.json"):
    with open(file.resolve().as_posix()) as json_file:
        data = json.load(json_file)
        system = System()
        system.read_system(data)
        systems.add_sys(system)
    return system
    # systems.print_systems()


def parser_hmm(hmm_path: Path) -> dict:
    """
    Recover information of hmm results in a tsv file

    :param hmm_path: Path to hmm results of tsv file

    :raise: Unexpected error when opening the list of hmm

    :return: hmm_res : dictionary of information hmm results
    """
    hmm_res = {}
    try:
        hmm_file = open(hmm_path.absolute(), 'r')
    except IOError:
        raise IOError
    except Exception:
        raise Exception("Unexpected error when opening the list of hmm")
    else:
        hmm = csv.reader(hmm_file, delimiter="\t")
        for line in hmm:
            if not (re.compile("^#").search(line[0])):
                target_name = line[0]
                if line[1] != "-":
                    accession = line[1]
                else:
                    accession = None
                query_name = line[2]
                evalue = line[4]
                if line[18] != "-":
                    system = line[18]
                else:
                    system = target_name
                hmm_res[query_name] = {"name": target_name,
                                       "accession": accession,
                                       "e-value": evalue,
                                       "description": system
                                       }
        return hmm_res


def annot_pangenome(pangenome: Pangenome, hmm: Path, system: Path, source: str = None, force: bool = False,
                    disable_bar: bool = False):
    """ Main function to add annotation to pangenome

    :param pangenome: Pangenome object to ppanggolin
    :param hmm: Path of hmm
    :param system: Path of system
    :param source: Source of hmm
    :param force: Boolean of force argument
    :param disable_bar: Disable bar
    """
    check_pangenome_annotation(pangenome, disable_bar=disable_bar)
    annot2fam = dict()
    for query_name, annot in parser_hmm(hmm).items():
        gene_fam = pangenome.get_gene_family(name=query_name)
        gene_fam.add_annotation(source=hmm.stem if source is None else source, annotation=annot, force=force)
        if not annot["description"] in annot2fam:
            annot2fam[annot["description"]] = set()
        annot2fam[annot["description"]].add(gene_fam)
    return annot2fam


def get_max_separation_sys(system: System):
    """
    Verify the maximum separation condition for genes family in the pangenome

    :param system: one system object
    """
    max_sep_annot = dict()
    for annot_name, fam_obj in system.families.items():
        if fam_obj.parameters is not None and fam_obj.parameters["max_separation"]:
            max_sep_annot[annot_name] = get_max_separation(fam_obj)
        else:
            for fu_name, fu_obj in system.func_units.items():
                if fam_obj.func_unit == fu_name:
                    if fu_obj.parameters is not None and fu_obj.parameters["max_separation"]:
                        max_sep_annot[annot_name] = get_max_separation(fu_obj)
                    elif system.parameters["max_separation"]:
                        max_sep_annot[annot_name] = get_max_separation(system)
    return max_sep_annot


def check_presence_family(system: System, annot2fam: dict, one_family: bool):
    if one_family is False or None:
        if set(system.families.keys()).intersection(annot2fam.keys()) == system.families.keys():
            return True
        else:
            return False
    else:
        for annot in system.families.keys():
            return len(annot2fam.get(annot))


# TODO pour vérifier si annot in sys, faire un break si y'a pas tl mandatory
def present_system(system: System, annot2fam: dict):
    """
    Look if system is present in the pangenome
    :param system:
    :param annot2fam:
    :return:
    """
    set_bool = set()
    for annot in system.families:
        if annot in annot2fam.keys():
            set_bool.add(True)
        else:
            set_bool.add(False)
    if all(boolean for boolean in set_bool):
        draw_graph(system, annot2fam)


def draw_graph(system: System, annot2fam: dict):
    g = nx.Graph()
    for annot in system.families:
        for fam_obj in annot2fam[annot]:
            g.add_node(fam_obj)
    for annot in system.families:
        for fam_obj in annot2fam[annot]:
            if fam_obj in g.nodes():
                add_node_edge(g, g.nodes(), get_max_separation_sys(system)[annot], i=0)
    nx.draw(g)
    plt.show()
    return g


def add_node_edge(g: nx.Graph, set_fam: set, max_sep: int, i: int):
    pass
#     if i != max_sep+1:
#         i += 1
#         for family in set_fam:
#             for edge in family.edges:
#                 if edge in g.node():
#                     g.add_edge(fam_obj, edge.target)
            # if family.target in g.nodes():
    #
    #         else:
    #             edge.add(edge_fam)
    #             g.add_node(edge.target)
    #             g.add_edge(fam_obj, edge.target)
    #             nx.draw(g)
    #             plt.show()
    #     i += 1
    #     add_node_edge(g, edge.target.edges, edge.target, max_sep, i)
    pass


def get_annot_annot2fam(annot2fam: dict, max_sep_sys: dict):
    annot_sys_pan = dict()
    for annot in max_sep_sys.keys():
        if annot in annot2fam.keys():
            annot_sys_pan[annot] = annot2fam.get(annot)
    return annot_sys_pan


def get_max_separation(obj: object):
    """
    Get the parameter max_separation in a object

    :param obj:
    :return:
    """
    return obj.parameters["max_separation"]


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    systems_to_path = args.systems.absolute()
    system = System()
    for file in systems_to_path.glob("*.json"):
        system = read_systems(file)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for k, v in pan_to_path.items():
        pangenome = Pangenome(name=k)
        pangenome.add_file(v)
        annot2fam = annot_pangenome(pangenome, args.hmm, args.systems, args.force, args.disable_prog_bar)
        present_system(system, annot2fam)
        print("finish")


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
    required.add_argument('-s', '--systems',
                          required=True,
                          type=Path,
                          help="Path to systems directory")
    required.add_argument('-p', '--pangenomes',
                          required=True,
                          type=Path,
                          nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument('--hmm',
                          required=True,
                          type=Path,
                          nargs='?',
                          help='Findings gene family of hmm in .tsv file')

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--source", required=False, type=str, nargs="?", default=None,
                          help='Name of the annotation source. Default use name of hmm result file')


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
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
