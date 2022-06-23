#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import csv
import json
import logging
from pathlib import Path
import tempfile

# installed libraries
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from ppanggolin.formats.readBinaries import check_pangenome_info

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.annotation.rules import System, Systems
from panorama.pangenomes import Pangenome
from panorama.annotation.hmm_search import launch_hmm_search

col_names = ['Gene family', 'Annotation', 'Accesion Id', 'e-value', 'score', 'overlap', 'Annotation description']


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


def filter_df(annotated_df: pd.DataFrame, eval_threshold: float = 0.0001):
    """
    Filter annotation from dataframe to keep the best annotation for each family.

    :param annotated_df: Dataframe with annotation and families
    :param eval_threshold: e-value threshold

    :return: Filtered dataframe
    """
    df_filter = annotated_df[annotated_df[col_names[3]] <= eval_threshold]
    df_eval = df_filter[df_filter[col_names[3]] == df_filter.groupby(col_names[0])[col_names[3]].transform(min)]
    return df_eval[df_eval[col_names[4]] == df_eval.groupby(col_names[0])[col_names[4]].transform(max)]


def annotation_to_families(annotation_df: pd.DataFrame, pangenome: Pangenome, source: str = None,
                           force: bool = False) -> dict:
    """ Add to gene families an annotation and create a dictionary with for each annotation a set of gene family

    :param annotation_df: Dataframe with for each family an annotation
    :param pangenome: Pangenome with gene families
    :param source: source of the annotation
    :param force: force to write the annotation in gene families

    :return: Dictionary with for each annotation a set of gene family
    """
    annot2fam = {}
    for index, row in annotation_df.iterrows():
        gene_fam = pangenome.get_gene_family(name=row['Gene family'])
        gene_fam.add_annotation(source=source, annotation=row['Annotation'], force=force)
        if row['Annotation'] not in annot2fam:
            annot2fam[row['Annotation']] = {gene_fam}
        else:
            annot2fam[row['Annotation']].add(gene_fam)
    return annot2fam


def annot_pangenome(pangenome: Pangenome, hmm: Path, tsv: Path, e_value: float = 0.0001, source: str = None,
                    tmpdir: Path = Path(tempfile.gettempdir()), threads: int = 1, force: bool = False,
                    disable_bar: bool = False) -> dict:
    """ Main function to add annotation to pangenome from tsv file

    :param pangenome: Pangenome object to ppanggolin
    :param hmm: Path to hmm file or directory
    :param tsv: Path to tsv with gene families annotation
    :param e_value: e-value threshold to associate annoation to gene family
    :param source: Source of annotation
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param force: Boolean of force argument
    :param disable_bar: Disable bar

    :return: Dictionnary with for each annotation a set of corresponding gene families
    """
    check_pangenome_annotation(pangenome, disable_bar=disable_bar)
    if tsv is not None:
        annotation_df = pd.read_csv(tsv, sep="\t", header=None, quoting=csv.QUOTE_NONE,
                                    names=col_names)
    elif hmm is not None:
        annotation_df = launch_hmm_search(pangenome, hmm, tmpdir=tmpdir, threads=threads, disable_bar=disable_bar)
    else:
        raise Exception("You did not provide tsv or hmm for annotation")
    annotation_fitlered = filter_df(annotation_df, e_value)
    pd.set_option('display.max_columns', None)
    print(annotation_fitlered[annotation_fitlered.duplicated(subset=col_names[0])])
    pd.reset_option('max_columns')
    return annotation_to_families(annotation_fitlered, pangenome, source, force)


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
        annot2fam = annot_pangenome(pangenome=pangenome, hmm=args.hmm, tsv=args.tsv, source=args.source,
                                    e_value=args.e_value, tmpdir=args.tmpdir, threads=args.threads,
                                    force=args.force, disable_bar=args.disable_prog_bar)
        present_system(system, annot2fam)
        logging.getLogger().info("Annotation Done")


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
    required.add_argument('-s', '--systems', required=True, type=Path,
                          help="Path to systems directory")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    exclusive = required.add_mutually_exclusive_group(required=True)
    exclusive.add_argument('--tsv', type=Path, nargs='?',
                           help='Gene families annotation in TSV file. See our github for more detail about format')
    exclusive.add_argument('--hmm', type=Path, nargs='?',
                           help="File with all HMM or a directory with one HMM by file")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--e_value", required=False, type=float, nargs='?', default=0.0001,
                          help='E-value threshold use to assign annotation to gene families')
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of av available threads")


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
