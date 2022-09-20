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
import re

from ppanggolin.genome import Organism
from tqdm import tqdm
# installed libraries
import pandas as pd
from ppanggolin.formats.readBinaries import check_pangenome_info

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.annotation.rules import System, Systems
from panorama.pangenomes import Pangenome
from panorama.annotation.hmm_search import annot_with_hmm, res_col_names
from panorama.annotation.system_search import launch_system_search


def check_parameter(args):
    if args.tsv is None and args.hmm is None:
        raise Exception("You did not provide tsv or hmm for annotation")
    if args.tsv is not None and args.e_value is None:
        raise Exception("An e-value is needed to keep the best annotation for each family")
    if args.hmm is not None and args.e_value is None and args.cutoffs is None:
        raise Exception("You did not provide e_value or cutoffs file to assign function with HMM")


def check_pangenome_annotation(pangenome: Pangenome, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param disable_bar: Disable bar
    """
    check_pangenome_info(pangenome, need_families=True, need_graph=True, need_annotations=True, disable_bar=disable_bar)


def read_systems(systems_path: Path, systems=Systems()):
    """ Read all json files systems in the directory

    :param systems_path: path of systems directory
    :param systems: class Systems with all systems
    """
    for file in systems_path.glob("*.json"):
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


def check_presence_family(system: System, annot2fam: dict, one_family: bool):
    if one_family is False or None:
        if set(system.families.keys()).intersection(annot2fam.keys()) == system.families.keys():
            return True
        else:
            return False
    else:
        for annot in system.families.keys():
            return len(annot2fam.get(annot))


def search_system(systems: Systems, annot2fam: dict, disable_bar: bool = False):
    """
    Search present system in the pangenome
    :param systems:
    :param annot2fam:
    :return:
    """
    pred = []
    org_pred = {}
    spot_pred = {}

    def org2pred(org_pred_dict: dict, pred_res: dict, system_name: str):
        for keys, value in pred_res.items():
            org_inter = set()
            for fam in value:
                if len(org_inter) == 0:
                    org_inter = fam.organisms
                else:
                    org_inter.intersection(fam.organisms)
            if len(org_inter) == 0:
                print(system_name)
            for org in org_inter:
                all_fam = {family: False for family in value}
                index = 0
                for contig in org.contigs:
                    for gene in contig.genes:
                        if gene.family in all_fam:
                            all_fam[gene.family] = True
                if all(x for x in all_fam.values()):
                    if org.name in org_pred_dict:
                        org_pred_dict[org.name].add(system_name)
                    else:
                        org_pred_dict[org.name] = {system_name}

    for system in tqdm(systems, total=systems.size, unit='system', disable=disable_bar):
        # if check_all_present_families(system, annot2fam) is True:
        pred_res = launch_system_search(system, annot2fam)
        if pred_res is not None:
            pred.append([system.name, max(pred_res.keys()) + 1])
            org2pred(org_pred, pred_res, system.name)
    proj = pd.DataFrame.from_dict(org_pred, orient='index')
    sys_df = pd.DataFrame(pred, columns=['System', 'Nb Detection']).sort_values('System').reset_index(drop=True)
    proj.to_csv("projection5.tsv", sep="\t", index=['Organisms'], index_label='Organisms',
                header=[f"System {i}" for i in range(1, proj.shape[1] + 1)])
    sys_df.to_csv("system5.tsv", sep="\t", header=['System', "Nb_detected"])


def filter_df(annotated_df: pd.DataFrame, eval_threshold: float = 0.0001) -> pd.DataFrame:
    """
    Filter annotation from dataframe to keep the best annotation for each family.

    :param annotated_df: Dataframe with annotation and families
    :param eval_threshold: e-value threshold

    :return: Filtered dataframe
    """
    df_filter = annotated_df[annotated_df[res_col_names[3]] <= eval_threshold]
    df_eval = df_filter[
        df_filter[res_col_names[3]] == df_filter.groupby(res_col_names[0])[res_col_names[3]].transform(min)]
    return df_eval[df_eval[res_col_names[4]] == df_eval.groupby(res_col_names[0])[res_col_names[4]].transform(max)]


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
    for gf in annotation_df[res_col_names[0]].unique():
        select_df = annotation_df.loc[annotation_df[res_col_names[0]] == gf]
        gene_fam = pangenome.get_gene_family(name=gf)
        gene_fam.add_annotation(source=source, annotation=list(select_df['protein_name']), force=force)
        for index, row in select_df.iterrows():
            if row['protein_name'] not in annot2fam:
                annot2fam[row['protein_name']] = {gene_fam}
            else:
                annot2fam[row['protein_name']].add(gene_fam)
            if not pd.isna(row['secondary_name']):  # Test if there is a second name
                for other_name in re.split('\|', row['secondary_name']):
                    if other_name not in annot2fam:
                        annot2fam[other_name] = {gene_fam}
                    else:
                        annot2fam[other_name].add(gene_fam)
    return annot2fam


def annot_pangenome(pangenome: Pangenome, hmm: Path, tsv: Path, meta: Path = None, e_value: float = 0.0001,
                    source: str = None, prediction_size: int = 1, tmpdir: Path = Path(tempfile.gettempdir()),
                    threads: int = 1, force: bool = False, disable_bar: bool = False) -> dict:
    """ Main function to add annotation to pangenome from tsv file

    :param pangenome: Pangenome object to ppanggolin
    :param hmm: Path to hmm file or directory
    :param tsv: Path to tsv with gene families annotation
    :param meta:
    :param e_value:
    :param source: Source of annotation
    :param prediction_size:
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param force: Boolean of force argument
    :param disable_bar: Disable bar

    :return: Dictionnary with for each annotation a set of corresponding gene families
    """
    check_pangenome_annotation(pangenome, disable_bar=disable_bar)
    if tsv is not None:
        annotation_df = pd.read_csv(tsv, sep="\t", header=None, quoting=csv.QUOTE_NONE, names=res_col_names)
        if e_value != 0:
            logging.getLogger().debug("Filter TSV annotation file")
            annotation_df = filter_df(annotation_df, e_value)
    elif hmm is not None:
        annotation_df = annot_with_hmm(pangenome, hmm, method="hmmsearch", meta=meta, prediction_size=prediction_size,
                                       tmpdir=tmpdir, threads=threads, disable_bar=disable_bar)
    else:
        raise Exception("You did not provide tsv or hmm for annotation")
    return annotation_to_families(annotation_df, pangenome, source, force)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    check_parameter(args)
    systems_to_path = args.systems.absolute()
    systems = read_systems(systems_to_path)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for pangenome_name, pangenome_file in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name)
        pangenome.add_file(pangenome_file)
        annot2fam = annot_pangenome(pangenome=pangenome, hmm=args.hmm, tsv=args.tsv, source=args.source,
                                    meta=args.meta, e_value=args.e_value, prediction_size=args.prediction_size,
                                    tmpdir=args.tmpdir, threads=args.threads, force=args.force,
                                    disable_bar=args.disable_prog_bar)
        search_system(systems, annot2fam, args.disable_prog_bar)
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
    required.add_argument("--meta", required=True, type=Path,
                          help="Metadata link to HMM with protein name, description and cutoff")
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--tsv', type=Path, nargs='?',
                                help='Gene families annotation in TSV file. See our github for more detail about format')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?',
                                help="File with all HMM or a directory with one HMM by file")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--prediction_size", required=False, type=int, default=1,
                          help="Number of prediction associate with gene families")
    optional.add_argument('--e_value', required=False, type=float, default=0,
                          help="Set the same e-value for all hmm to filter TSV")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of av available threads")


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
