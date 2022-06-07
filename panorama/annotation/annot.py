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
from ppanggolin.formats.readBinaries import check_pangenome_info

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.annotation.rules import System, Systems
from panorama.pangenomes import Pangenome


def check_pangenome_annotation(pangenome: Pangenome, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome:
    :param disable_bar:
    :return:
    """
    check_pangenome_info(pangenome, need_families=True, disable_bar=disable_bar)


def read_files(systems_path, systems=Systems()):
    """
    Read all json files in the directory

    :param systems_path: path of systems directory
    :param systems: class Systems with all systems
    """
    for file in systems_path.glob("*.json"):
        with open(file.resolve().as_posix()) as json_file:
            data = json.load(json_file)
            system = System()
            system.read_system(data)
            systems.add_sys(system)
    systems.print_systems()


def parser_pangenomes(pangenomes_path):
    """
    Read all pangenome in a list tsv file

    :param pangenomes_path: Path to the list tsv file
    """
    try:
        file = open(pangenomes_path.absolute(), 'r')
        pangenomes = csv.reader(file, delimiter="\t")
    except IOError:
        raise IOError
    except Exception:
        raise Exception("Unexpected error when opening the list of pangenomes")
    else:
        for line in pangenomes:
            print(line[1])
            pangenome_path = Path(line[1])


def parser_hmmer(hmmer_path: Path) -> dict:
    """
    Recover information of hmmer results in a tsv file

    :param hmmer_path: Path to hmmer results of tsv file

    :raise:

    :return:
    """
    hmm_res = {}
    try:
        hmm_file = open(hmmer_path.absolute(), 'r')
    except IOError:
        raise IOError
    except Exception:
        raise Exception("Unexpected error when opening the list of hmmer")
    else:
        hmmer = csv.reader(hmm_file, delimiter="\t")
        for line in hmmer:
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
                hmm_res[query_name] = {"name": {target_name},
                                       "accession": {accession},
                                       "e-value": {evalue},
                                       "description": {system}
                                       }
        return hmm_res


def annot_pangenome(pangenome: Pangenome, hmm: Path, system: Path, source: str = None, force: bool = False,
                    disable_bar: bool = False):
    """ Main function to add annotation to pangenome

    :param pangenome:
    :param hmm:
    :param system:
    :param source:
    :param force:
    :param disable_bar:
    """
    check_pangenome_annotation(pangenome, disable_bar=disable_bar)

    for gene_fam_name, hmm_res in parser_hmmer(hmm).items():
        gene_fam = pangenome.get_gene_family(name=gene_fam_name)
        gene_fam.add_annotation(source=hmm.stem if source is None else source, annotation=hmm_res, force=force)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for k, v in pan_to_path.items():
        pangenome = Pangenome(name=k)
        pangenome.add_file(v)
        annot_pangenome(pangenome, args.hmm, args.systems, args.force, args.disable_prog_bar)


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
                          help='A list of pangenome .h5 files')
    required.add_argument('--hmm',
                          required=True,
                          type=Path,
                          nargs='?',
                          help='Findings systems of hmmer in tsv file')

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
