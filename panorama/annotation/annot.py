from __future__ import annotations

import argparse
import csv
import json
import os
import re

from pathlib import Path
from rules import System, Systems


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
        pangenomes = csv.reader(file, delimiter=" ")
    except IOError:
        raise IOError
    except Exception:
        raise Exception("Unexpected error when opening the list of pangenomes")
    else:
        for line in pangenomes:
            print(line[2])
            pangenome_path = Path(line[2])


def parser_hmmer(hmmer_path):
    """
    Recover information of hmmer results in a tsv file

    :param hmmer_path: Path to hmmer results of tsv file
    """
    try:
        file = open(hmmer_path.absolute(), 'r')
        hmmer = csv.reader(file, delimiter="\t")
    except IOError:
        raise IOError
    except Exception:
        raise Exception("Unexpected error when opening the list of hmmer")
    else:
        for line in hmmer:
            if not(re.compile("^#").search(line[0])):
                target_name = line[0]
                if not(line[1] == "-"):
                    accession = line[1]
                else:
                    accession = None
                query_name = line[2]
                evalue = line[4]
                if not(line[18] == "-"):
                    system = line[18]
                else:
                    system = target_name
                print(f'name : {target_name}, accession : {accession}, query : {query_name}, e-value : {evalue}, '
                      f'system : {system}')


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    # systems_path = Path(args.systems)
    # read_files(systems_path)
    list_pangenomes_path = Path(args.pangenomes)
    parser_pangenomes(list_pangenomes_path)
    sys_hmm_path = Path(args.hmmer)
    parser_hmmer(sys_hmm_path)


def parser_annot(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-s', '--systems',
                          required=True,
                          type=str,
                          help="Path to systems directory")
    required.add_argument('-p', '--pangenomes',
                          required=True,
                          type=str,
                          help='A list of pangenome .h5 files')
    required.add_argument('-hm', '--hmmer',
                          required=True,
                          type=str,
                          help='Findings systems of hmmer in tsv file')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser_annot(parser)
    args = parser.parse_args()
    launch(args)
