#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import argparse
from typing import Union

# installed libraries
import ppanggolin.metadata

# local libraries
from panorama.annotate.hmm_search import profile_gfs
from panorama.format.read_binaries import check_pangenome_info, load_pangenomes
from panorama.format.write_proksee import write_proksee
from panorama.utils import check_tsv_sanity
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenomes
from panorama.systems.system_association import *
from panorama.format.conserved_spot import *
from panorama.alignment.align import all_against_all

# For Quentin
from panorama.format.figure import *
from multiprocessing import Manager
import tempfile

need_annotations = False
need_families = False
need_graph = False
need_partitions = False
need_spots = False
need_rgp = False
need_modules = False
need_gene_sequences = False
need_metadata = False
need_systems = False
bool_rgp = False
bool_modules = False
bool_spots = False


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    print("writing systems...")
    # check_flat_parameters(args)
    # pan_to_path = check_tsv_sanity(args.pangenomes)
    # mkdir(args.output, force=args.force)
    # write_flat_files(pan_to_path, args.pangenomes, output=args.output, annotation=args.annotations,
    #                  systems=args.systems,
    #                  systems_asso=args.systems_asso, conserved_spot=args.conserved_spot, draw_spot=args.draw_spot,
    #                  models=args.models, sources=args.sources, proksee=args.proksee,
    #                  proksee_template=args.proksee_template,
    #                  organisms_list=args.organisms, hmm=args.hmm, msa_path=args.msa, msa_format=args.msa_format,
    #                  threads=args.threads, force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write_systems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    required.add_argument('-m', '--models', required=False, type=Path, default=None, nargs="+",
                          help="Path to model list file.")
    required.add_argument("-s", "--sources", required=False, type=str, nargs="+", default=None,
                          help='Name of the annotation source where panorama as to select in pangenomes')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("-a", "--association", required=False, type=str, default=None,
                          choices=["all", "rgp-modules", "rgp-spots", "modules-spots", "modules", "rgp"],
                          help="Write association between systems and others pangenomes elements")
    optional.add_argument("--proksee", required=False, type=str, default=None, nargs='+',
                          choices=["all", "base", "modules", "rgp", "spots", "annotations", "systems"])
    optional.add_argument("--proksee_template", required=False, type=Path, default=None, nargs='?')
    optional.add_argument("--organisms", required=False, type=str, default=None, nargs='+')
    optional.add_argument("--threads", required=False, type=int, default=1)
