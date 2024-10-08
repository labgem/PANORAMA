#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from pathlib import Path
import tempfile
from typing import Dict, Union
from multiprocessing import Manager, Lock
import subprocess
from time import time
from shutil import rmtree

# installed libraries
from tqdm import tqdm
import pandas as pd
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import search_gene_context_in_pangenome, check_pangenome_for_context_search

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.region import GeneContext
from panorama.alignment.align import parser_mmseqs2_align
from panorama.compare.utils import parser_comparison, common_launch


def check_compare_systems_args(args):
    pass


def check_pangenome_for_systems_comparison():
    pass


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    print("hello")
    # check_compare_systems_args(args)
    # pangenomes, tmpdir, _, lock = common_launch(args, check_pangenome_for_systems_comparison,
    #                                             {"need_families": True,
    #                                              'need_annotations': False if args.context_results else True})
    #
    # output = mkdir(args.output, force=args.force)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("compare_systems",
                                   description='Comparison of systems among pangenomes')

    parser_comparison_context(parser)
    return parser


def parser_comparison_context(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required, compare_opt, optional = parser_comparison(parser)
