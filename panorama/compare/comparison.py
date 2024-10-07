#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import sys
import argparse
from pathlib import Path
import tempfile
from multiprocessing import Manager

# installed libraries

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.compare.context import context_comparison, context_comparison_parser


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    # check_parameter(args)

    pan_to_path = check_tsv_sanity(args.pangenomes)
    manager = Manager()
    lock = manager.Lock()
    if args.context:

        run_context_args = {"sequences": args.sequences, "families": args.family, "transitive": args.transitive,
                            "identity": args.identity, "coverage": args.coverage, "jaccard_threshold": args.jaccard,
                            "window_size": args.window_size,
                            "no_defrag": args.no_defrag, 'graph_format': args.graph_format}

        context_comparison(pan_to_path, args.context_results, args.family_clusters,
                           ppanggolin_context_args=run_context_args,
                           synteny_score=args.synteny_score, lock=lock, output=args.output,
                           tmpdir=args.tmpdir, task=args.task, threads_per_task=args.threads_per_task,
                           disable_bar=args.disable_prog_bar, force=args.force)

    elif args.compare_subcommand == "module":
        raise NotImplementedError(f"Module comparison is not yet implemented.")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("compare", description='Comparison of modules and gene contexts among pangenomes')

    parser_comparison(parser)

    if '--context' in sys.argv or '-h' in sys.argv or '--help' in sys.argv:
        context_comparison_parser(parser)
    #
    # module_parser = compare_subparser.add_parser("module")
    # parser_comparison(module_parser)
    # module_comparison_parser(module_parser)
    return parser


def parser_comparison(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')

    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")
    comparison = parser.add_argument_group(title="Pangenome comparison set object",
                                           description="Choose on which pangenome object you want to compare pangenomes")
    comparison.add_argument("--context", action='store_true',
                            help="Launch context comparison")
    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--task", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    optional.add_argument("--threads_per_task", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
