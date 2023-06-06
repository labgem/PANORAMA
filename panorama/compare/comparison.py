#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import sys
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
        # run_context_args = {"sequences": args.sequences, "families": args.family, "transitive": args.transitive,
        #                 "identity": args.identity, "coverage": args.coverage, "jaccard": args.jaccard,
        #                 "no_defrag": args.no_defrag}
        run_context_args = {} # Not implemented yet
        context_comparison(pan_to_path, args.context_results, args.family_clusters, lock, args.output, args.tmpdir, args.task, args.threads_per_task,
                           args.disable_prog_bar, args.force, **run_context_args)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("compare", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_comparison(parser)


    if '--context' in sys.argv:
        # remove 'Pangenome comparison set object' groups
        parser._action_groups = [action_grp for action_grp in parser._action_groups if action_grp.title != 'Pangenome comparison set object']

        context_comparison_parser(parser)
        
    if "--modules" in sys.argv:
        # remove 'Pangenome comparison set object' groups
        parser._action_groups = [action_grp for action_grp in parser._action_groups if action_grp.title != 'Pangenome comparison set object']

        raise NotImplementedError("Functionnality still in progress")
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
    
    comparison.add_argument("--modules", action='store_true',
                            help="Launch modules comparison")
    

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--task", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    optional.add_argument("--threads_per_task", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_comparison(main_parser)
    if '--context' in sys.argv:
        context_comparison_parser(main_parser)
    common = main_parser._action_groups.pop(1)  # get the 'optional arguments' action group.
    common.title = "Common arguments"
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    main_parser._action_groups.append(common)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
