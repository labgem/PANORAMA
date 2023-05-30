#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
from pathlib import Path
import tempfile
from multiprocessing import Manager, Lock

# installed libraries

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.compare.context import context_comparison


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
        context_comparison(pan_to_path, lock, args.cpu, args.disable_prog_bar)


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
    context = parser.add_argument_group(title="Context comparison arguments")
    context.add_argument("--context", action='store_true',
                         help="Launch context comparison")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--cpu", required=False, nargs='?', type=int, default=1,
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
