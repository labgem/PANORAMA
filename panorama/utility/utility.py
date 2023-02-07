#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
from pathlib import Path
import tempfile

# install libraries

# local libraries
from panorama.utils import mkdir
from panorama.utility.translate import launch_translate


def check_parameters(args):
    if args.source == "padloc" and args.meta is None:
        raise argparse.ArgumentError(argument=None, message="If you're using padloc as source, "
                                                            "it's necessary to give metadata file associate")


def launch(args):
    """
    Launch utilities function for pangenomes

    :param args: Argument given
    """
    check_parameters(args)
    outdir = mkdir(output=args.output, force=args.force)
    if args.translate:
        launch_translate(models=args.models, source=args.source, output=outdir,
                         meta_data=args.meta, tmpdir=args.tmpdir)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("utility", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_utils(parser)
    return parser


def parser_utils(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-o', '--output', required=True, type=Path, nargs='?',
                          help='Path to output directory.')
    translate = parser.add_argument_group(title="Translate arguments",
                                          description="Arguments to translate systems models from different sources")
    translate.add_argument("--translate", required=False, action='store_true',
                           help="Flag to launch the transate function.")
    translate.add_argument("--source", required=False, type=str, choices=["padloc", "defense-finder"], nargs='?',
                           help="Available sources that we know how to translate. "
                                "The directory will be read recursively to catch all models.")
    translate.add_argument("--models", required=False, type=Path, nargs="?",
                           help="Path to models directory")
    translate.add_argument("--meta", required=False, type=Path,
                           help='For Padloc we use the metadata file to translate models.')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")


if __name__ == "__main__":
    from overall import set_verbosity_level, check_log
    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_utils(main_parser)
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
