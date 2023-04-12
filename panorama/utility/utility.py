#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from pathlib import Path
import tempfile

import panorama
# install libraries

# local libraries
from panorama.utils import mkdir
from panorama.utility.translate import launch_translate
from panorama.models import Models


def check_parameters(args):
    if args.translate is None and args.models is None:
        raise argparse.ArgumentError(argument=None, message="You did not provide any action to do")
    if args.translate is not None:
        if args.output is None:
            raise argparse.ArgumentError(argument=args.output,
                                         message="Required to give an output directory to translate")
        if args.source == "padloc" and args.meta is None:
            raise argparse.ArgumentError(argument=args.meta,
                                         message="If you're using padloc as source, "
                                                 "it's necessary to give metadata file associate")
        if args.source == "defense-finder" and args.hmm is None:
            raise argparse.ArgumentError(argument=args.hmm,
                                         message="If you're using defense-finder as source, "
                                                 "it's necessary to give hmm path to get hidden information")


def check_models(models_path: Path, disable_bar: bool = False) -> Models:
    """Read to check all json files models in the directory

    :param models_path: path of models directory
    :param disable_bar: Disable progress bar

    :raise KeyError: One or more keys are missing or non-acceptable
    :raise TypeError: One or more value are not with good presence
    :raise ValueError: One or more value are not non-acceptable
    :raise Exception: Manage unexpected error
    """
    try:
        logging.getLogger().info("Check models translation...")
        models = Models()
        models.read(models_path, disable_bar)
    except Exception:
        raise Exception("Problem with translated models. Check that you give correct input and option. "
                        "If nothing wrong please report an issue on our github.")


def launch(args):
    """
    Launch utilities function for pangenomes

    :param args: Argument given
    """
    check_parameters(args)
    if args.translate is not None:
        outdir = mkdir(output=args.output, force=args.force)
        launch_translate(models=args.translate, source=args.source, output=outdir, hmms_path=args.hmm,
                         meta_data=args.meta, tmpdir=args.tmpdir, disable_bar=args.disable_prog_bar)
        check_models(models_path=outdir, disable_bar=args.disable_prog_bar)

    if args.models is not None:
        check_models(models_path=args.models, disable_bar=args.disable_prog_bar)


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
    translate = parser.add_argument_group(title="Translate arguments",
                                          description="Arguments to translate systems models from different sources")
    translate.add_argument("--translate", required=False, type=Path, default=None,
                           help="path to models to be translated")
    translate.add_argument("--source", required=False, type=str, choices=["padloc", "defense-finder"], nargs='?',
                           help="Available sources that we know how to translate. "
                                "The directory will be read recursively to catch all models.")
    translate.add_argument("--meta", required=False, type=Path,
                           help='For Padloc we use the metadata file to translate models. '
                                'A padloc_meta.tsv will be generated for PANORAMA to annotate in the output directory')
    translate.add_argument("--hmm", required=False, type=Path,
                           help='For Defense finder we need HMM to get all hidden information to translate models.')
    check_model = parser.add_argument_group(title="Check models arguments",
                                            description="Arguments necessary to check if models are readable")
    check_model.add_argument("--models", required=False, type=Path, nargs="?", default=None,
                             help="Path to models directory")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=Path, nargs='?', default=None,
                          help='Path to output directory.')
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")


if __name__ == "__main__":
    from panorama.utils import set_verbosity_level, check_log

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
