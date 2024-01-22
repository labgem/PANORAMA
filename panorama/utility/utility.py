#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from typing import List
from pathlib import Path
import tempfile

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from panorama.utils import mkdir
from panorama.utility.genInput import create_hmm_list_file, read_metadata
from panorama.utility.translate import launch_translate
from panorama.models import Models


def check_parameters(args: argparse.Namespace) -> None:
    """
    Checks the provided arguments to ensure that they are valid.

    Args:
        args: The parsed arguments.

    Raises:
        argparse.ArgumentError: If any required arguments are missing or invalid.
    """
    if not any(arg is not None for arg in [args.hmm, args.translate, args.models]):
        raise argparse.ArgumentError(argument=None, message="You did not provide any action to do")
    if args.hmm is not None and args.output is None:
        raise argparse.ArgumentError(argument=args.output, message="Required to give an output directory to translate")
    if args.models is not None:
        if args.output is None:
            raise argparse.ArgumentError(argument=args.output,
                                         message="Required to give an output directory to translate")
        if args.meta is not None:
            logging.getLogger("PANORAMA").warning("Metadata are not supported for models."
                                                  "Report an issue if you want it in future version")
    if args.translate is not None:
        if args.output is None:
            raise argparse.ArgumentError(argument=args.output,
                                         message="Required to give an output directory to translate")
        if args.source is None:
            raise argparse.ArgumentError(argument=args.source, message="Required to know how to translate models")


def create_models_list(models_path: List[Path], output:Path, recursive: bool = False,
                       disable_bar: bool = False) -> None:
    """
    Create a file that listing models and path to them. Also, models are checked

    Args:
        models_path: List of paths to models
        output: Directory to write models list file
        recursive: Flag to read models directory recursively (default: False)
        disable_bar: Flag to disable progress bar (default: False)

    Returns:
        None
    """
    logging.getLogger("PANORAMA").info("Begin to create model list file...")
    model_list = []
    models = Models()
    for path in tqdm(models_path, unit="models", disable=disable_bar):
        if path.is_file():
            models.read(path)
            model_list.append([path.stem, path])
        elif path.is_dir():
            for model_file in path.rglob("*.hmm") if recursive else path.glob("*.hmm"):
                models.read(model_file)
                model_list.append([model_file.stem, model_file])
        else:
            if not path.exists():
                raise FileNotFoundError(f"The given path is not find: {path}")
            else:
                raise Exception("Unexpected error")
    model_df = pd.DataFrame(model_list, columns=['name', 'path'])
    model_df = model_df.sort_values('name')
    model_df.to_csv(output/"models_list.tsv", sep="\t", header=False, index=False)


def check_models(models_list: Path, disable_bar: bool = False):
    """
    Checks all JSON files in the given directory to ensure that they are valid models.

    Args:
        models_list (Path): paths to the models_list.tsv.
        disable_bar (bool, optional): Whether to disable the progress bar. Defaults to False.

    Raises:
        KeyError: If any keys are missing or non-acceptable.
        TypeError: If any values are not with good presence.
        ValueError: If any values are not non-acceptable.
        Exception: If an unexpected error occurs.
    """
    logging.info("Check models translation...")
    models_df = pd.read_csv(models_list, sep="\t", names=['name', 'path'])
    models = Models()
    for idx, model_file in tqdm(models_df["path"].items(), total=models_df.shape[0],
                                desc="Check models", unit="model", disable=disable_bar):
        try:
            models.read(Path(model_file))
        except Exception:
            raise Exception(f"Problem with model {model_file}. Check that you give correct input and option. "
                            "If nothing wrong please report an issue on our github.")


def launch(args: argparse.Namespace):
    """
    Launches the utilities function for pangenomes.

    Args:
        args: The parsed arguments.

    Returns:
        None
    """
    check_parameters(args)
    if args.hmm:
        outdir = mkdir(output=args.output, force=args.force)
        if args.meta:
            logging.getLogger("PANORAMA").info("Read metadata file...")
            metadata_df = read_metadata(args.meta)
        else:
            metadata_df = None
        create_hmm_list_file(args.hmm, outdir, metadata_df, args.recursive, args.disable_prog_bar)
    if args.translate is not None:
        outdir = mkdir(output=args.output, force=args.force)
        launch_translate(db=args.translate, source=args.source, output=outdir,
                         tmpdir=args.tmp, force=args.force, disable_bar=args.disable_prog_bar)
        check_models(models_list=outdir/"models_list.tsv", disable_bar=args.disable_prog_bar)

    if args.models is not None:
        outdir = mkdir(output=args.output, force=args.force)
        create_models_list(models_path=args.models, output=outdir, recursive=args.recursive, disable_bar=args.disable_prog_bar)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in the command line.

    Args:
        sub_parser (argparse._SubParsersAction): The subparser for the utils command.

    Returns:
        argparse.ArgumentParser: The parser for the utils command.
    """
    parser = sub_parser.add_parser("utils", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_utils(parser)
    return parser


def parser_utils(parser: argparse.ArgumentParser):
    """
    Parser for the specific arguments of the utils command.

    Args:
        parser (argparse.ArgumentParser): The parser for the utils command.

    Returns:
        None
    """
    panorama_input = parser.add_argument_group("Create input files arguments",
                                               description="Create some input files used by PANORAMA")
    panorama_input.add_argument("--hmm", required=False, type=Path, default=None, nargs='+',
                                help="Path to HMM files or directory containing HMM")
    panorama_input.add_argument("--models", required=False, type=Path, nargs="+", default=None,
                                help="Create a models_list.tsv file from the given models and check them.")
    panorama_input.add_argument("--meta", required=False, type=Path, default=None, nargs='?',
                                help="Path to metadata file to add some to list file")
    panorama_input.add_argument("--recursive", required=False, action="store_true", default=False,
                                help="Flag to indicate if directories should be read recursively")
    translate = parser.add_argument_group(title="Translate arguments",
                                          description="Arguments to translate systems models from different sources")
    translate.add_argument("--translate", required=False, type=Path, default=None,
                           help="Path to models to be translated. Give the directory with models, hmms and other files."
                                "PANORAMA will take care of everything it needs to translate.")
    translate.add_argument("--source", required=False, type=str, choices=["padloc", "defense-finder", "macsy-finder"],
                           nargs='?',
                           help="Available sources that we know how to translate. "
                                "The directory will be read recursively to catch all models.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('-o', '--output', required=False, type=Path, nargs='?', default=None,
                          help='Path to output directory.')
    optional.add_argument("--tmp", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
