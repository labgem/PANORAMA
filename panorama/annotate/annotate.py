#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from typing import Dict, Tuple
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Manager, Lock

# installed libraries
from tqdm import tqdm
import pandas as pd
from ppanggolin.meta.meta import check_metadata_format, assign_metadata

# local libraries
from panorama.utils import init_lock
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.format.read_binaries import load_pangenomes
from panorama.annotate.hmm_search import read_hmms, annot_with_hmm
from panorama.pangenomes import Pangenome, Pangenomes


def check_parameter(args):
    """
    Checks the provided arguments to ensure that they are valid.

    Args:
        args: The parsed arguments.

    Raises:
        argparse.ArgumentError: If any required arguments are missing or invalid.
    """
    if args.table is None and args.hmms is None:
        raise argparse.ArgumentError(argument=None, message="Please provide either a table with annotation for gene "
                                                            "families or a hmms to annotate them.")
    if args.table is not None:
        if args.mode is not None:
            logging.getLogger('PANORAMA').error("You cannot specify both --table and --mode in the same command")
            raise argparse.ArgumentError(argument=None, message="--table is incompatible option with '--mode'.")
        if args.k_best_hit is not None:
            logging.getLogger('PANORAMA').error("You cannot specify both --table and --k_best_hit in the same command")
            raise argparse.ArgumentError(argument=None, message="--table is Incompatible option with '--k_best_hit'.")
        if args.only_best_hit:
            logging.getLogger('PANORAMA').error("You cannot specify both --table and --only_best_hit in the same command")
            raise argparse.ArgumentError(argument=None, message="--table is Incompatible option with '--only_best_hit'.")

    else:  # args.hmm is not None
        if args.mode is None:
            args.mode = "fast"
        if args.k_best_hit is not None:
            if args.k_best_hit < 1:
                raise argparse.ArgumentTypeError("k_best_hit must be greater than 1.")
            if args.only_best_hit:
                if args.k_best_hit == 1:
                    logging.getLogger("PANORAMA").warning("'--only_best_hit' is an alias for '--k_best_hit 1'. "
                                                          "You can use only one of them.")
                else:
                    logging.getLogger("PANORAMA").error(f"You set the maximum of hit at {args.k_best_hit} but "
                                                        f"you also use '--only_best_hit' which is an alias for '--k_best_hit 1'. "
                                                        f"Please use only one of those.")
                    raise argparse.ArgumentError(argument=None, message="--k_best_hit is incompatible with --only_best_hit")
        else:
            if args.only_best_hit:
                args.k_best_hit = 1

def check_pangenome_annotation(pangenome: Pangenome, source: str, force: bool = False):
    """
    Check pangenome information before adding annotation

    Args:
        pangenome: pangenome object that will be checked
        source: source of annotation to check if already in pangenome
        force: Flag to allow overwriting/erasing annotation
    Raises:
        KeyError: if a source with the same name already exists and force is False
    """
    if pangenome.status["metadata"]["families"] == "inFile" and source in pangenome.status["metasources"]["families"]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source)
        else:
            raise KeyError(f"A metadata corresponding to the source : '{source}' already exist in pangenome."
                            f"Add the option --force to erase")


def read_families_metadata(pangenome: Pangenome, metadata: Path) -> Tuple[pd.DataFrame, str]:
    """
    Read families metadata for one pangenome

    Args:
        pangenome: pangenome object for which metadata will be associated
        metadata:  Path to metadata file containing metadata to add to pangenome

    Returns:
        Tuple[pd.DataFrame, str]: the metadata dataframe and the name of the pangenome
    """
    metadata_df = check_metadata_format(metadata, "families")
    return metadata_df, pangenome.name


def read_families_metadata_mp(pangenomes: Pangenomes, table: Path, threads: int = 1,
                              lock: Lock = None, disable_bar: bool = False) -> Dict[str, pd.DataFrame]:
    """
    Read families metadata for multiple pangenomes in multiprocessing

    Args:
        pangenomes: Pangenomes object containing all the pangenome to annotate
        table: Path to metadata file for gene families
        threads: Number of available threads
        lock: Lock for multiprocessing execution
        disable_bar: Flag to disable progress bar

    Returns:
        Dict[str, pd.DataFrame]: Dictionary with metadata link to pangenome by its name
    """
    path_to_metadata = pd.read_csv(table, delimiter="\t", names=["Pangenome", "path"])
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome in pangenomes:
                logging.getLogger("PANORAMA").debug(f"read metadata for pangenome {pangenome.name}")
                metadata_file = path_to_metadata.loc[path_to_metadata["Pangenome"] == pangenome.name]["path"].squeeze()
                future = executor.submit(read_families_metadata, pangenome, metadata_file)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            results = {}
            for future in futures:
                res = future.result()
                results[res[1]] = res[0]
    return results


def keep_best_hit(metadata: pd.DataFrame, k_best_hit: int) -> pd.DataFrame:
    """
    Keep the k best hit for a given metadata

    Args:
        metadata: metadata dataframe with multiple annotation for gene families
        k_best_hit: number of best hits to keep

    Returns:
        pd.DataFrame: Filtered metadata dataframe with only the k best hit
    """
    get_best_hit = lambda group: group.nlargest(k_best_hit, columns=['score', 'e_value', 'bias'])
    return metadata.groupby('families', group_keys=False).apply(get_best_hit)



def write_annotations_to_pangenome(pangenome: Pangenome, metadata: pd.DataFrame, source: str, k_best_hit: int = None,
                                   force: bool = False, disable_bar: bool = False):
    """
    Write gene families annotation for one pangenome

    Args:
        pangenome: Pangenome link to metadata
        metadata: Metadata dataframe
        source: Metadata source
        k_best_hit: number of best hits to keep
        force: Boolean to allow force write in pangenomes
        disable_bar: Allow to disable progress bar
    """
    meta_df = keep_best_hit(metadata, k_best_hit) if k_best_hit is not None else metadata
    assign_metadata(meta_df, pangenome, source, "families", omit=True, disable_bar=disable_bar)
    write_pangenome(pangenome, pangenome.file, force=force, disable_bar=disable_bar)


def write_annotations_to_pangenomes(pangenomes: Pangenomes, pangenomes2metadata: Dict[str, pd.DataFrame],
                                    source: str, k_best_hit: int = None, threads: int = 1, lock: Lock = None,
                                    force: bool = False, disable_bar: bool = False):
    """
    Write gene families annotation for pangenomes in multiple processing

    Args:
        pangenomes:  Pangenomes object containing all the pangenome to annotate
        pangenomes2metadata: Dictionary with for each pangenomes the metadata dataframe associated
        source: Metadata source
        k_best_hit: number of best hits to keep
        threads: Number of available threads
        lock: Lock for multiprocessing execution
        force: Boolean to allow force write in pangenomes
        disable_bar: Allow to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome_name, metadata in pangenomes2metadata.items():
                pangenome = pangenomes.get_pangenome(pangenome_name)
                logging.getLogger("PANORAMA").debug(f"Write annotation for pangenome {pangenome.name}")
                future = executor.submit(write_annotations_to_pangenome, pangenome,
                                         metadata, source, k_best_hit, force, disable_bar)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def annot_pangenomes_with_hmm(pangenomes: Pangenomes, hmm: Path = None, mode: str = "fast", bit_cutoffs: str = None,
                              threads: int = 1, disable_bar: bool = False) -> Dict[str, pd.DataFrame]:
    """
    Main function to add annotation to pangenome from tsv file

    Args:
        pangenomes: Pangenomes object containing all the pangenome to annotate
        hmm: Path to hmm list file
        mode: Which mode to use to annotate gene families with HMM
        bit_cutoffs: Bit cutoff threshold mode for hmmsearch
        threads: Number of available threads
        disable_bar: Flag to disable progress bar

    Returns:
        Dict[str, pd.DataFrame]: Dictionary with for each pangenome a dataframe containing families metadata given by HMM
    """
    logging.getLogger("PANORAMA").info("Begin HMM searching")
    # Get list of HMM with Plan7 data model
    pangenome2annot = {}
    hmms, hmm_df = read_hmms(hmm, disable_bar=disable_bar)
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Align gene families to HMM for {pangenome.name}")
        pangenome2annot[pangenome.name] = annot_with_hmm(pangenome, hmms, hmm_df, mode, bit_cutoffs,
                                                         threads, disable_bar)

    return pangenome2annot


def annot_pangenomes(pangenomes: Pangenomes, source: str = None, table: Path = None,
                     hmm: Path = None, mode: str = "fast", k_best_hit: int = None, bit_cutoffs: str = None,
                     threads: int = 1, lock: Lock = None, force: bool = False, disable_bar: bool = False):
    """
    Gene families annotation with HMM or TSV files for multiple pangenomes in multiprocessing

    Args:
        pangenomes: Pangenomes object containing all the pangenome to annotate
        source: Name of the annotation source
        table: Path to metadata file for gene families annotation
        hmm: Path to hmm list file
        mode: Which mode to use to annotate gene families with HMM
        k_best_hit: number of best hits to keep
        bit_cutoffs: Bit cutoff threshold mode for hmmsearch
        threads: Number of available threads
        lock: Lock for multiprocessing
        force: Flag to allow force overwrite in pangenomes
        disable_bar: Flag to disable progress bar

    Raises:
        AssertionError: If neither HMM nor TSV are provided
    """
    assert table is not None or hmm is not None, 'Must provide either table or hmm'
    if table is not None:
        pangenomes2metadata = read_families_metadata_mp(pangenomes, table, threads, disable_bar)
    else:  # hmm is not None:
        pangenomes2metadata = annot_pangenomes_with_hmm(pangenomes, hmm, mode, bit_cutoffs, threads, disable_bar)
    write_annotations_to_pangenomes(pangenomes, pangenomes2metadata, source, k_best_hit, threads, lock, force, disable_bar)


def launch(args: argparse.Namespace) -> None:
    """
    Launch functions to annotate pangenomes

    Args:
        args: argument given in CLI
    """
    manager = Manager()
    lock = manager.Lock()
    check_parameter(args)

    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info={"need_annotations": True,
                                                                            "need_families": True},
                                 check_function=check_pangenome_annotation, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, source=args.source, force=args.force)

    annot_pangenomes(pangenomes=pangenomes, source=args.source, table=args.table,
                     hmm=args.hmm, mode=args.mode, k_best_hit=args.k_best_hit, bit_cutoffs=args.bit_cutoffs,
                     threads=args.threads, lock=lock, force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for annot command

    Returns:
        argparse.ArgumentParser: parser arguments for annot command
    """
    parser = sub_parser.add_parser("annotation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_annot(parser)
    return parser


def parser_annot(parser):
    """
    Add argument to parser for annot command

    Args:
        parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome.h5 files in .tsv file')
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source.')
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--table', type=Path, nargs='+', default=None,
                                help='A list of tab-separated file, containing annotation of gene families.'
                                     'Expected format is pangenome name in first column '
                                     'and path to the TSV with annotation in second column.')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?', default=None,
                                help="A tab-separated file with HMM information and path."
                                     "Note: Use panorama utils --hmm to create the HMM list file")
    hmm_param = parser.add_argument_group(title="HMM arguments",
                                          description="All of the following arguments are required,"
                                                      " if you're using HMM mode :")
    hmm_param.add_argument("--mode", required=False, type=str, default=None, choices=['fast', 'profile'],
                           help="Choose the mode use to align HMM database and gene families. "
                                "Fast will align the reference sequence of gene family against HMM."
                                "Profile will create an HMM profile for each gene family and "
                                "this profile will be aligned")
    hmm_param.add_argument("--k_best_hit", required=False, type=int, default=None,
                           help="Keep the k best annotation hit per gene family."
                                "If not specified, all hit will be kept.")
    hmm_param.add_argument("-b", "--only_best_hit", required=False, action="store_true",
                           help="alias to keep only the best hit for each gene family.")
    hmm_param.add_argument("--bit_cutoffs", required=False, type=str, default=None,
                           choices=["noise", "gathering", "trusted"],
                           help='Define the model-specific thresholding option to use for pyHMMer. '
                                'Will be automatically set by pyHMMer if not specified. '
                                'Look at pyHMMer and HMMer documentation for more information'
                           )
    # hmm_param.add_argument("--msa", required=False, type=Path, default=None,
    #                        help=argparse.SUPPRESS)
    # # help="To create a HMM profile for families, you can give a msa of each gene in families."
    # #      "This msa could be gotten from ppanggolin (See ppanggolin msa). "
    # #      "If no msa provide Panorama will launch one.")
    # hmm_param.add_argument("--msa-format", required=False, type=str, default="afa",
    #                        choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
    #                                 "clustal", "clustallike", "phylip", "phylips"],
    #                        help=argparse.SUPPRESS)
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads.")
