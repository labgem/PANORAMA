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
from ppanggolin.meta.meta import check_metadata_format, assign_metadata, write_pangenome

# local libraries
from panorama.utils import init_lock
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.format.read_binaries import load_multiple_pangenomes
from panorama.annotate.hmm_search import read_metadata, read_hmm, annot_with_hmm
from panorama.pangenomes import Pangenome, Pangenomes


def check_pangenome_annotation(pangenome: Pangenome, source: str, force: bool = False):
    """ Check pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: source of the annotation
    :param force: erase if an annotation for the provide source already exist
    """
    if pangenome.status["metadata"]["families"] == "inFile" and source in pangenome.status["metasources"]["families"]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source)
        else:
            raise Exception(f"A metadata corresponding to the source : '{source}' already exist in pangenome."
                            f"Add the option --force to erase")


def read_families_metadata(pangenome: Pangenome, metadata: Path) -> Tuple[pd.DataFrame, str]:
    """Read families metadata for one pangenome

    :param pangenome: Pangenome link to metadata
    :param metadata: Metadata file

    :return: Dataframe with metadata link to pangenome name
    """
    metadata_df = check_metadata_format(metadata, "families")
    return metadata_df, pangenome.name


def read_families_metadata_mp(pangenomes: Pangenomes, metadata_list: Path, threads: int = 1,
                              lock: Lock = None, disable_bar: bool = False) -> Dict[str, pd.DataFrame]:
    """Read families metadata for multiple pangenomes in multiprocessing

    :param pangenomes: Pangenomes object containing all the pangenome to annotate
    :param metadata_list: Path to list file with for each pangenome an associated metadata file for gene families
    :param threads: Number of available threads
    :param lock: Lock for multiprocessing execution
    :param disable_bar: Disable bar

    :return: Dataframe with metadata link to pangenome name
    """
    path_to_metadata = pd.read_csv(metadata_list, delimiter="\t", names=["Pangenome", "path"])
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []

            for pangenome in pangenomes:
                logging.debug(f"read metadata for pangenome {pangenome.name}")
                metadata_file = path_to_metadata.loc[path_to_metadata["Pangenome"] == pangenome.name]["path"].squeeze()
                future = executor.submit(read_families_metadata, pangenome, metadata_file)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            results = {}
            for future in futures:
                res = future.result()
                results[res[1]] = res[0]
    return results


def write_annotations_to_pangenome(pangenome: Pangenome, metadata: pd.DataFrame, source: str,
                                   force: bool = False, disable_bar: bool = False):
    """Write gene families annotation for one pangenome

    :param pangenome: Pangenome link to metadata
    :param metadata: Metadata dataframe
    :param source: Metadata source
    :param force: Boolean to allow force write in pangenomes
    :param disable_bar: Allow to disable progress bar
    """
    assign_metadata(metadata, pangenome, source, "families", omit=True, disable_bar=disable_bar)
    write_pangenome(pangenome, pangenome.file, force=force, disable_bar=disable_bar)


def write_annotations_to_pangenomes(pangenomes: Pangenomes, pangenomes2metadata: Dict[str, pd.DataFrame],
                                    source: str, threads: int = 1, lock: Lock = None,
                                    force: bool = False, disable_bar: bool = False):
    """Write gene families annotation for pangenomes in multiple processing

    :param pangenomes: Pangenomes object containing all the pangenome to annotate
    :param pangenomes2metadata: Dictionnary with for each pangenomes the metadata dataframe associated
    :param source: Metadata source
    :param threads: Number of available threads
    :param lock: Lock for multiprocessing execution
    :param force: Boolean to allow force write in pangenomes
    :param disable_bar: Allow to disable progress bar
    :return:
    """
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []

            for pangenome_name, metadata in pangenomes2metadata.items():
                pangenome = pangenomes.get_pangenome(pangenome_name)
                logging.debug(f"Write annotation for pangenome {pangenome.name}")
                future = executor.submit(write_annotations_to_pangenome, pangenome,
                                         metadata, source, force, disable_bar)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def annot_pangenomes_with_hmm(pangenomes: Pangenomes, hmm: Path = None, meta: Path = None,
                              threads: int = 1, task: int = 1, lock: Lock = None,
                              disable_bar: bool = False) -> Dict[str, pd.DataFrame]:
    """ Main function to add annotation to pangenome from tsv file

    :param pangenomes: Pangenomes obejct containing all the pangenome to annotate
    :param hmm: Path to hmm file or directory
    :param meta: Metadata file to annotate with HMM
    :param threads: Number of available threads
    :param lock: Lock for multiprocessing execution
    :param task: Number of task to split processes
    :param disable_bar: Disable bar

    :return: Dictionnary with for each pangenome a dataframe containing families metadata given by HMM
    """
    logging.info("Begin HMM searching")
    if meta is not None:
        metadata = read_metadata(meta)
    else:
        metadata = None
    # Get list of HMM with Plan7 data model
    hmms = read_hmm(hmm_dir=hmm.iterdir(), disable_bar=disable_bar)
    with ThreadPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []

            for pangenome in pangenomes:
                logging.debug(f"Align gene families to HMM for {pangenome.name} with {threads // task} threads...")
                future = executor.submit(annot_with_hmm, pangenome, hmms, metadata, threads // task, disable_bar)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            results = {}
            for future in futures:
                res = future.result()
                results[res[1]] = res[0]
    return results


def annot_pangenomes(pangenomes: Pangenomes, hmm: Path, metadata_list: Path, meta: Path = None, source: str = None,
                     threads: int = 1, task: int = 1, lock: Lock = None,
                     force: bool = False, disable_bar: bool = False):
    """Gene families annotation with HMM or TSV files for multiple pangenomes in multiprocessing


    :param pangenomes: Pangenomes obejct containing all the pangenome to annotate
    :param hmm: Path to HMM directory
    :param metadata_list: Path to list file with for each pangenome an associated metadata file for gene families
    :param meta: Path to metadata related to HMM
    :param source: Name of the annotation source
    :param threads: Number of available threads
    :param lock: Lock for multiprocessing execution
    :param task: Number of task to split processes
    :param force: Boolean to allow force write in pangenomes
    :param disable_bar: Allow to disable progress bar
    """
    assert not all(x is not None for x in [metadata_list, hmm]), "TSV and HMM given to assign metadata. " \
                                                                 "Should be only one !"

    if metadata_list is not None:
        pangenomes2metadata = read_families_metadata_mp(pangenomes, metadata_list, threads, disable_bar)
    elif hmm is not None:
        pangenomes2metadata = annot_pangenomes_with_hmm(pangenomes, hmm, meta=meta, threads=threads,
                                                        task=task, lock=lock, disable_bar=disable_bar)
    else:
        raise Exception("You did not provide tsv or hmm for annotation")
    write_annotations_to_pangenomes(pangenomes, pangenomes2metadata, source, threads, lock, force, disable_bar)


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    manager = Manager()
    lock = manager.Lock()
    pangenomes = load_multiple_pangenomes(pangenome_list=args.pangenomes, need_info={"need_annotations": True,
                                                                                     "need_families": True},
                                          lock=lock, max_workers=args.threads,
                                          check_function=check_pangenome_annotation,
                                          source=args.source, force=args.force, disable_bar=args.disable_prog_bar)

    annot_pangenomes(pangenomes=pangenomes, hmm=args.hmm, metadata_list=args.tsv, source=args.source, meta=args.meta,
                     task=args.task, threads=args.threads, lock=lock, disable_bar=args.disable_prog_bar)


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
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--tsv', type=Path, nargs='?', default=None,
                                help='List of Gene families annotation in TSV format. One TSV per pangenome.'
                                     'See our github for more detail about format')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?', default=None,
                                help="File with all HMM or a directory with one HMM by file")
    hmm_param = parser.add_argument_group(title="HMM arguments",
                                          description="All of the following arguments are required,"
                                                      " if you're using HMM mode :")
    hmm_param.add_argument("--meta", required=False, type=Path, default=None,
                           help="Metadata link to HMM. If no one is given information in HMM will be used")
    hmm_param.add_argument("--mode", required=False, type=str, default='fast', choices=['fast', 'profile'],
                           help="Choose the mode use to align HMM database and gene families. "
                                "Fast will align the reference sequence of gene family against HMM."
                                "Profile will create an HMM profile for each gene family and "
                                "this profile will be aligned")
    hmm_param.add_argument("--msa", required=False, type=Path, default=None,
                           help=argparse.SUPPRESS)
    # help="To create a HMM profile for families, you can give a msa of each gene in families."
    #      "This msa could be gotten from ppanggolin (See ppanggolin msa). "
    #      "If no msa provide Panorama will launch one.")
    hmm_param.add_argument("--msa-format", required=False, type=str, default="afa",
                           choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
                                    "clustal", "clustallike", "phylip", "phylips"],
                           help=argparse.SUPPRESS)
    # help="Format of the input MSA.")
    optional = parser.add_argument_group(title="Optional arguments")
    # optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
    #                       help="directory for storing temporary files")
    optional.add_argument("--task", required=False, nargs='?', type=int, default=1,
                          help="Number of simultaneous task.")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads. If task is more than one,"
                               " threads will be divides by task.")


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
