#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from typing import Dict
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import Manager, Lock, Process, Pool

# installed libraries
from tqdm import tqdm
import pandas as pd
from ppanggolin.meta.meta import check_metadata_format, assign_metadata, write_pangenome

# local libraries
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.utils import init_lock
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


def write_annotations_to_pangenome(pangenome: Pangenome, metadata: pd.DataFrame, source: str,
                                   force: bool = False, disable_bar: bool = False):
    assign_metadata(metadata, pangenome, source, "families", omit=True, disable_bar=disable_bar)
    write_pangenome(pangenome, pangenome.file, force=force, disable_bar=disable_bar)


def write_annotations_to_pangenomes(pangenomes: Pangenomes, pangenomes2metadata: Dict[str, pd.DataFrame],
                                    source: str, threads: int = 1, lock: Lock = None,
                                    force: bool = False, disable_bar: bool = False):
    lock = Lock()
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            _ = []
            for pangenome_name, metadata in pangenomes2metadata.items():
                pangenome = pangenomes.get_pangenome(pangenome_name)
                logging.debug(f"Write annotation for pangenome {pangenome.name}")
                future = executor.submit(write_annotations_to_pangenome, pangenome, metadata, source, force,
                                         disable_bar)
                future.add_done_callback(lambda p: progress.update())
                _.append(future)


def annot_pangenomes_with_hmm(pangenomes: Pangenomes, hmm: Path = None, meta: Path = None,
                              threads: int = 1, task: int = 1, lock: Lock = None,
                              disable_bar: bool = False) -> Dict[str, pd.DataFrame]:
    """ Main function to add annotation to pangenome from tsv file

    :param pangenome: Pangenome object to ppanggolin
    :param hmm: Path to hmm file or directory
    :param tsv: Path to tsv with gene families annotation
    :param meta: Metadata file to annotate with HMM
    :param mode: Which way will be used to annotate families
    :param msa: Path to MSA files
    :param source: Source of annotation
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param disable_bar: Disable bar
    """
    logging.info("Begin HMM searching")
    if meta is not None:
        metadata = read_metadata(meta)
    else:
        metadata = None
    # Get list of HMM with Plan7 data model
    hmms = read_hmm(hmm_dir=hmm.iterdir(), disable_bar=disable_bar)

    with Pool(processes=task) as pool:
        jobs = []
        for pangenome in pangenomes:
            logging.debug(f"Align gene families to HMM for {pangenome.name} with {threads // task} threads...")
            jobs.append(pool.apply_async(func=annot_with_hmm,
                                         args=(pangenome, hmms, metadata, threads // task, disable_bar)))
        results = {}
        for job in tqdm(jobs, unit='pangenome', disable=disable_bar):
            res = job.get()
            results[res[1]] = res[0]
    return results


def annot_pangenomes(pangenomes: Pangenomes, hmm: Path, tsv: Path, meta: Path = None, source: str = None,
                     lock: Lock = None, threads: int = 1, task: int = 1,
                     force: bool = False, disable_bar: bool = False):
    assert not all(x is not None for x in [tsv, hmm]), "TSV and HMM given to assign metadata. Should be only one !"

    if tsv is not None:
        pangenomes2metadata = check_metadata_format(tsv, "families")
    elif hmm is not None:
        pangenomes2metadata = annot_pangenomes_with_hmm(pangenomes, hmm, meta=meta, threads=threads,
                                                        task=task, disable_bar=disable_bar)
    else:
        raise Exception("You did not provide tsv or hmm for annotation")
    write_annotations_to_pangenomes(pangenomes, pangenomes2metadata, source, threads, lock, force, disable_bar)


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    pangenomes = Pangenomes()
    manager = Manager()
    lock = manager.Lock()
    pan_list = load_multiple_pangenomes(pangenome_list=args.pangenomes, need_info={"need_families": True},
                                        lock=lock, max_workers=args.threads, check_function=check_pangenome_annotation,
                                        source=args.source, force=args.force, disable_bar=args.disable_prog_bar)
    pangenomes.add_list_pangenomes(pan_list)
    annot_pangenomes(pangenomes=pangenomes, hmm=args.hmm, tsv=args.tsv, source=args.source, meta=args.meta,
                     lock=lock, task=args.task, threads=args.threads, disable_bar=args.disable_prog_bar)


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
                                help='Gene families annotation in TSV file. See our github for more detail about format')
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
