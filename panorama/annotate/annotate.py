#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations

import argparse
import logging
import shutil
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock, Manager
from pathlib import Path
from typing import Any, Dict, Tuple

# installed libraries
import pandas as pd
from numpy import nan
from ppanggolin.meta.meta import assign_metadata, check_metadata_format
from tqdm import tqdm

# local libraries
from panorama.annotate.hmm_search import annot_with_hmm, read_hmms
from panorama.format.read_binaries import load_pangenomes
from panorama.format.write_binaries import erase_pangenome, write_pangenome
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.utils import init_lock, is_empty, mkdir


def check_annotate_args(
    args: argparse.Namespace, silence_warning: bool = False
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Checks the provided arguments to ensure that they are valid.

    Args:
        args (argparse.Namespace): The parsed arguments.
        silence_warning (bool, optional): Flag to silence warning messages. Defaults to False.
            This option is used for pansystems workflow to not have unwanted warnings.

    Returns:
        Tuple[Dict[str, Any], Dict[str, Any]]:
        Two dictionaries containing necessary information and HMM keyword arguments.

    Raises:
        argparse.ArgumentError: If any required arguments are missing or invalid.
    """
    if args.table is None and args.hmm is None:
        raise argparse.ArgumentError(
            argument=None,
            message="Please provide either a table with annotation for gene families or a hmms to annotate them.",
        )
    need_info = {"need_families": True}
    hmm_kwgs = {}
    if args.table is not None:
        if args.mode is not None:
            logging.getLogger("PANORAMA").error("You cannot specify both --table and --mode in the same command")
            raise argparse.ArgumentError(argument=None, message="--table is incompatible option with '--mode'.")
        if is_empty(args.table):
            raise IOError(f"File: {args.table} is empty.")
        if args.k_best_hit is not None:
            logging.getLogger("PANORAMA").error("You cannot specify both --table and --k_best_hit in the same command")
            raise argparse.ArgumentError(
                argument=None,
                message="--table is Incompatible option with '--k_best_hit'.",
            )
        if args.only_best_hit:
            logging.getLogger("PANORAMA").error(
                "You cannot specify both --table and --only_best_hit in the same command"
            )
            raise argparse.ArgumentError(
                argument=None,
                message="--table is Incompatible option with '--only_best_hit'.",
            )
        if args.output and not silence_warning:
            logging.getLogger("PANORAMA").warning("--output option is incompatible with --table.")

    else:  # args.hmm is not None
        # todo: See to make this the default value in argparse
        args.mode = "fast" if args.mode is None else args.mode

        hmm_kwgs["mode"] = args.mode

        if not args.hmm.exists():
            raise IOError(f"File: {args.hmm} not exists.")
        if is_empty(args.hmm):
            raise IOError(f"File: {args.hmm} is empty.")
        if args.mode == "fast":
            need_info["need_families_sequences"] = True

        else:  # args.mode == "profile" or "sensitive"
            need_info["need_annotations"] = True
            need_info["need_gene_sequences"] = True

        if args.mode == "fast" and args.keep_tmp:
            logging.getLogger("PANORAMA").warning("--keep_tmp is not working with --mode fast")
        hmm_kwgs["tmp"] = Path(tempfile.mkdtemp(prefix="panorama_tmp", dir=args.tmp))

        if args.msa is not None and args.mode != "profile":
            raise argparse.ArgumentError(argument=None, message="--msa is working only with --profile")

        hmm_kwgs["msa"] = args.msa
        hmm_kwgs["msa_format"] = args.msa_format
        if args.Z is not None and not isinstance(args.Z, int):
            raise TypeError("Z arguments must be an integer")
        hmm_kwgs["Z"] = args.Z
        if args.k_best_hit is not None:
            if args.k_best_hit < 1:
                raise argparse.ArgumentTypeError("k_best_hit must be greater than 1.")
            if args.only_best_hit:
                if args.k_best_hit == 1:
                    logging.getLogger("PANORAMA").warning(
                        "'--only_best_hit' is an alias for '--k_best_hit 1'. You can use only one of them."
                    )
                else:
                    logging.getLogger("PANORAMA").error(
                        f"You set the maximum of hit at {args.k_best_hit} but you also "
                        f"use '--only_best_hit' which is an alias for '--k_best_hit 1'."
                        f" Please use only one of those."
                    )
                    raise argparse.ArgumentError(
                        argument=None,
                        message="--k_best_hit is incompatible with --only_best_hit",
                    )
        else:
            if args.only_best_hit:
                args.k_best_hit = 1

        if args.save_hits is not None:
            if args.output is None:
                raise argparse.ArgumentError(argument=None, message="--output is required to save hits results.")
            else:
                hmm_kwgs["output"] = mkdir(args.output, force=args.force, erase=False)
            if "tblout" in args.save_hits:
                hmm_kwgs["tblout"] = True
            if "domtblout" in args.save_hits:
                hmm_kwgs["domtblout"] = True
            if "pfamtblout" in args.save_hits:
                hmm_kwgs["pfamtblout"] = True
        else:
            if args.output and not silence_warning:
                logging.getLogger("PANORAMA").warning("--output option is compatible only with --save_hits.")
    return need_info, hmm_kwgs


def check_pangenome_annotation(pangenome: Pangenome, source: str, force: bool = False):
    """
    Check pangenome information before adding annotation.

    Args:
        pangenome (Pangenome): Pangenome object that will be checked.
        source (str): Source of annotation to check if already in pangenome.
        force (bool, optional): Flag to allow overwriting/erasing annotation. Defaults to False.

    Raises:
        KeyError: If a source with the same name already exists, and force is False.
    """
    if pangenome.status["metadata"]["families"] == "inFile" and source in pangenome.status["metasources"]["families"]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source)
        else:
            raise KeyError(
                f"A metadata corresponding to the source : '{source}' already exist in pangenome."
                f"Add the option --force to erase"
            )


def read_families_metadata(pangenome: Pangenome, metadata: Path) -> Tuple[pd.DataFrame, str]:
    """
    Read gene families metadata for one pangenome.

    Args:
        pangenome (Pangenome): Pangenome object for which metadata will be associated.
        metadata (Path): Path to metadata file containing metadata to add to pangenome.

    Returns:
        Tuple[pd.DataFrame, str]: The metadata dataframe and the name of the pangenome.
    """
    metadata_df = check_metadata_format(metadata, "families")
    return metadata_df, pangenome.name


def read_families_metadata_mp(
    pangenomes: Pangenomes,
    table: Path,
    threads: int = 1,
    lock: Lock = None,
    disable_bar: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Read gene families metadata for multiple pangenomes in multiprocessing.

    Args:
        pangenomes (Pangenomes): Pangenomes object containing all the pangenome to annotate.
        table (Path): Path to the metadata file for gene families.
        threads (int, optional): Number of available threads. Defaults to 1.
        lock (Lock, optional): Lock for multiprocessing execution. Defaults to None.
        disable_bar (bool, optional): Flag to disable progress bar. Defaults to False.

    Returns:
        Dict[str, pd.DataFrame]: Dictionary with the metadata linked to pangenome by its name.
    """
    t0 = time.time()
    path_to_metadata = pd.read_csv(table, delimiter="\t", names=["Pangenome", "path"])
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit="pangenome", disable=disable_bar) as progress:
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
    logging.getLogger("PANORAMA").info(f"Pangenomes annotation ridden from TSV file in {time.time() - t0:2f} seconds")
    return results


def get_k_best_hit(group, k_best_hit: int) -> pd.DataFrame:
    """
    Get the K best hits for a given group in a dataframe.

    Args:
        group: Dataframe group.
        k_best_hit (int): Number of best hits to keep.

    Returns:
        pd.DataFrame: K best hits per group.
    """
    return group.nlargest(k_best_hit, columns=["score", "bias", "e_value"])


def remove_redundant_annotation(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Remove redundant annotation based on score, e-value, and bias.

    Args:
        metadata (pd.DataFrame): Metadata dataframe containing annotations.

    Returns:
        pd.DataFrame: Dataframe with redundant annotations removed.
    """
    logging.getLogger("PANORAMA").debug("Remove duplicate hits and keep the best score")
    metadata_df = metadata.sort_values(by=["score", "e_value", "bias"], ascending=[False, True, False])
    group = metadata_df.groupby(["families", "protein_name"])
    metadata_df = group.first().assign(
        secondary_names=group.agg({"secondary_names": lambda x: ",".join(set(x.dropna())) or nan})
    )
    metadata_df = metadata_df.reset_index()
    return metadata_df


def keep_best_hit(metadata: pd.DataFrame, k_best_hit: int) -> pd.DataFrame:
    """
    Keep the k best hit for a given metadata.

    Args:
        metadata (pd.DataFrame): Metadata dataframe with multiple annotations for gene families.
        k_best_hit (int): Number of best hits to keep.

    Returns:
        pd.DataFrame: Filtered metadata dataframe with only the k best hits.
    """
    logging.getLogger("PANORAMA").debug(f"keep the {k_best_hit} best hits")
    # Sort by the priority criteria
    metadata_sorted = metadata.sort_values(
        [
            "score",  # Highest score first (descending)
            "bias",  # Lowest bias first (ascending)
            "e_value",  # Lowest e_value first (ascending)
            "i_e_value",  # Lowest i_e_value first (ascending)
        ],
        ascending=[False, True, True, True],
    )
    return metadata_sorted.groupby("families").head(k_best_hit).reset_index(drop=True)


def write_annotations_to_pangenome(
    pangenome: Pangenome,
    metadata: pd.DataFrame,
    source: str,
    k_best_hit: int = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Write gene families annotation for one pangenome.

    Args:
        pangenome (Pangenome): Pangenome linked to metadata.
        metadata (pd.DataFrame): Metadata dataframe.
        source (str): Metadata source.
        k_best_hit (int, optional): Number of best hits to keep. Defaults to None.
        force (bool, optional): Boolean to allow force writing in pangenomes. Defaults to False.
        disable_bar (bool, optional): Allow disabling the progress bar. Defaults to False.
    """

    logging.getLogger("PANORAMA").debug("Remove duplicate hits and keep the best score")
    meta_df = remove_redundant_annotation(metadata)
    if k_best_hit is not None:
        meta_df = keep_best_hit(meta_df, k_best_hit)
    meta_df = meta_df.sort_values(by=["score", "e_value", "bias"], ascending=[False, True, False])
    assign_metadata(meta_df, pangenome, source, "families", omit=True, disable_bar=disable_bar)
    write_pangenome(pangenome, pangenome.file, force=force, disable_bar=disable_bar)


def write_annotations_to_pangenomes(
    pangenomes: Pangenomes,
    pangenomes2metadata: Dict[str, pd.DataFrame],
    source: str,
    k_best_hit: int = None,
    threads: int = 1,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Write gene families annotation for pangenomes in multiple processing.

    Args:
        pangenomes (Pangenomes): Pangenomes object containing all the pangenome to annotate.
        pangenomes2metadata (Dict[str, pd.DataFrame]): Dictionary with for each pangenome
        the metadata dataframe associated.
        source (str): Metadata source.
        k_best_hit (int, optional): Number of best hits to keep. Defaults to None.
        threads (int, optional): Number of available threads. Defaults to 1.
        lock (Lock, optional): Lock for multiprocessing execution. Defaults to None.
        force (bool, optional): Boolean to allow force to write in pangenomes. Defaults to False.
        disable_bar (bool, optional): Allow disabling the progress bar. Defaults to False.
    """
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit="pangenome", disable=disable_bar) as progress:
            futures = []
            for pangenome_name, metadata in pangenomes2metadata.items():
                pangenome = pangenomes.get(pangenome_name)
                logging.getLogger("PANORAMA").debug(f"Write annotation for pangenome {pangenome.name}")
                future = executor.submit(
                    write_annotations_to_pangenome,
                    pangenome,
                    metadata,
                    source,
                    k_best_hit,
                    force,
                    disable_bar,
                )
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def annot_pangenomes_with_hmm(
    pangenomes: Pangenomes,
    hmm: Path = None,
    source: str = "",
    mode: str = "fast",
    threads: int = 1,
    disable_bar: bool = False,
    **hmm_kwgs,
) -> Dict[str, pd.DataFrame]:
    """
    Main function to add annotation to pangenome from tsv file.

    Args:
        pangenomes (Pangenomes): Pangenomes object containing all the pangenome to annotate.
        hmm (Path, optional): Path to hmm list file. Defaults to None.
        source (str, optional): Name of the annotation source. Defaults to "".
        mode (str, optional): Which mode to use to annotate gene families with HMM. Defaults to "fast".
        threads (int, optional): Number of available threads. Defaults to 1.
        disable_bar (bool, optional): Flag to disable progress bar. Defaults to False.
        **hmm_kwgs: Arbitrary keyword arguments for HMM annotation.

    Returns:
        Dict[str, pd.DataFrame]: Dictionary with for each pangenome a dataframe
        containing gene families metadata given by HMM.
    """
    t0 = time.time()
    logging.getLogger("PANORAMA").info("Begin HMM searching")
    # Get the list of HMM with Plan7's data model
    pangenome2annot = {}
    hmms, hmm_df = read_hmms(hmm, disable_bar=disable_bar)
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit="pangenome", disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Align gene families to HMM for {pangenome.name}")
        pangenome2annot[pangenome.name] = annot_with_hmm(
            pangenome,
            hmms,
            hmm_df,
            source,
            mode,
            threads=threads,
            disable_bar=disable_bar,
            **hmm_kwgs,
        )
    logging.getLogger("PANORAMA").info(f"Pangenomes annotation with HMM done in {time.time() - t0:2f} seconds")
    return pangenome2annot


def annot_pangenomes(
    pangenomes: Pangenomes,
    source: str = None,
    table: Path = None,
    hmm: Path = None,
    threads: int = 1,
    k_best_hit: int = None,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
    **hmm_kwgs: Any,
):
    """
    Gene families annotation with HMM or TSV files for multiple pangenomes in multiprocessing.

    Args:
        pangenomes (Pangenomes):
            Pangenomes object containing all the pangenome to annotate.
        source (str, optional):
            Name of the annotation source. Defaults to None.
        table (Path, optional):
            Path to the metadata file for gene families annotation. Defaults to None.
        hmm (Path, optional):
            Path to hmm list file. Defaults to None.
        threads (int, optional):
            Number of available threads. Defaults to 1.
        k_best_hit (int, optional):
            Number of best hits to keep. Defaults to None.
        lock (Lock, optional):
            Lock for multiprocessing. Defaults to None.
        force (bool, optional):
            Flag to allow force overwrite in pangenomes. Defaults to False.
        disable_bar (bool, optional):
            Flag to disable progress bar. Defaults to False.
        **hmm_kwgs (Any):
            Arbitrary keyword arguments for hmm alignment.

    Raises:
        AssertionError: If neither HMM nor TSV are provided.
    """
    assert table is not None or hmm is not None, "Must provide either table or hmm"
    if table is not None:
        pangenomes2metadata = read_families_metadata_mp(pangenomes, table, threads, disable_bar)
    else:  # hmm is not None:
        pangenomes2metadata = annot_pangenomes_with_hmm(
            pangenomes,
            hmm,
            source,
            threads=threads,
            disable_bar=disable_bar,
            **hmm_kwgs,
        )
    try:
        write_annotations_to_pangenomes(
            pangenomes,
            pangenomes2metadata,
            source,
            k_best_hit,
            threads,
            lock,
            force,
            disable_bar,
        )
    except Exception as e:
        for pangenome in pangenomes:
            erase_pangenome(pangenome, metadata=True, source=source)
            logging.getLogger("PANORAMA").debug(f"erase annotation from {source} in {pangenome.name}")
        raise Exception(f"Annotation failed from : {e}") from e


def launch(args: argparse.Namespace) -> None:
    """
    Launch functions to annotate pangenomes

    Args:
        args  (argparse.Namespace): argument given in CLI
    """
    manager = Manager()
    lock = manager.Lock()
    need_info, hmm_kwgs = check_annotate_args(args)
    try:
        pangenomes = load_pangenomes(
            pangenome_list=args.pangenomes,
            need_info=need_info,
            check_function=check_pangenome_annotation,
            max_workers=args.threads,
            lock=lock,
            disable_bar=args.disable_prog_bar,
            source=args.source,
            force=args.force,
        )
        annot_pangenomes(
            pangenomes=pangenomes,
            source=args.source,
            table=args.table,
            hmm=args.hmm,
            threads=args.threads,
            k_best_hit=args.k_best_hit,
            lock=lock,
            force=args.force,
            disable_bar=args.disable_prog_bar,
            **hmm_kwgs,
        )
    except Exception as e:
        raise Exception(f"Annotation failed from : {e}") from e
    finally:
        if not args.keep_tmp:
            shutil.rmtree(hmm_kwgs["tmp"])
        else:
            logging.getLogger("PANORAMA").info(f"Temporary file has been saved here: {hmm_kwgs['tmp'].as_posix()}")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA Command line

    Args:
        sub_parser: sub_parser for annot command

    Returns:
        argparse.ArgumentParser: parser arguments for annot command
    """
    parser = sub_parser.add_parser("annotation")
    parser_annot(parser)
    return parser


def parser_annot_hmm(parser):
    """
    Add argument to parser for HMM annotation

    Args:
        parser: parser for annot argument
    """
    hmm_param = parser.add_argument_group(
        title="HMM arguments",
        description="All of the following arguments are required, if you're using HMM mode :",
    )
    hmm_param.add_argument(
        "--mode",
        required=False,
        type=str,
        default=None,
        choices=["fast", "profile", "sensitive"],
        help="Choose the mode use to align HMM database and gene families. "
        "Fast will align the reference sequence of gene family against HMM."
        "Profile will create an HMM profile for each gene family and "
        "this profile will be aligned."
        "Sensitive will align HMM to all genes in families.",
    )
    hmm_param.add_argument(
        "--k_best_hit",
        required=False,
        type=int,
        default=None,
        help="Keep the k best annotation hit per gene family.If not specified, all hit will be kept.",
    )
    hmm_param.add_argument(
        "-b",
        "--only_best_hit",
        required=False,
        action="store_true",
        help="alias to keep only the best hit for each gene family.",
    )
    hmm_param.add_argument(
        "--msa",
        required=False,
        type=Path,
        default=None,
        help="To create a HMM profile for families, you can give a msa of each gene in families."
        "This msa could be gotten from ppanggolin (See ppanggolin msa). "
        "If no msa provide Panorama will launch one.",
    )
    hmm_param.add_argument(
        "--msa-format",
        required=False,
        type=str,
        default="afa",
        choices=[
            "stockholm",
            "pfam",
            "a2m",
            "psiblast",
            "selex",
            "afa",
            "clustal",
            "clustallike",
            "phylip",
            "phylips",
        ],
        help=argparse.SUPPRESS,
    )
    hmm_param.add_argument(
        "--save_hits",
        required=False,
        type=str,
        default=None,
        nargs="*",
        choices=["tblout", "domtblout", "pfamtblout"],
        help="Save HMM alignment results in tabular format. Option are the same than in HMMSearch.",
    )
    hmm_param.add_argument(
        "-Z",
        "--Z",
        required=False,
        type=int,
        nargs="?",
        default=None,
        help="From HMMER: Assert that the total number of targets in your searches is <x>, "
        "for the purposes of per-sequence E-value calculations, "
        "rather than the actual number of targets seen.",
    )
    return hmm_param


def parser_annot(parser):
    """
    Add argument to parser for annot command

    Args:
        parser: parser for annot argument
    """
    required = parser.add_argument_group(
        title="Required arguments",
        description="All of the following arguments are required :",
    )
    required.add_argument(
        "-p",
        "--pangenomes",
        required=True,
        type=Path,
        nargs="?",
        help="A list of pangenome.h5 files in .tsv file",
    )
    required.add_argument(
        "-s",
        "--source",
        required=True,
        type=str,
        nargs="?",
        help="Name of the annotation source.",
    )
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument(
        "--table",
        type=Path,
        default=None,
        help="A list of tab-separated file, containing annotation of gene families."
        "Expected format is pangenome name in first column "
        "and path to the TSV with annotation in second column.",
    )
    exclusive_mode.add_argument(
        "--hmm",
        type=Path,
        nargs="?",
        default=None,
        help="A tab-separated file with HMM information and path."
        "Note: Use panorama utils --hmm to create the HMM list file",
    )
    hmm_param = parser_annot_hmm(parser)
    hmm_param.add_argument(
        "-o",
        "--output",
        required=False,
        type=Path,
        nargs="?",
        help="Output directory to write HMM results",
    )
    hmm_param.add_argument(
        "--keep_tmp",
        required=False,
        action="store_true",
        help="Keep the temporary files. Useful for debugging in sensitive or profile mode.",
    )
    hmm_param.add_argument(
        "--tmp",
        required=False,
        nargs="?",
        type=Path,
        default=Path(tempfile.gettempdir()),
        help=f"Path to temporary directory, defaults path is {Path(tempfile.gettempdir()) / 'panorama'}",
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "--threads",
        required=False,
        nargs="?",
        type=int,
        default=1,
        help="Number of available threads.",
    )
