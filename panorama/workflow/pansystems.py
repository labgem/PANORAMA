#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to detect biological systems in pangenomes all in one workflow command
"""

# default libraries
from __future__ import annotations
import logging
import argparse
from typing import List
from pathlib import Path
from multiprocessing import Manager
import tempfile
from typing import Any, Dict, Tuple
import time

# installed libraries
from tqdm import tqdm
from pyhmmer.plan7 import HMM
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.format.write_binaries import erase_pangenome
from panorama.annotate.annotate import (check_annotate_args, check_pangenome_annotation, read_families_metadata,
                                        write_annotations_to_pangenome, parser_annot_hmm)
from panorama.annotate.hmm_search import read_hmms, annot_with_hmm
from panorama.systems.models import Models
from panorama.systems.detection import check_detection_args, search_systems, write_systems_to_pangenome
from panorama.systems.write_systems import write_flat_systems_to_pangenome
from panorama.utility.utility import check_models, mkdir


def check_pansystems_parameters(args: argparse.Namespace) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Checks and validates the parameters for the pansystems function.

    Args:
        args (argparse.Namespace): The arguments passed to the function.

    Returns:
        Tuple[Dict[str, Any], Dict[str, Any]]: A tuple containing the necessary information and HMM keyword arguments.

    Raises:
        argparse.ArgumentError: If no type of systems writing is chosen.
    """
    need_info, hmm_kwgs = check_annotate_args(args)
    if 'output' not in hmm_kwgs:
        hmm_kwgs["output"] = mkdir(args.output, force=args.force, erase=False)
    args.annotation_sources = None
    need_info.update(check_detection_args(args))
    need_info["need_metadata"] = False
    need_info["need_families_info"] = True
    if not any(arg for arg in [args.projection, args.partition, args.association, args.proksee]):
        raise argparse.ArgumentError(argument=None, message="You should at least choose one type of systems writing "
                                                            "between: projection, partition, association or proksee.")

    if "all" in args.association:
        args.association = ["RGPs", "spots", "modules"]

    for asso in args.association:
        if asso == "modules":
            need_info["need_modules"] = True
        if asso in ["RGPs", "spots"]:
            need_info["need_rgp"] = True
        if asso == "spots":
            need_info["need_spots"] = True
    return need_info, hmm_kwgs


def check_pangenome_pansystems(pangenome: Pangenome, source: str, force: bool = False) -> None:
    """
    Checks the annotation of a pangenome and its systems.

    Args:
        pangenome (Pangenome): The pangenome to check.
        source (str): The source of the annotation.
        force (bool, optional): Whether to force the erased of already computed systems. Defaults to False.

    Raises:
        ValueError: If systems are already detected based on the source and force is False.
    """
    check_pangenome_annotation(pangenome, source, force=force)
    if pangenome.status["systems"] == "inFile" and source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=source)
        else:
            raise ValueError(f"Systems are already detected based on the source : {source}. "
                             f"Use the --force option to erase the already computed systems.")


def pansystems_pangenome(pangenome: Pangenome, source: str, models: Models, table: Path = None,
                         hmms: Dict[str, List[HMM]] = None, k_best_hit: int = None, jaccard_threshold: float = 0.8,
                         sensitivity: int = 1, projection: bool = False, association: List[str] = None,
                         partition: bool = False, proksee: str = None, threads: int = 1, force: bool = False,
                         disable_bar: bool = False, **hmm_kwgs: Any) -> None:
    """
    Detects systems in a single pangenome.

    Args:
        pangenome (Pangenome): The pangenome to analyze.
        source (str): The source of the annotation.
        models (Models): The models to use for systems detection.
        table (Path, optional): The path to a table with annotation information. Defaults to None.
        hmms (Dict[str, List[HMM]], optional): The HMMs to use for annotation. Defaults to None.
        k_best_hit (int, optional): The number of best annotation hits to keep per gene family. Defaults to None.
        jaccard_threshold (float, optional): The minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        projection (bool, optional): Whether to project the systems on organisms. Defaults to False.
        association (List[str], optional): The type of association to write between systems and other pangenome elements. Defaults to None.
        partition (bool, optional): Whether to write a heatmap file with for each organism, partition of the systems. Defaults to False.
        proksee (str, optional): Whether to write a proksee file with systems. Defaults to None.
        threads (int, optional): The number of available threads. Defaults to 1.
        force (bool, optional): Whether to force the erased of already computed systems. Defaults to False.
        disable_bar (bool, optional): Whether to disable the progress bar. Defaults to False.
        **hmm_kwgs (Any): Additional keyword arguments for HMM annotation.
    """
    if table is not None:
        metadata_df, _ = read_families_metadata(pangenome, table)
    else:
        t0 = time.time()
        metadata_df = annot_with_hmm(pangenome, hmms, source=source, threads=threads,
                                     disable_bar=disable_bar, **hmm_kwgs)
        logging.getLogger("PANORAMA").info(f"Pangenomes annotation with HMM done in {time.time() - t0:2f} seconds")
    write_annotations_to_pangenome(pangenome, metadata_df, source, k_best_hit, force, disable_bar)
    search_systems(models, pangenome, source, [source], jaccard_threshold, sensitivity,
                   threads=threads, disable_bar=disable_bar)
    write_systems_to_pangenome(pangenome, source, disable_bar=disable_bar)
    write_flat_systems_to_pangenome(pangenome, hmm_kwgs["output"], projection, association, partition, proksee,
                                    threads=threads, force=force, disable_bar=disable_bar)


def pansystems(pangenomes: Pangenomes, source: str, models: Models, hmm: Dict[str, List[HMM]] = None,
               table: pd.DataFrame = None, k_best_hit: int = None, jaccard_threshold: float = 0.8, sensitivity: int = 1,
               projection: bool = False, association: List[str] = None, partition: bool = False, proksee: str = None,
               threads: int = 1, force: bool = False, disable_bar: bool = False, **hmm_kwgs: Any) -> None:
    """
    Detects systems in multiple pangenomes.

    Args:
        pangenomes (Pangenomes): The pangenomes to analyze.
        source (str): The source of the annotation.
        models (Models): The models to detect systems.
        table (pd.Dataframe, optional): Dataframe containing for each pangenome a path to a table with annotation information. Defaults to None.
        hmm (Dict[str, List[HMM]], optional): A dictionary to identify which cutoff use to align HMM . Defaults to None.
        k_best_hit (int, optional): The number of best annotation hits to keep per gene family. Defaults to None.
        jaccard_threshold (float, optional): The minimum Jaccard similarity used to filter edges between gene families. Defaults to 0.8.
        projection (bool, optional): Whether to project the systems on organisms. Defaults to False.
        association (List[str], optional): The type of association to write between systems and other pangenome elements. Defaults to None.
        partition (bool, optional): Whether to write a heatmap file with for each organism, partition of the systems. Defaults to False.
        proksee (str, optional): Whether to write a proksee file with systems. Defaults to None.
        threads (int, optional): The number of available threads. Defaults to 1.
        force (bool, optional): Whether to force the erased of already computed systems. Defaults to False.
        disable_bar (bool, optional): Whether to disable the progress bar. Defaults to False.
        **hmm_kwgs (Any): Additional keyword arguments for HMM annotation.

    Raises:
        AssertionError: If neither table nor hmm is provided.
    """
    assert table is not None or hmm is not None, 'Must provide either table or hmm'
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        if table is not None:
            metadata_file = table.loc[table["Pangenome"] == pangenome.name]["path"].squeeze()
        else:
            metadata_file = None
        pansystems_pangenome(pangenome, source, models, metadata_file, hmm, k_best_hit, jaccard_threshold, sensitivity,
                             projection, association, partition, proksee, threads, force, disable_bar, **hmm_kwgs)


def check_input_files(models_path: Path, table: Path = None, hmm: Path = None,
                      disable_bar: bool = False) -> Tuple[pd.DataFrame, Dict[str, List[HMM]], pd.DataFrame, Models]:
    """
    Check the metadta table, the hmm and the models, in order to stop program before to read pangenome.

    Args:
        models_path (Path): The path to the models list file.
        table (Path, optional): The path to a table with annotation information. Defaults to None.
        hmm (Path, optional): The path to a tab-separated file with HMM information and path. Defaults to None.
        disable_bar (bool, optional): Whether to disable the progress bar. Defaults to False.

    Returns:
        pd.Dataframe, optional: Dataframe containing for each pangenome a path to a table with annotation information.
        Dict[str, List[HMM]]: A dictionary to identify which cutoff use to align HMM.
        pd.Dataframe: Dataframe with hmm metadata information
        Models: The models to detect systems.

    Raises:
        AssertionError: If neither table nor hmm is provided.
    """
    assert table is not None or hmm is not None, 'Must provide either table or hmm'
    if table is not None:
        path_to_metadata = pd.read_csv(table, delimiter="\t", names=["Pangenome", "path"])
        hmms, hmm_info = None, None
    else:
        hmms, hmm_info = read_hmms(hmm, disable_bar=disable_bar)
        path_to_metadata = None
    models = check_models(models_path, disable_bar=disable_bar)
    return path_to_metadata, hmms, hmm_info, models


def launch(args):
    """
    Launch functions to detect systems in pangenomes

    Args:
        args: argument given in CLI
    """
    from panorama.format.read_binaries import load_pangenomes

    need_info, hmm_kwgs = check_pansystems_parameters(args)
    manager = Manager()
    lock = manager.Lock()
    table, hmms, hmm_kwgs["meta"], models = check_input_files(args.models, args.table, args.hmm, args.disable_prog_bar)
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_pansystems, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, source=args.source, force=args.force)
    pansystems(pangenomes, args.source, models, hmms, table, args.k_best_hit, args.jaccard, args.sensitivity,
               args.projection, args.association, args.partition, args.proksee, args.threads,
               args.force, args.disable_prog_bar, **hmm_kwgs)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for systems command

    Returns:
        argparse.ArgumentParser: parser arguments for align command
    """
    parser = sub_parser.add_parser("pansystems")
    parser_pansystems(parser)
    return parser


def parser_pansystems(parser):
    """
    Add argument to parser for systems command

    Args:
        parser: parser for systems argument

    TODO:
        - add an option to write projection
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required:")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source where panorama as to select in pangenomes')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    annotate = parser.add_argument_group(title="Annotation arguments",
                                         description="All of the following arguments are used for annotation step:")
    exclusive_mode = annotate.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--table', type=Path, default=None,
                                help='A list of tab-separated file, containing annotation of gene families.'
                                     'Expected format is pangenome name in first column '
                                     'and path to the TSV with annotation in second column.')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?', default=None,
                                help="A tab-separated file with HMM information and path."
                                 "Note: Use panorama utils --hmm to create the HMM list file")
    parser_annot_hmm(annotate)
    detection = parser.add_argument_group(title="Systems detection arguments",
                                          description="All of the following arguments are used "
                                                      "for systems detection step:")
    detection.add_argument('-m', '--models', required=True, type=Path, nargs='?',
                           help="Path to model list file."
                                "Note: Use panorama utils --models to create the models list file")
    detection.add_argument('--jaccard', required=False, type=float, default=0.8,
                           help="minimum jaccard similarity used to filter edges between gene families. "
                                "Increasing it will improve precision but lower sensitivity a lot.")
    detection.add_argument('--sensitivity', required=False, type=int, default=3, choices=[1, 2, 3],
                           help=argparse.SUPPRESS)
    write = parser.add_argument_group(title="Optional arguments")
    write.add_argument("--projection", required=False, action="store_true",
                       help="Project the systems on organisms. If organisms are specified, "
                            "projection will be done only for them.")
    write.add_argument("--partition", required=False, action="store_true",
                       help="Write a heatmap file with for each organism, partition of the systems. "
                            "If organisms are specified, heatmap will be write only for them.")
    write.add_argument("--association", required=False, type=str, default=[], nargs='+',
                       choices=["all", "modules", "RGPs", "spots"],
                       help="Write association between systems and others pangenomes elements")
    write.add_argument("--proksee", required=False, type=str, default=None, nargs='+',
                       choices=["all", "base", "modules", "RGP", "spots", "annotate"],
                       help="Write a proksee file with systems. "
                            "If you only want the systems with genes, gene families and partition, use base value."
                            "Write RGPs, spots or modules -split by `,'- if you want them.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    optional.add_argument("--keep_tmp", required=False, action='store_true',
                          help="Keep the temporary files. Useful for debugging in sensitive or profile mode.")
    optional.add_argument("--tmp", required=False, nargs='?', type=Path, default=None,
                          help=f"Path to temporary directory, defaults path is {Path(tempfile.gettempdir()) / 'panorama'}")
