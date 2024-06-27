#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to detect biological systems in pangenomes
"""

# default libraries
from __future__ import annotations
import argparse
import logging
from typing import List
from pathlib import Path
from multiprocessing import Manager
import tempfile
from typing import Any, Dict, Tuple

# installed libraries
from tqdm import tqdm
from pyhmmer.plan7 import HMM
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.format.write_binaries import erase_pangenome
from panorama.annotate.annotate import (check_annotate_args, check_pangenome_annotation, read_families_metadata,
                                        write_annotations_to_pangenome)
from panorama.annotate.hmm_search import read_hmms, annot_with_hmm
from panorama.systems.models import Models
from panorama.systems.detection import check_detection_args, search_systems, write_systems_to_pangenome
from panorama.systems.write_systems import write_flat_systems_to_pangenome
from panorama.utility.utility import check_models


def check_pansystems_parameters(args: argparse.Namespace) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    need_info, hmm_kwgs = check_annotate_args(args)
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
    check_pangenome_annotation(pangenome, source, force=force)
    if pangenome.status["systems"] == "inFile" and source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=source)
        else:
            raise ValueError(f"Systems are already detected based on the source : {source}. "
                             f"Use the --force option to erase the already computed systems.")


def pansystems_pangenome(pangenome: Pangenome, source: str, models: Models, table: Path = None,
                         hmms: Dict[str, List[HMM]] = None, k_best_hit: int = None, jaccard_threshold: float = 0.8,
                         projection: bool = False, association: List[str] = None, partition: bool = False,
                         proksee: str = None, threads: int = 1, force: bool = False, disable_bar: bool = False,
                         **hmm_kwgs: Any) -> None:
    if table is not None:
        metadata_df, _ = read_families_metadata(pangenome, table)
    else:
        metadata_df = annot_with_hmm(pangenome, hmms, source=source, threads=threads,
                                     disable_bar=disable_bar, **hmm_kwgs)
    write_annotations_to_pangenome(pangenome, metadata_df, source, k_best_hit, force, disable_bar)
    search_systems(models, pangenome, source, [source], jaccard_threshold, threads=threads, disable_bar=disable_bar)
    write_systems_to_pangenome(pangenome, source, disable_bar=disable_bar)
    write_flat_systems_to_pangenome(pangenome, hmm_kwgs["output"], projection, association, partition, proksee,
                                    threads=threads, force=force, disable_bar=disable_bar)


def pansystems(pangenomes: Pangenomes, source: str, models_path: Path, table: Path = None, hmm: Path = None,
               k_best_hit: int = None, jaccard_threshold: float = 0.8, projection: bool = False,
               association: List[str] = None, partition: bool = False, proksee: str = None, threads: int = 1,
               force: bool = False, disable_bar: bool = False, **hmm_kwgs: Any) -> None:
    assert table is not None or hmm is not None, 'Must provide either table or hmm'
    if table is not None:
        path_to_metadata = pd.read_csv(table, delimiter="\t", names=["Pangenome", "path"])
        hmms = None
    else:
        hmms, hmm_kwgs["meta"] = read_hmms(hmm, disable_bar=disable_bar)
        path_to_metadata = None
    models = check_models(models_path, disable_bar=disable_bar)
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        if path_to_metadata is not None:
            metadata_file = path_to_metadata.loc[path_to_metadata["Pangenome"] == pangenome.name]["path"].squeeze()
        else:
            metadata_file = None
        pansystems_pangenome(pangenome, source, models, metadata_file, hmms, k_best_hit, jaccard_threshold,
                             projection, association, partition, proksee, threads, force, disable_bar, **hmm_kwgs)


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
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_pansystems, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, source=args.source, force=args.force)
    pansystems(pangenomes, args.source, args.models, args.table, args.hmm, args.k_best_hit, args.jaccard,
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
    parser = sub_parser.add_parser("pansystems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    annotations = parser.add_argument_group(title="Annotation arguments",
                                            description="All of the following arguments are used for annotation step:")
    exclusive_mode = annotations.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--table', type=Path, default=None,  # nargs='+',
                                help='A list of tab-separated file, containing annotation of gene families.'
                                     'Expected format is pangenome name in first column '
                                     'and path to the TSV with annotation in second column.')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?', default=None,
                                help="A tab-separated file with HMM information and path."
                                     "Note: Use panorama utils --hmm to create the HMM list file")
    hmm_param = parser.add_argument_group(title="HMM arguments",
                                          description="All of the following arguments are required,"
                                                      " if you're using HMM mode :")
    hmm_param.add_argument("--mode", required=False, type=str, default=None, choices=['fast', 'profile', 'sensitive'],
                           help="Choose the mode use to align HMM database and gene families. "
                                "Fast will align the reference sequence of gene family against HMM."
                                "Profile will create an HMM profile for each gene family and "
                                "this profile will be aligned."
                                "Sensitive will align HMM to all genes in families.")
    hmm_param.add_argument("--k_best_hit", required=False, type=int, default=None,
                           help="Keep the k best annotation hit per gene family."
                                "If not specified, all hit will be kept.")
    hmm_param.add_argument("-b", "--only_best_hit", required=False, action="store_true",
                           help="alias to keep only the best hit for each gene family.")
    hmm_param.add_argument("--msa", required=False, type=Path, default=None,
                           help="To create a HMM profile for families, you can give a msa of each gene in families."
                                "This msa could be gotten from ppanggolin (See ppanggolin msa). "
                                "If no msa provide Panorama will launch one.")
    hmm_param.add_argument("--msa-format", required=False, type=str, default="afa",
                           choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
                                    "clustal", "clustallike", "phylip", "phylips"],
                           help=argparse.SUPPRESS)
    hmm_param.add_argument("--save_hits", required=False, type=str, default=None, nargs='*',
                           choices=['tblout', 'domtblout', 'pfamtblout'],
                           help='Save HMM alignment results in tabular format. Option are the same than in HMMSearch.')
    detection = parser.add_argument_group(title="Systems detection arguments",
                                          description="All of the following arguments are used "
                                                      "for systems detection step:")
    detection.add_argument('-m', '--models', required=True, type=Path, nargs='?',
                           help="Path to model list file."
                                "Note: Use panorama utils --models to create the models list file")
    detection.add_argument('--jaccard', required=False, type=float, default=0.8,
                           help="minimum jaccard similarity used to filter edges between gene families. "
                                "Increasing it will improve precision but lower sensitivity a lot.")
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
                       choices=["all", "base", "modules", "RGP", "spots", "annotations"],
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
