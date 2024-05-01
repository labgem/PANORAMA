#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import argparse
import logging
from typing import Any, Dict, List
from multiprocessing import Manager, Lock
from pathlib import Path

import pandas as pd
# installed libraries
from tqdm import tqdm

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.systems.systems_projection import project_pangenome_systems, write_projection_systems
from panorama.systems.systems_partitions import systems_partition
from panorama.systems.system_association import association_pangenome_systems


def check_write_systems_args(args: argparse.Namespace) -> Dict[str, Any]:
    """ Checks the provided arguments to ensure that they are valid.

    Args:
        args: The parsed arguments.

    Raises:
        argparse.ArgumentTypeError: If number of sources is not the same as models
        argparse.ArgumentTypeError: if annotation are given and their number is not the same as systems sources.
    """
    need_info = {"need_annotations": True, "need_families": True, "need_families_info": True, "need_graph": True,
                 "need_metadata": True, "metatypes": ["families"], "need_systems": True, "systems_sources": args.sources}
    if not any(arg for arg in [args.projection, args.partition, args.association, args.proksee]):
        raise argparse.ArgumentError(argument=None, message="You should at least choose one type of systems writing "
                                                            "between: projection, partition, association or proksee.")
    if len(args.sources) != len(args.models):
        raise argparse.ArgumentError(argument=None, message="Number of sources and models are different.")
    if args.annotation_sources is not None:
        if len(args.annotation_sources) != len(args.sources):
            raise argparse.ArgumentError(argument=None, message="Number of annotation sources is different from "
                                                                "number of systems sources.")
    else:
        args.annotation_sources = args.sources

    need_info["sources"] = args.annotation_sources

    if args.association is not None:
        for asso in args.association:
            if asso == "modules":
                need_info["need_modules"] = True
            if asso == "RGPs":
                need_info["need_rgp"] = True
            if asso == "spots":
                need_info["need_regions"] = True
                need_info["need_spots"] = True
    return need_info


def check_pangenome_write_systems(pangenome: Pangenome, sources: List[str]) -> None:
    """
     Check and load pangenome information before adding annotation

    Args:
        pangenome:  Pangenome object
        sources: Sources used to detect systems

    Raises:
        KeyError: Provided systems source is not in the pangenome
        Exception: Systems have not been detected in pangenome
        AttributeError: If there is no metadata associated to families
    """
    if pangenome.status["systems"] != "inFile":
        raise AttributeError("Systems have not been detected."
                             "Use 'panorama detect' subcommand to detect systems in pangenomes.")
    else:
        for systems_source in sources:
            if systems_source not in pangenome.status["systems_sources"]:
                logging.getLogger("PANORAMA").error(f"Systems in pangenome {pangenome.name} are: "
                                                    f"{pangenome.status['systems_sources']}")
                raise KeyError(f"There is no systems in pangenome {pangenome.name}, for the source: {systems_source}."
                               f"Look at 'panorama detect' subcommand to detect systems in {pangenome.name}.")


def write_pangenomes_systems(pangenomes: Pangenomes, output: Path, annotation_sources: List[str],
                             projection: bool = False, association: List[str] = None,
                             partition: bool = False, proksee: str = None, organisms: List[str] = None,
                             threads: int = 1, lock: Lock = None, force: bool = False, disable_bar: bool = False):
    """
    Write flat files about systems for all pangenomes

    Args:
        pangenomes: Pangenome objects with all pangenome
        output: Path to write flat files about systems
        projection: Flag to enable/disable pangenome projection (default: False)
        association: write systems association to the given pangenome object (default: False)
        partition: Flag to enable write system partition (default: False)
        proksee: write proksee with the systems and the given pangenome object (default: False)
        organisms: List of organism names to write (default: all organisms)
        threads: Number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        force: Flag to allow overwriting files (default: False)
        disable_bar: Flag to disable the progress bar (default: False)
    """
    pangenomes_proj = pd.DataFrame()
    for pangenome in tqdm(pangenomes, total=len(pangenomes), unit='pangenome', disable=disable_bar):
        logging.getLogger("PANORAMA").debug(f"Begin write systems for {pangenome.name}")
        for system_source in pangenome.systems_sources:
            logging.getLogger("PANORAMA").debug(
                f"Begin write systems for {pangenome.name} on system source {system_source}. Based on annotation sources {annotation_sources} ")
            pangenome_proj, organisms_proj = project_pangenome_systems(pangenome, system_source, annotation_sources,
                                                                       association=association, threads=threads,
                                                                       lock=lock, disable_bar=disable_bar)
            if partition:
                logging.getLogger("PANORAMA").debug(f"Write partition systems for {pangenome.name}")
                systems_partition(pangenome.name, pangenome_proj, output)
            if association:
                logging.getLogger("PANORAMA").debug(f"Write systems association for {pangenome.name}")
                association_pangenome_systems(pangenome, association, output)
            if proksee:
                raise NotImplementedError("Proksee not implemented")
            if projection:
                logging.getLogger("PANORAMA").debug(f"Write projection systems for {pangenome.name}")
                write_projection_systems(pangenome.name, output, system_source, pangenome_proj, organisms_proj,
                                         organisms, force)
            pangenome_proj.insert(0, "pangenome name", pangenome.name)
            pangenomes_proj = pd.concat([pangenomes_proj, pangenome_proj])


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    from panorama.format.read_binaries import load_pangenomes
    from panorama.utility.utility import check_models

    need_info = check_write_systems_args(args)
    models_list = []
    for models in args.models:
        models_list.append(check_models(models, disable_bar=args.disable_prog_bar))

    need_info["models"] = models_list

    outdir = mkdir(args.output, force=args.force)
    manager = Manager()
    lock = manager.Lock()

    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, check_function=check_pangenome_write_systems,
                                 need_info=need_info, sources=args.sources, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar)

    write_pangenomes_systems(pangenomes, outdir, args.annotation_sources, projection=args.projection,
                             proksee=args.proksee, association=args.association, partition=args.partition,
                             organisms=args.organisms, threads=args.threads, lock=lock, force=args.force,
                             disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write_systems", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    required.add_argument('-m', '--models', required=True, type=Path, nargs="+",
                          help="Path to model list file. You can specify multiple models from different source. "
                               "For that separate the model list files by a space and "
                               "make sure you give them in the same order as the sources.")
    required.add_argument("-s", "--sources", required=True, type=str, nargs="+",
                          help="Name of the systems sources. You can specify multiple sources. "
                               "For that separate names by a space and "
                               "make sure you give them in the same order as the sources.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--projection", required=False, action="store_true",
                          help="Project the systems on organisms. If organisms are specified, "
                               "projection will be done only for them.")
    optional.add_argument("--partition", required=False, action="store_true",
                          help="Write a heatmap file with for each organism, partition of the systems. "
                               "If organisms are specified, heatmap will be write only for them.")
    optional.add_argument("--association", required=False, type=str, default=None, nargs='+',
                          choices=["all", "modules", "RGPs", "spots"],
                          help="Write association between systems and others pangenomes elements")
    optional.add_argument("--proksee", required=False, type=str, default=None, nargs='+',
                          choices=["all", "base", "modules", "RGP", "spots", "annotations"],
                          help="Write a proksee file with systems. "
                               "If you only want the systems with genes, gene families and partition, use base value."
                               "Write RGPs, spots or modules -split by `,'- if you want them.")
    required.add_argument("--annotation_sources", required=False, type=str, nargs="+", default=None,
                          help="Name of the annotation sources if different from systems. "
                               "You can specify multiple sources. For that separate names by a space and "
                               "make sure you give them in the same order as the sources.")
    optional.add_argument("--organisms", required=False, type=str, default=None, nargs='+')
    optional.add_argument("--threads", required=False, type=int, default=1)
