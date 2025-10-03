#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to write information into the pangenome file
"""

# default libraries
from __future__ import annotations
import argparse
import logging
import time
from typing import Any, Dict, List
from multiprocessing import Manager, Lock
from pathlib import Path

# installed libraries
from tqdm import tqdm

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.systems.systems_projection import (
    project_pangenome_systems,
    write_projection_systems,
)
from panorama.systems.systems_partitions import systems_partition
from panorama.systems.systems_association import create_pangenome_system_associations


def check_write_systems_args(args: argparse.Namespace) -> Dict[str, Any]:
    """Checks the provided arguments to ensure that they are valid.

    Args:
        args (argparse.Namespace): The parsed arguments.

    Returns:
        Dict[str, Any]: A dictionary containing the necessary information for further processing.

    Raises:
        argparse.ArgumentTypeError: If the number of sources is different from models, or if annotations are given, and
        their number is different from systems sources.
    """
    need_info = {
        "need_annotations": True,
        "need_families": True,
        "need_families_info": True,
        "need_graph": True,
        "need_metadata": True,
        "metatypes": ["families"],
        "need_systems": True,
        "systems_sources": args.sources,
        "read_canonical": args.canonical,
    }
    if not any([args.projection, args.partition, args.association, args.proksee]):
        raise argparse.ArgumentError(
            argument=None,
            message="You should at least choose one type of systems writing "
            "between: projection, partition, association or proksee.",
        )
    if len(args.sources) != len(args.models):
        raise argparse.ArgumentError(
            argument=None, message="Number of sources and models are different."
        )

    if "all" in args.association:
        args.association = ["RGPs", "spots", "modules"]

    for asso in args.association:
        if asso == "modules":
            need_info["need_modules"] = True
        if asso in ["RGPs", "spots"]:
            need_info["need_rgp"] = True
        if asso == "spots":
            need_info["need_spots"] = True
    return need_info


def check_pangenome_write_systems(pangenome: Pangenome, sources: List[str]) -> None:
    """
    Check and load pangenome information before adding annotation.

    Args:
        pangenome (Pangenome): The Pangenome object.
        sources (List[str]): Sources used to detect systems.

    Raises:
        KeyError: If the provided systems source is not in the pangenome.
        Exception: If systems have not been detected in pangenome.
        AttributeError: If there is no metadata associated with families.
    """
    if pangenome.status["systems"] != "inFile":
        raise AttributeError(
            "Systems have not been detected."
            "Use 'panorama systems' subcommand to detect systems in pangenomes."
        )
    else:
        for systems_source in sources:
            if systems_source not in pangenome.status["systems_sources"]:
                logging.getLogger("PANORAMA").error(
                    f"Systems in pangenome {pangenome.name} are: "
                    f"{pangenome.status['systems_sources']}"
                )
                raise KeyError(
                    f"There is no systems in pangenome {pangenome.name}, for the source: {systems_source}."
                    f"Look at 'panorama systems' subcommand to detect systems in {pangenome.name}."
                )


def write_flat_systems_to_pangenome(
    pangenome: Pangenome,
    output: Path,
    projection: bool = False,
    association: List[str] = None,
    partition: bool = False,
    proksee: str = None,
    output_formats: List[str] = None,
    organisms: List[str] = None,
    canonical: bool = False,
    threads: int = 1,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Write detected systems from a pangenome to an output directory in a flat format.

    Args:
        pangenome (Pangenome): The pangenome object containing the detected systems.
        output (Path): The directory where the systems will be written.
        projection (bool, optional): If True, write projection systems. Defaults to False.
        association (List[str], optional): List of associations to be considered. Defaults to None.
        partition (bool, optional): If True, write partition systems. Defaults to False.
        proksee (str, optional): A placeholder for future Proksee integration. Defaults to None.
        output_formats (List[str]): A list of output formats for visualization. Defaults to None.
        organisms (List[str], optional): List of organisms to be considered for projection. Defaults to None.
        canonical (bool, optional): If True, write the canonical version of systems too. Defaults to False.
        threads (int, optional): Number of threads to use for parallel processing. Defaults to 1.
        lock (Lock, optional): A multiprocessing lock to synchronize access. Defaults to None.
        force (bool, optional): If True, overwrite existing files. Defaults to False.
        disable_bar (bool, optional): If True, disable the progress bar. Defaults to False.

    Raises:
        NotImplementedError: If Proksee integration is requested but not implemented.

    """
    logging.getLogger("PANORAMA").info(f"Begin write systems for {pangenome.name}")
    begin = time.time()
    pangenome_res_output = mkdir(output / f"{pangenome.name}", force=force)
    for system_source in pangenome.systems_sources:
        logging.getLogger("PANORAMA").debug(
            f"Begin write systems for {pangenome.name} "
            f"on system source: {system_source}"
        )
        pangenome_proj, organisms_proj = project_pangenome_systems(
            pangenome,
            system_source,
            association=association,
            canonical=canonical,
            threads=threads,
            lock=lock,
            disable_bar=disable_bar,
        )
        source_res_output = mkdir(
            pangenome_res_output / f"{system_source}", force=force
        )
        if projection:
            logging.getLogger("PANORAMA").debug(
                f"Write projection systems for {pangenome.name}"
            )
            write_projection_systems(
                source_res_output,
                pangenome_proj,
                organisms_proj,
                organisms,
                threads,
                force,
                disable_bar,
            )
        if partition:
            logging.getLogger("PANORAMA").debug(
                f"Write partition systems for {pangenome.name}"
            )
            systems_partition(pangenome.name, pangenome_proj, source_res_output)
        if association:
            logging.getLogger("PANORAMA").debug(
                f"Write systems association for {pangenome.name}"
            )
            create_pangenome_system_associations(
                pangenome,
                association,
                source_res_output,
                output_formats=output_formats,
                threads=threads,
                disable_bar=disable_bar,
            )
        if proksee:
            raise NotImplementedError("Proksee not implemented")
    logging.getLogger("PANORAMA").info(
        f"Done write system for {pangenome.name} in {time.time() - begin:2f} seconds"
    )


def write_pangenomes_systems(
    pangenomes: Pangenomes,
    output: Path,
    projection: bool = False,
    association: List[str] = None,
    partition: bool = False,
    proksee: str = None,
    output_formats: List[str] = None,
    organisms: List[str] = None,
    canonical: bool = False,
    threads: int = 1,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Write flat files about systems for all pangenomes.

    Args:
        pangenomes (Pangenomes): Pangenome objects with all pangenome.
        output (Path): Path to write flat files about systems.
        projection (bool, optional): Flag to enable/disable pangenome projection. Defaults to False.
        association (List[str], optional): Write systems association to the given pangenome object. Defaults to None.
        partition (bool, optional): Flag to enable write system partition. Defaults to False.
        proksee (str, optional): Write proksee with the systems and the given pangenome object. Defaults to None.
        output_formats (List[str], optional): List of output formats for visualization. Defaults to None.
        organisms (List[str], optional): List of organism names to write. Defaults to all organisms.
        canonical (bool, optional): If True, write the canonical version of systems too. Defaults to False.
        threads (int, optional): Number of available threads. Defaults to 1.
        lock (Lock, optional): Global lock for multiprocessing execution. Defaults to None.
        force (bool, optional): Flag to allow overwriting files. Defaults to False.
        disable_bar (bool, optional): Flag to disable the progress bar. Defaults to False.
    """
    t0 = time.time()
    for pangenome in tqdm(
        pangenomes, total=len(pangenomes), unit="pangenome", disable=disable_bar
    ):
        write_flat_systems_to_pangenome(
            pangenome,
            output,
            projection,
            association,
            partition,
            proksee,
            output_formats,
            organisms,
            canonical,
            threads,
            lock,
            force,
            disable_bar,
        )
    logging.getLogger("PANORAMA").info(
        f"Done write system for all pangenomes in {time.time() - t0:2f} seconds"
    )


def launch(args):
    """
    Launch functions to read systems.

    Args:
        args: Argument given.
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

    pangenomes = load_pangenomes(
        pangenome_list=args.pangenomes,
        check_function=check_pangenome_write_systems,
        need_info=need_info,
        sources=args.sources,
        max_workers=args.threads,
        lock=lock,
        disable_bar=args.disable_prog_bar,
    )

    write_pangenomes_systems(
        pangenomes,
        outdir,
        projection=args.projection,
        proksee=args.proksee,
        association=args.association,
        partition=args.partition,
        output_formats=args.output_formats,
        organisms=args.organisms,
        canonical=args.canonical,
        threads=args.threads,
        lock=lock,
        force=args.force,
        disable_bar=args.disable_prog_bar,
    )


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in the Command line.

    Args:
        sub_parser (argparse.ArgumentParser): Subparser for align command.

    Returns:
        argparse.ArgumentParser: Parser arguments for align command.
    """
    parser = sub_parser.add_parser("write_systems")
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Parser for specific arguments of the write_systems command.

    Args:
        parser (argparse.ArgumentParser): Parser for annot argument.
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
        help="A list of pangenome .h5 files in .tsv file",
    )
    required.add_argument(
        "-o", "--output", required=True, type=Path, nargs="?", help="Output directory"
    )
    required.add_argument(
        "-m",
        "--models",
        required=True,
        type=Path,
        nargs="+",
        help="Path to model list file. You can specify multiple models from different source. "
        "For that separate the model list files by a space and "
        "make sure you give them in the same order as the sources.",
    )
    required.add_argument(
        "-s",
        "--sources",
        required=True,
        type=str,
        nargs="+",
        help="Name of the systems sources. You can specify multiple sources. "
        "For that separate names by a space and "
        "make sure you give them in the same order as the sources.",
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument(
        "--projection",
        required=False,
        action="store_true",
        help="Project the systems on organisms. If organisms are specified, "
        "projection will be done only for them.",
    )
    optional.add_argument(
        "--partition",
        required=False,
        action="store_true",
        help="Write a heatmap file with for each organism, partition of the systems. "
        "If organisms are specified, heatmap will be write only for them.",
    )
    optional.add_argument(
        "--association",
        required=False,
        type=str,
        default=[],
        nargs="+",
        choices=["all", "modules", "RGPs", "spots"],
        help="Write association between systems and others pangenomes elements",
    )
    optional.add_argument(
        "--proksee",
        required=False,
        type=str,
        default=None,
        nargs="+",
        choices=["all", "base", "modules", "RGP", "spots", "annotations"],
        help="Write a proksee file with systems. "
        "If you only want the systems with genes, gene families and partition, use base value."
        "Write RGPs, spots or modules -split by `,'- if you want them.",
    )
    optional.add_argument(
        "--canonical",
        required=False,
        action="store_true",
        help="Write the canonical version of systems too.",
    )
    optional.add_argument(
        "--organisms",
        required=False,
        type=str,
        default=None,
        nargs="+",
        help="List of organisms to write. If not specified, all organisms will be written.",
    )
    optional.add_argument(
        "--output_formats",
        required=False,
        choices=["html", "png"],
        nargs="+",
        default=["html"],
        help="Visualization output format customization.",
    )
    optional.add_argument("--threads", required=False, type=int, default=1)
