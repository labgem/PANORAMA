#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations

import argparse
import logging
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock, Manager
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

# Third-party imports
import numpy as np
import pandas as pd
import ppanggolin.metadata
from tqdm import tqdm

# local libraries
from panorama.annotate.hmm_search import profile_gf
from panorama.format.read_binaries import load_pangenomes
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.utils import init_lock, mkdir


def check_flat_parameters(
    args: argparse.Namespace,
) -> Tuple[Dict[str, Union[bool, Any]], Dict[str, Union[bool, Any]]]:
    """
    Checks if given command argument are legit if so return a dictionary with information needed to load pangenomes.

    Args:
        args: Argument in the command line

    Returns:
        Dictionary needed to load pangenomes information
    """
    if not args.hmm and not args.annotations and not args.conserved_spots:
        raise argparse.ArgumentError(
            argument=None,
            message="You need to provide at least --annotation, --hmm or --conserved_spots",
        )
    else:
        need_info = {}
        kwargs = {}
        if args.annotations:
            if args.sources is None:
                raise argparse.ArgumentError(
                    argument=None,
                    message="You need to provide at least one annotation source to write annotation",
                )
            else:
                need_info.update(
                    {
                        "need_families": True,
                        "need_metadata": True,
                        "metatypes": ["families"],
                        "sources": args.sources,
                    }
                )
                kwargs["sources"] = args.sources
        if args.hmm:
            if args.msa is None or args.msa_format is None:
                raise argparse.ArgumentError(None, "To write HMM you need to give msa files and format")
            else:
                need_info["need_families"] = True
                kwargs["msa"] = args.msa
                kwargs["msa_format"] = args.msa_format

        return need_info, kwargs


def check_pangenome_write_flat_annotations(func):
    """
    Decorator to check pangenome to write annotations

    Args:
        func: Function to decorate

    Returns:
        wrapped function
    """

    def wrapper(pangenome, *args, **kwargs):
        """Wrapper function

        Args:
            pangenome: pangenome to check
            *args: all possible necessary args
            **kwargs: all possible necessary kwargs

        Returns:
            the function wrapped
        """
        if kwargs.get("check_annotations", False):
            if pangenome.status["metadata"]["families"] not in [
                "inFile",
                "Computed",
                "Loaded",
            ]:
                raise ValueError("Pangenome families are not associated to any metadata/annotation")
            else:
                for source in kwargs["sources"]:
                    if source not in pangenome.status["metasources"]["families"]:
                        raise KeyError(
                            f"There is non metadata corresponding to the source : '{source}' "
                            f"in pangenome: {pangenome.name}"
                        )
        return func(pangenome, *args, **kwargs)

    return wrapper


def check_pangenome_write_flat_hmm(func):
    """
    Decorator to check pangenome to write annotations

    Args:
        func: Function to decorate

    Returns:
        wrapped function
    """

    def wrapper(pangenome, *args, **kwargs):
        """Wrapper function

        Args:
            pangenome: pangenome to check
            *args: all possible necessary args
            **kwargs: all possible necessary kwargs

        Returns:
            the function wrapped
        """
        if kwargs.get("check_hmm", False):
            if pangenome.status["geneFamilySequences"] == "No":
                raise AttributeError("There is no sequences associated to the gene families.")
        return func(pangenome, *args, **kwargs)

    return wrapper


@check_pangenome_write_flat_annotations
@check_pangenome_write_flat_hmm
def check_pangenome_write_flat(pangenome: Pangenome, **kwargs) -> None:
    """
    Function to check and load one pangenome in the loading process
    Args:
        pangenome:
    """
    if pangenome.status["genesClustered"] == "No":
        raise AttributeError(
            "You did not cluster the genes. See the 'ppanggolin cluster' command if you want to do that."
        )


def write_pangenome_families_annotations(
    pangenome: Pangenome, output: Path, sources: List[str], disable_bar: bool = False
):
    """
    Write a tsv file with all annotations and sources present in pangenome

    Args:
        pangenome: Pangenome with annotation loaded
        output: Output directory to save the tsv file
        sources: sources to write
        disable_bar: Flag to disable the progress bar (default: False)
    """

    source_column_name = [f"Annotation_{source},Accession_{source},Secondary_names_{source}" for source in sources]
    column_name = np.array(f"Pangenome,families,{','.join(source_column_name)}".split(","))
    array_list = []
    for gf in tqdm(pangenome.gene_families, unit="gene families", disable=disable_bar):
        if any(source in sources for source in gf.sources):
            annot_array = np.empty((gf.max_metadata_by_source()[1], 2 + len(sources) * 3), dtype=object)
            if annot_array.shape[0] > 0:
                annot_array[:, 0] = pangenome.name
                annot_array[:, 1] = gf.name
                index_source = 2
                for source in sources:
                    index_annot = 0
                    if source in gf.sources:
                        for annotation in gf.get_metadata_by_source(source).values():
                            annotation: ppanggolin.metadata.Metadata
                            annot_array[index_annot, index_source] = annotation.protein_name
                            annot_array[index_annot, index_source + 1] = annotation.Accession
                            if "secondary_name" in annotation.__dict__.keys() and (
                                annotation.secondary_name is not None or annotation.secondary_name != pd.NA
                            ):
                                annot_array[index_annot, index_source + 2] = annotation.secondary_name
                            else:
                                annot_array[index_annot, index_source + 2] = "-"
                            index_annot += 1
                    index_source += 3
                array_list.append(annot_array)
    out_df = pd.DataFrame(np.concatenate(array_list), columns=column_name)
    out_df = out_df.sort_values(by=["Pangenome", "families"] + list(column_name[range(3, len(column_name), 2)]))
    out_df.to_csv(
        output / pangenome.name / "families_annotations.tsv",
        sep="\t",
        columns=out_df.columns[1:],
        header=True,
        index=False,
    )


def write_pangenomes_families_annotations(
    pangenomes: Pangenomes,
    output: Path,
    sources: List[str],
    threads: int = 1,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Function to write annotations from multiple pangenomes

    Args:
        pangenomes: Pangenomes object containing all pangenome
        output: Path to the output directory
        sources: List of sources to write annotations for families
        threads: Number of available threads (default = 1)
        lock: Global lock for multiprocessing execution (default: None)
        force: Flag to indicate if a path can be overwritten (default: False)
        disable_bar: Disable progress bar (default: False)
    """

    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit="pangenome", disable=disable_bar) as progress:
            futures = []
            for pangenome in pangenomes:
                if Path(output / pangenome.name / "families_annotations.tsv").exists():
                    if force:
                        logging.getLogger("PANORAMA").warning(
                            "Family annotations files already exist and will be overwritten"
                        )
                    else:
                        raise FileExistsError(
                            "Family annotations files already exist. Please use --force to overwrite them."
                        )
                logging.getLogger("PANORAMA").debug(f"Write annotation for pangenome {pangenome.name}")
                future = executor.submit(
                    write_pangenome_families_annotations,
                    pangenome,
                    output,
                    sources,
                    disable_bar,
                )
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def write_hmm(
    family: GeneFamily,
    output: Path,
    msa_file_path: Path = None,
    msa_format: str = "afa",
):
    """
    Write an HMM profile for a gene family

    Args:
        family: A pangenome gene family
        output: Path to the directory to save the HMM
        msa_file_path: Path to the msa file to compute the HMM
        msa_format: format of the msa file to read it correctly (default: "afa")
    """
    if family.HMM is None:
        if msa_file_path is None:
            raise AssertionError(
                "Your gene family is not associated to a HMM and it could not be computed without a MSA file."
            )
        profile_gf(family, msa_file_path, msa_format=msa_format)
    with open(output / f"{family.name}.hmm", "wb") as hmm_file:
        family.HMM.write(hmm_file)


def write_hmm_profile(
    pangenomes: Pangenomes,
    msa_tsv_path: Path,
    output: Path,
    msa_format: str = "afa",
    threads: int = 1,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Write an HMM profile for all gene families in pangenomes

    Args:
        pangenomes: Pangenomes object with all pangenome
        msa_tsv_path: Path to the tsv file with msa
        output: Path to the output directory
        msa_format: format of the msa file to read it correctly (default: "afa")
        threads: Number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        force: Flag to indicate if a path can be overwritten (default: False)
        disable_bar: Disable progress bar (default: False)
    """
    msa_df = pd.read_csv(
        msa_tsv_path,
        sep="\t",
        names=["path"],
        index_col=0,
    )
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        total_families = sum(pangenome.number_of_gene_families for pangenome in pangenomes)
        with tqdm(
            total=total_families,
            unit="gene families",
            desc="write gene families hmm/profile",
            disable=disable_bar,
        ) as progress:
            futures = []
            for pangenome in pangenomes:
                output_path = mkdir(output / f"{pangenome.name}/HMM", force)
                msa_dir_path = Path(msa_df.loc[pangenome.name].to_list()[0])
                for family in pangenome.gene_families:
                    msa_file_path = msa_dir_path / f"{family.name}.aln"
                    if msa_file_path.exists():
                        future = executor.submit(write_hmm, family, output_path, msa_file_path, msa_format)
                    else:
                        future = executor.submit(write_hmm, family, output_path, None, msa_format)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)

            for future in futures:
                future.result()


def write_flat_files(
    pangenomes: Pangenomes,
    output: Path,
    annotation: bool = False,
    hmm: bool = False,
    threads: int = 1,
    lock: Lock = None,
    force: bool = False,
    disable_bar: bool = False,
    **kwargs,
):
    """
    Global function to write all flat files from pangenomes.

    Args:
        pangenomes: Pangenomes object containing pangenome
        output: Path to the output directory for all flat files
        annotation: Flag to indicate whether annotations should be written or not (default: False)
        hmm: Flag to indicate whether hmms should be written or not (default: False)
        threads: Number of available threads (default: 1)
        lock: Global lock for multiprocessing execution (default: None)
        force: Flag to indicate if a path can be overwritten (default: False)
        disable_bar: Disable progress bar (default: False)
        **kwargs: Additional keyword needed according to expected flat file
    """
    for pangenome in pangenomes:
        mkdir(output / pangenome.name, force=force)
    if annotation:
        assert "sources" in kwargs, "No sources were given to write families annotations."
        write_pangenomes_families_annotations(pangenomes, output, kwargs["sources"], threads, lock, force, disable_bar)
    if hmm:
        write_hmm_profile(
            pangenomes,
            msa_tsv_path=kwargs["msa_tsv_path"],
            msa_format=kwargs["msa_format"],
            output=output,
            threads=threads,
            lock=lock,
            force=force,
            disable_bar=disable_bar,
        )


def launch(args):
    """
    Launch functions to write flat files from pangenomes

    Args:
        args: argument given in CLI
    """
    need_info, kwargs = check_flat_parameters(args)
    flat_outdir = mkdir(args.output, force=args.force)
    manager = Manager()
    lock = manager.Lock()
    pangenomes = load_pangenomes(
        pangenome_list=args.pangenomes,
        need_info=need_info,
        check_function=check_pangenome_write_flat,
        max_workers=args.threads,
        lock=lock,
        disable_bar=args.disable_prog_bar,
        check_annotations=args.annotations,
        sources=args.sources,
        check_hmm=args.hmm,
        check_conserved_spots=args.conserved_spots,
    )
    write_flat_files(
        pangenomes,
        output=flat_outdir,
        annotation=args.annotations,
        hmm=args.hmm,
        threads=args.threads,
        lock=lock,
        force=args.force,
        disable_bar=args.disable_prog_bar,
        **kwargs,
    )


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for annot command

    Returns:
        argparse.ArgumentParser: parser arguments for annot command
    """
    parser = sub_parser.add_parser("write")
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Add argument to parser for write command

    Args:
        parser: parser for write argument
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
    required.add_argument("-o", "--output", required=True, type=Path, nargs="?", help="Output directory")
    annotation = parser.add_argument_group(
        title="Write annotation",
        description="Arguments for writing annotation/metadata assigned to gene families in pangenomes",
    )
    annotation.add_argument(
        "--annotations",
        required=False,
        action="store_true",
        help="Write all the annotations from families for the given sources",
    )
    annotation.add_argument(
        "--sources",
        required=False,
        type=str,
        nargs="+",
        default=None,
        help="Name of the annotation source where panorama as to select in pangenomes",
    )
    hmm = parser.add_argument_group(title="Write HMM", description="Arguments for writing gene families HMM")
    hmm.add_argument(
        "--hmm",
        required=False,
        action="store_true",
        help="Write an hmm for each gene families in pangenomes",
    )
    hmm.add_argument(
        "--msa",
        required=False,
        type=Path,
        default=None,
        help="To create a HMM profile for families, you can give a msa of each gene in families."
        "This msa could be get from ppanggolin (See ppanggolin msa). "
        "Should be a 2 column tsv file with pangenome name in first and path to MSA in second",
    )
    hmm.add_argument(
        "--msa_format",
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
        help="Format of the input MSA.",
    )
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--threads", required=False, type=int, default=1)
