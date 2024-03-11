#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from typing import Dict, Union
from multiprocessing import Manager, Lock
from concurrent.futures import ThreadPoolExecutor

# installed libraries
from tqdm import tqdm
import numpy as np
import ppanggolin.metadata

# local libraries
from panorama.annotate.hmm_search import profile_gfs
from panorama.format.read_binaries import check_pangenome_info, load_pangenomes
from panorama.format.write_proksee import write_proksee
from panorama.utils import check_tsv_sanity, mkdir, init_lock
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenomes
from panorama.format.conserved_spot import *
from panorama.alignment.align import all_against_all


def check_flat_parameters(args: argparse.Namespace) -> Dict[str, Union[bool, List[str]]]:
    """
    Checks if given command argument are legit if so return a dictionary with information needed to load pangenomes.

    Args:
        args: Argument in the command line

    Returns:
        Dictionary needed to load pangenomes information
    """
    # if args.hmm and (args.msa is None or args.msa_format is None):
    #     raise Exception("To write HMM you need to give msa files and format")
    need_info = {}
    if args.annotations is not None:
        if args.sources is None:
            raise argparse.ArgumentError(argument=None, message="You need to provide at least "
                                                                "one annotation source to write annotation")
        else:
            need_info.update({"need_families": True, "need_metadata": True,
                              "metatypes": ["families"], "sources": args.sources})
    return need_info


def check_pangenome_write_flat(pangenome: Pangenome, sources: List[str]) -> None:
    """
    Function to check and load one pangenome in the loading process
    Args:
        pangenome:
        sources:
    """
    if pangenome.status["metadata"]["families"] != "inFile":
        raise ValueError("Pangenome families are not associated to any metadata/annotation")
    else:
        for source in sources:
            if source not in pangenome.status["metasources"]["families"]:
                raise KeyError(f"There is non metadata corresponding to the source : '{source}' "
                               f"in pangenome: {pangenome.name}")


def write_pangenome_families_annotations(pangenome: Pangenome, output: Path, sources: List[str],
                                         disable_bar: bool = False):
    """
    Write a tsv file with all annotations and sources present in pangenome

    Args:
        pangenome: Pangenome with annotation loaded
        output: Output directory to save the tsv file
        sources: sources to write
        disable_bar: Flag to disable the progress bar (default: False)
    """

    source_column_name = [f'Annotation_{source},Accession_{source},Secondary_names_{source}' for source in sources]
    column_name = np.array(f"Pangenome,families,{','.join(source_column_name)}".split(','))
    array_list = []
    for gf in tqdm(pangenome.gene_families, unit='gene families', disable=disable_bar):
        if any(source in sources for source in gf.sources):
            annot_array = np.empty((gf.max_metadata_by_source()[1], 2 + len(sources) * 3), dtype=object)
            if annot_array.shape[0] > 0:
                annot_array[:, 0] = pangenome.name
                annot_array[:, 1] = gf.name
                index_source = 2
                for source in sources:
                    index_annot = 0
                    if source in gf.sources:
                        for annotation in gf.get_metadata_by_source(source):
                            annotation: ppanggolin.metadata.Metadata
                            annot_array[index_annot, index_source] = annotation.protein_name
                            annot_array[index_annot, index_source + 1] = annotation.Accession
                            if ('secondary_name' in annotation.__dict__.keys() and
                                    (annotation.secondary_name is not None or annotation.secondary_name != pd.NA)):
                                annot_array[index_annot, index_source + 2] = annotation.secondary_name
                            else:
                                annot_array[index_annot, index_source + 2] = '-'
                            index_annot += 1
                    index_source += 3
                array_list.append(annot_array)
    out_df = pd.DataFrame(np.concatenate(array_list), columns=column_name)
    out_df = out_df.sort_values(by=['Pangenome', 'families'] + list(column_name[range(3, len(column_name), 2)]))
    out_df.to_csv(output / pangenome.name / "families_annotations.tsv", sep="\t",
                  columns=out_df.columns[1:], header=True, index=False)


def write_pangenomes_families_annotations(pangenomes: Pangenomes, output: Path, sources: List[str], threads: int = 1,
                                          lock: Lock = None, force: bool = False, disable_bar: bool = False):
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
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome in pangenomes:
                if Path(output / pangenome.name / "families_annotations.tsv").exists():
                    if force:
                        logging.getLogger("PANORAMA").warning("Family annotations files already exist "
                                                              "and will be overwritten")
                    else:
                        raise FileExistsError("Family annotations files already exist. "
                                              "Please use --force to overwrite them.")
                logging.getLogger("PANORAMA").debug(f"Write annotation for pangenome {pangenome.name}")
                future = executor.submit(write_pangenome_families_annotations, pangenome, output, sources, disable_bar)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def write_flat_files(pangenomes, output: Path, annotation: bool = False, hmm: bool = False, threads: int = 1,
                     lock: Lock = None, force: bool = False, disable_bar: bool = False, **kwargs):
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
        mkdir(output/pangenome.name, force=force)
    if annotation:
        assert 'sources' in kwargs, "No sources were given to write families annotations."
        write_pangenomes_families_annotations(pangenomes, output, kwargs["sources"], threads, lock, force, disable_bar)
    if hmm:
        raise NotImplementedError("Hidden Markov Model for gene families is not implemented yet.")


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    need_info = check_flat_parameters(args)
    flat_outdir = mkdir(args.output, force=args.force)
    manager = Manager()
    lock = manager.Lock()
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_write_flat, max_workers=args.threads, lock=lock,
                                 disable_bar=args.disable_prog_bar, sources=args.sources)
    write_flat_files(pangenomes, output=flat_outdir, annotation=args.annotations, sources=args.sources,
                     force=args.force, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--annotations", required=False, action="store_true",
                          help="Write all the annotations from families for the given sources")
    optional.add_argument("--sources", required=False, type=str, nargs="+", default=None,
                          help='Name of the annotation source where panorama as to select in pangenomes')
    optional.add_argument("--hmm", required=False, action="store_true",
                          help="Write an hmm for each gene families in pangenomes")
    # optional.add_argument("--msa", required=False, type=Path, default=None,
    #                       help="To create a HMM profile for families, you can give a msa of each gene in families."
    #                            "This msa could be get from ppanggolin (See ppanggolin msa). "
    #                            "If no msa provide Panorama will launch one.")
    # optional.add_argument("--msa_format", required=False, type=str, default="afa",
    #                       choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
    #                                "clustal", "clustallike", "phylip", "phylips"],
    #                       help="Format of the input MSA.")
    optional.add_argument("--threads", required=False, type=int, default=1)
