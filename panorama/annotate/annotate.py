#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import csv
import logging
from pathlib import Path
import tempfile

# installed libraries
import pandas as pd

# local libraries
from panorama.utils import check_tsv_sanity
from panorama.annotate.hmm_search import annot_with_hmm, res_col_names
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.format.read_binaries import check_pangenome_info
from annotation import Annotation
from panorama.pangenomes import Pangenome


def check_parameter(args):
    if args.tsv is None and args.hmm is None:
        raise Exception("You did not provide tsv or hmm for annotation")
    if args.hmm is not None and args.meta is None:
        raise Exception("You did not provide metadata file to assign function with HMM")


def check_pangenome_annotation(pangenome: Pangenome, source: str, force: bool = False, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param disable_bar: Disable bar
    """
    try:
        _ = pangenome.status["annotation_source"]
    except KeyError:
        check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)
    else:
        if source in pangenome.status["annotation_source"]:
            if not force:
                raise Exception(f"An annotation corresponding to the source : '{source}' already exist in pangenome."
                                f"Change the source name or add the option --force to erase")
            else:
                erase_pangenome(pangenome, annotation=True, source=source)
        if pangenome.status['annotation'] == 'inFile':  # Source annotation was removed but other sources could be loaded
            check_pangenome_info(pangenome, need_annotations=True, need_families=True,
                                 need_annotation_fam=True, disable_bar=disable_bar)
        else:
            check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)


def annotation_to_families(annotation_df: pd.DataFrame, pangenome: Pangenome, source: str = None,
                           max_prediction: int = None) -> dict:
    """ Add to gene families an annotation and create a dictionary with for each annotation a set of gene family

    :param annotation_df: Dataframe with for each family an annotation
    :param pangenome: Pangenome with gene families
    :param source: source of the annotation
    :param force: force to write the annotation in gene families

    :return: Dictionary with for each annotation a set of gene family
    """
    # annot2fam = {}
    for row in annotation_df.itertuples(index=False):
        gene_fam = pangenome.get_gene_family(name=row.Gene_family)
        if gene_fam is not None:
            annotation = Annotation(source=source, accession=row.Accession, name=row.protein_name,
                                    secondary_names=row.secondary_name, description=row.Description,
                                    score=row.score, e_val=row.e_value, bias=row.bias)
            gene_fam.add_annotation(source=source, annotation=annotation, max_prediction=max_prediction)
        else:
            logging.getLogger().warning(
                f"Family {row['Gene_family']} does not exist in pangenome. If you give a tsv check name."
                f"If you're using hmm annotation, please post an issue on our github.")
            # if search_sys:
            #     for index, row in select_df.iterrows():
            #         if row['protein_name'] not in annot2fam:
            #             annot2fam[row['protein_name']] = {gene_fam}
            #         else:
            #             annot2fam[row['protein_name']].add(gene_fam)
            #         if not pd.isna(row['secondary_name']):  # Test if there is a second name
            #             for other_name in re.split('\|', row['secondary_name']):
            #                 if row['secondary_name'] != '-':
            #                     if other_name not in annot2fam:
            #                         annot2fam[other_name] = {gene_fam}
            #                     else:
            #                         annot2fam[other_name].add(gene_fam)
    pangenome.status["annotation"] = "Computed"
    # return annot2fam


def annot_pangenome(pangenome: Pangenome, hmm: Path, tsv: Path, meta: Path = None, mode: str = 'fast', msa: Path = None,
                    source: str = None, max_prediction: int = 1, tmpdir: Path = Path(tempfile.gettempdir()),
                    threads: int = 1, disable_bar: bool = False) -> dict:
    """ Main function to add annotation to pangenome from tsv file

    :param pangenome: Pangenome object to ppanggolin
    :param hmm: Path to hmm file or directory
    :param tsv: Path to tsv with gene families annotation
    :param meta:
    :param mode:
    :param msa:
    :param source: Source of annotation
    :param prediction_size:
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param force: Boolean of force argument
    :param disable_bar: Disable bar

    :return: Dictionnary with for each annotation a set of corresponding gene families
    """
    if tsv is not None:
        annotation_df = pd.read_csv(tsv, sep="\t", header=None, quoting=csv.QUOTE_NONE, names=res_col_names)
    elif hmm is not None:
        annotation_df = annot_with_hmm(pangenome, hmm, meta=meta, tmpdir=tmpdir, threads=threads, mode=mode, msa=msa,
                                       disable_bar=disable_bar)
    else:
        raise Exception("You did not provide tsv or hmm for annotation")
    return annotation_to_families(annotation_df, pangenome, source, max_prediction)


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    check_parameter(args)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for pangenome_name, pangenome_info in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
        pangenome.add_file(pangenome_info["path"])
        check_pangenome_annotation(pangenome, source=args.source, force=args.force, disable_bar=args.disable_prog_bar)
        annot2fam = annot_pangenome(pangenome=pangenome, hmm=args.hmm, tsv=args.tsv, source=args.source,
                                    meta=args.meta, max_prediction=args.max_prediction, mode=args.mode, msa=args.msa,
                                    tmpdir=args.tmpdir, threads=args.threads, disable_bar=args.disable_prog_bar)
        logging.getLogger().info("Annotation Done")
        logging.getLogger().info(f"Write Annotation in pangenome {pangenome_name}")
        write_pangenome(pangenome, pangenome_info["path"], source=args.source, disable_bar=args.disable_prog_bar)


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
    required.add_argument("--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--tsv', type=Path, nargs='?',
                                help='Gene families annotation in TSV file. See our github for more detail about format')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?',
                                help="File with all HMM or a directory with one HMM by file")
    hmm_param = parser.add_argument_group(title="HMM arguments",
                                          description="All of the following arguments are required,"
                                                      " if you're using HMM mode :")
    hmm_param.add_argument("--meta", required=False, type=Path, default=None,
                           help="Metadata link to HMM with protein name, description and cutoff")
    hmm_param.add_argument("--mode", required=False, type=str, default='fast',  choices=['fast', 'profile'],
                           help="Choose the mode use to align HMM database and gene families. "
                                "Fast will align the reference sequence of gene family against HMM."
                                "Profile will create an HMM profile for each gene family and "
                                "this profile will be aligned")
    hmm_param.add_argument("--msa", required=False, type=Path, default=None,
                           help="To create a HMM profile for families, you can give a msa of each gene in families."
                                "This msa could be get from ppanggolin (See ppanggolin msa). "
                                "If no msa provide Panorama will launch one.")
    hmm_param.add_argument("--msa-format", required=False, type=str, default="afa",
                           choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
                                    "clustal", "clustallike", "phylip", "phylips"],
                           help="Format of the input MSA.")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--max_prediction", required=False, type=int, default=1,
                          help="Maximum number of prediction associate to each gene family")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


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
