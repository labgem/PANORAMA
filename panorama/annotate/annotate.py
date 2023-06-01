#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from pathlib import Path
import tempfile

# installed libraries
from ppanggolin.meta.meta import check_metadata_format, assign_metadata

# local libraries
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.format.read_binaries import check_pangenome_info
from panorama.utils import check_tsv_sanity
from panorama.annotate.hmm_search import annot_with_hmm
from panorama.pangenomes import Pangenome


def check_pangenome_annotation(pangenome: Pangenome, source: str, force: bool = False, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: source of the annotation
    :param force: erase if an annotation for the provide source already exist
    :param disable_bar: Disable bar
    """
    if pangenome.status["metadata"]["families"] == "inFile" and source in pangenome.status["metasources"]["families"]:
        if force:
            erase_pangenome(pangenome, metadata=True, source=source)
        else:
            raise Exception(f"A metadata corresponding to the source : '{source}' already exist in pangenome."
                            f"Add the option --force to erase")
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=disable_bar)


def annot_pangenome(pangenome: Pangenome, hmm: Path, tsv: Path, meta: Path = None, mode: str = 'fast', msa: Path = None,
                    source: str = None, tmpdir: Path = Path(tempfile.gettempdir()),
                    threads: int = 1, disable_bar: bool = False):
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
    if tsv is not None:
        metadata_df = check_metadata_format(tsv, "families")
        # TODO check column name
    elif hmm is not None:
        metadata_df = annot_with_hmm(pangenome, hmm, meta=meta, tmpdir=tmpdir, threads=threads, mode=mode, msa=msa,
                                     disable_bar=disable_bar)
    else:
        raise Exception("You did not provide tsv or hmm for annotation")
    assign_metadata(metadata_df, pangenome, source, "families", omit=True, disable_bar=disable_bar)


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    # check_parameter(args)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for pangenome_name, pangenome_info in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
        pangenome.add_file(pangenome_info["path"])
        check_pangenome_annotation(pangenome, source=args.source, force=args.force, disable_bar=args.disable_prog_bar)
        annot_pangenome(pangenome=pangenome, hmm=args.hmm, tsv=args.tsv, source=args.source,
                        meta=args.meta, mode=args.mode, msa=args.msa,
                        tmpdir=args.tmpdir, threads=args.threads, disable_bar=args.disable_prog_bar)
        logging.info("Metadata assignation Done")
        logging.info(f"Write Metadata in pangenome {pangenome_name}")
        write_pangenome(pangenome, pangenome_info["path"], force=args.force, disable_bar=args.disable_prog_bar)


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
    exclusive_mode.add_argument('--tsv', type=Path, nargs='?',
                                help='Gene families annotation in TSV file. See our github for more detail about format')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?',
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
                           #      "This msa could be get from ppanggolin (See ppanggolin msa). "
                           #      "If no msa provide Panorama will launch one.")
    hmm_param.add_argument("--msa-format", required=False, type=str, default="afa",
                           choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
                                    "clustal", "clustallike", "phylip", "phylips"],
                           help=argparse.SUPPRESS)
                           # help="Format of the input MSA.")
    optional = parser.add_argument_group(title="Optional arguments")
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
