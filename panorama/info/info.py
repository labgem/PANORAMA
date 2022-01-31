#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from tqdm import tqdm
# installed libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats.readBinaries import checkPangenomeInfo
from multiprocessing import get_context


# def mp_info(pangenome_file: str):
#     pass


def loop_info(pangenomes: list, need_info: dict, disable_bar: bool = False):
    for pangenome_file in pangenomes:
        pangenome = Pangenome()
        pangenome.addFile(pangenome_file)
        checkPangenomeInfo(pangenome=pangenome, **need_info, disable_bar=disable_bar)


def check_info(args) -> dict:
    check_dict = {"needAnnotations": False,
                  "needFamilies": False,
                  "needGraph": False,
                  "needPartitions": False,
                  "needRGP": False,
                  "needSpots": False,
                  "needGeneSequences": False,
                  "needModules": False
                  }
    if args.modules is not None:
        check_dict['needFamilies'] = True
        check_dict['needModules'] = True

    return check_dict


def launch(args: argparse.Namespace):
    logging.getLogger().debug("launch info command")
    if len(args.pangenomes) == 0:
        raise Exception("Not one pangenome was found. Please check path to pangenome files.")
    if not any(arg for arg in [args.modules]):
        raise Exception("You did not indicate which information you want.")
    need_info = check_info(args)
    if args.cpu > 1 and len(args.pangenomes) > 1:
        with get_context('fork').Pool(args.cpu) as p:
            for pangenome_file in tqdm(p.imap_unordered(mp_info, args.pangenomes), unit='pangenome',
                                       total=len(pangenomes), disable=args.disable_prog_bar):
                logging.getLogger().debug(f"{pangenome_file} Done")
    else:
        loop_info(args.pangenomes, need_info, args.disable_prog_bar)

    logging.getLogger().info("All the genomes fluidity were computed")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Parser arguments specific to fluidity command
    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """

    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=str, nargs='+',
                          help="A list of pangenome .h5 files")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('--modules', required=False, action="store_true",
                        help="Compute the pangenome genomic fluidity")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")
    return parser
