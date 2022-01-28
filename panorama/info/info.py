# coding:utf-8

# default libraries
import argparse
import logging


def launch(args):
    logging.getLogger().debug("launch info command")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Parser arguments specific to fluidity command
    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """

    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=str, help="list of all pangenome")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('--modules', required=False, action="store_true",
                        help="Compute the pangenome genomic fluidity")
    return parser
