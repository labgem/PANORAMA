#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse

def subparser_info(sub_parser):
    """
    Parser arguments specific to metrics command

    :param sub_parser : sub_parser for align command
    :type sub_parser : argparse._SubParsersAction

    :return : parser arguments for align command
    :rtype : argparse.ArgumentParser
    """
    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=str, help="list of all pangenome")
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument('--fluidity', required=False, action="store_true",
                        help="Compute the pangenome genomic fluidity")
    return parser
