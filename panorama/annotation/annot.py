#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import sys
import json
from rules import Rule


# Read json file given and print the rule
def read_json(json_file: str):
    first_rule = Rule()
    first_rule.read_json(json_file)
    first_rule.print_rule()


def launch(args):
    # print les arguments required et les optionals
    if args.verbose:
        print("Verbose")
    print(args.json_file)
    read_json(args.json_file)


# Define parser arguments
def parser_annot(parser):
    # required arguments
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-j', '--json_file',
                          required=True,
                          type=str,
                          help="to read json file with rule",
                          nargs='?')
    # optional arguments
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument('-v', '--verbose',
                          required=False,
                          help="verbose",
                          action="store_true")


if __name__ == "__main__":
    # from panorama.main import check_log
    # definition du parser
    parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser_annot(parser)
    args = parser.parse_args()
    launch(args)
