#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import json
from pathlib import Path

# installed libraries

# local libraries
from panorama.annotation.rules import System, Systems


def check_rule(json_file: Path):
    """
    Check conditions of content rule file

    :param json_file: path of the json_file rule

    :return: data: dictionary of file rule content
    """
    with open(json_file) as my_file:
        data = json.load(my_file)
    try:
        sub_dict = data["func_unit"]
    except KeyError:
        raise KeyError(f"Requires func_unit in {json_file.stem}")
    else:
        f_keys = ['families', 'parameters']
        for key in f_keys:
            try:
                s_dict = sub_dict[key]
            except KeyError:
                raise KeyError(f"Requires {key} in {json_file.stem}")
            else:
                if key == 'families':
                    fam_keys = ['mandatory', 'accessory', 'forbidden']
                    for f_key in fam_keys:
                        try:
                            _ = s_dict[f_key]
                        except KeyError:
                            raise KeyError(f"Requires {f_key} of families in {json_file.stem}")
                elif key == 'parameters':
                    param_keys = ['min_mandatory', 'min_total', 'max_separation', 'max_forbidden']
                    for p_key in param_keys:
                        try:
                            _ = s_dict[p_key]
                        except KeyError:
                            raise KeyError(f"Requires {p_key} of parameters in {json_file.stem}")
    data['name'] = json_file.stem
    return data


def read_files(systems_path, systems=Systems()):
    for file in systems_path.glob("*.json"):
        with open(file.resolve().as_posix()) as json_file:
            system = System()
            data = json.load(json_file)
            system.read_system(data)
            system.check_param_fu()
            systems.add_sys(system)
    systems.print_systems()


def launch(args: argparse.Namespace):
    """
    Launch function with arguments

    :param args: argument given
    """
    systems_path = Path(args.systems)
    read_files(systems_path)


def subparser(sub_parser: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """
    Parser arguments specific to annotation command

    :param sub_parser : sub_parser for annot command

    :return : parser arguments for annot command
    """
    parser = sub_parser.add_parser("annot", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_annot(parser)
    return parser


def parser_annot(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")

    required.add_argument("-s", "--systems", required=True, type=str, help="Path to rules directory")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument('--mandatory',
                          required=False,
                          help="To read values mandatory",
                          action="store_true")
    optional.add_argument('--accessory',
                          required=False,
                          help="To read values accessory",
                          action="store_true")
    optional.add_argument('--forbidden',
                          required=False,
                          help="To read values forbidden",
                          action="store_true")
    optional.add_argument('--neutral',
                          required=False,
                          help="To read values neutral",
                          action="store_true")


if __name__ == "__main__":
    # from panorama.main import check_log
    # definition du parser
    parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser_annot(parser)
    args = parser.parse_args()
    launch(args)
