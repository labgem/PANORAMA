#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import json
from pathlib import Path

# installed libraries

# local libraries
from panorama.annotation.rules import Rule, Rules


def fill_attr(content_file: dict, rule: Rule = Rule()) -> Rule:
    """
    From a rule file, give values in attributes of object Rule, calculate and print parameters not respected

    :param content_file: dictionary of file rule content
    :param rule : json file rule in object Rule

    :return: rule: json file rule
    """
    rule.fill_attr(content_file)
    return rule


def read_rules_directory(rules_dir, rules=Rules()):
    """
    Create directory path, read all rules in the directory, store rules respected in dictionary attribute of object
    Rules and print all rules in the dictionary

    :param rules_dir: directory path
    :param rules: all rules respected in object Rules

    :return: rules: all rules respected
    """
    rules_path = Path(rules_dir)
    for rule_file in rules_path.glob("*.json"):  # Find all file in the directory path
        json_dict = check_rule(rule_file)
        rule = fill_attr(json_dict)
        rules.add_rule(rule)
    rules.show_rules()
    return rules


def check_rule(json_file: Path):
    """
    Check conditions of content rule file

    :param json_file: path of the json_file rule

    :return: data: dictionary of file rule content
    """
    with open(json_file) as my_file:
        data = json.load(my_file)
    f_keys = ['families', 'parameters']
    for key in f_keys:
        try:
            sub_dict = data[key]
        except KeyError:
            raise KeyError(f"Requires {key} in {json_file.stem}")
        else:
            if key == 'families':
                fam_keys = ['mandatory', 'accessory', 'forbidden']
                for f_key in fam_keys:
                    try:
                        _ = sub_dict[f_key]
                    except KeyError:
                        raise KeyError(f"Requires {f_key} of families in {json_file.stem}")
            elif key == 'parameters':
                param_keys = ['min_mandatory', 'min_total', 'max_forbidden', 'max_free']
                for p_key in param_keys:
                    try:
                        _ = sub_dict[p_key]
                    except KeyError:
                        raise KeyError(f"Requires {p_key} of parameters in {json_file.stem}")
    data['name'] = json_file.stem
    return data


def launch(args: argparse.Namespace):
    """
    Launch function with arguments

    :param args: argument given
    """
    rules = read_rules_directory(args.rules)
    attr_dict = {'mandatory': args.mandatory, 'accessory': args.accessory,
                 'forbidden': args.forbidden, 'free': args.free}
    if not any(x for x in attr_dict.values()):
        attr_dict = {key: True for key in attr_dict.keys()}


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

    required.add_argument("-r", "--rules", required=True, type=str, help="Path to rules directory")
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
    optional.add_argument('--free',
                          required=False,
                          help="To read values free",
                          action="store_true")


if __name__ == "__main__":
    # from panorama.main import check_log
    # definition du parser
    parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser_annot(parser)
    args = parser.parse_args()
    launch(args)
