#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import re

import yaml
import json
from pathlib import Path
from lxml import etree as et

# installed libraries
import pandas as pd


# local libraries


def read_yaml(system_file):
    """
    Open yaml file and return data file

    :param system_file: yaml file
    :return: data : data yaml file
    """
    with open(system_file, "r") as stream:
        data = yaml.safe_load(stream)
        return data


def translate_yaml(data_yaml: dict, data_json: dict, meta: pd.DataFrame = None, canonical: list = None):
    """
    Translate the yaml data in json dictionary

    :param data_yaml: data yaml file to translate
    :param data_json: json dictionary translated
    """
    data_json['func_units'] = [{'name': data_json['name'], 'type': 'mandatory', 'parameters': {'max_forbidden': 0}}]
    data_json['parameters'] = {"max_forbidden": 0, "max_separation": 1, "min_mandatory": 1, "min_total": 1}
    if canonical is not None:
        data_json['canonical'] = canonical
    secondary_names = meta['secondary_name'].dropna().unique().tolist()
    family_list = list()
    for key, value in data_yaml.items():
        if key == 'maximum_separation':
            data_json['func_units'][0]['parameters']['max_separation'] = value
        elif key == 'minimum_core':
            data_json['func_units'][0]['parameters']['min_mandatory'] = value
        elif key == 'minimum_total':
            data_json['func_units'][0]['parameters']['min_total'] = value
        elif key == 'core_genes':
            for fam_name in value:
                if not (fam_name == 'NA'):
                    if fam_name in secondary_names:
                        filter_df = meta.loc[meta['secondary_name'] == fam_name]
                        prime_names = filter_df['protein_name'].dropna().unique().tolist()
                        for prime in prime_names:
                            family_list.append({'name': prime,
                                                'type': 'mandatory',
                                                'relation': None,
                                                'parameters': None})
                    family_list.append({'name': fam_name,
                                        'type': 'mandatory',
                                        'relation': None,
                                        'parameters': None})
        elif key == 'optional_genes':
            for fam_name in value:
                if not (fam_name == 'NA'):
                    if fam_name in secondary_names:
                        filter_df = meta.loc[meta['secondary_name'] == fam_name]
                        prime_names = filter_df['protein_name'].dropna().unique().tolist()
                        for prime in prime_names:
                            family_list.append({'name': prime,
                                                'type': 'accessory',
                                                'relation': None,
                                                'parameters': None})
                    family_list.append({'name': fam_name,
                                        'type': 'accessory',
                                        'relation': None,
                                        'parameters': None})
        elif key == 'prohibited_genes':
            for fam_name in value:
                if not (fam_name == 'NA'):
                    if fam_name in secondary_names:
                        filter_df = meta.loc[meta['secondary_name'] == fam_name]
                        prime_names = filter_df['protein_name'].dropna().unique().tolist()
                        for prime in prime_names:
                            family_list.append({'name': prime,
                                                'type': 'forbidden',
                                                'relation': None,
                                                'parameters': None})
                    family_list.append({'name': fam_name,
                                        'type': 'forbidden',
                                        'relation': None,
                                        'parameters': None})
    data_json['func_units'][0]['families'] = family_list
    # list_fam = list()
    # for k, v in data_json['func_units'].items():
    #     for fam in v["families"].keys():
    #         list_fam.append(fam)
    #     data_json['func_units'][list_fam[0]] = data_json['func_units'].pop(k)


def read_xml(system_file):
    """
    Open xml file and return the data root

    :param system_file: xml file
    :return: data root
    """
    my_tree = et.parse(system_file.absolute().as_posix())
    my_root = my_tree.getroot()
    return my_root


def translate_xml(root, data_json: dict):
    """
    Translate the xml data in json dictionary

    :param root: root of xml data
    :param data_json: json dictionary translated
    :param name_file: name of xml file
    """
    if root.attrib is not None:  # Read attributes root
        for parameter, number in root.attrib.items():
            if parameter == 'inter_gene_max_space':
                data_json['parameters']['max_separation'] = int(number)
            elif parameter == 'min_mandatory_genes_required':
                data_json['parameters']['min_mandatory'] = int(number)
            elif parameter == 'min_genes_required':
                data_json['parameters']['min_total'] = int(number)
    data_json['parameters']['max_forbidden'] = 0
    for fu_field in root:
        try:
            name_fu = "".join(fu_field.get('name').split('__')[1:])
        except KeyError:
            Exception(f"In {data_json['name']} of MacSyFinder, one gene doesn't have a name")
        else:
            dict_fu = {'name': name_fu, 'type': str(), 'families': [],
                       'parameters': {'max_forbidden': 0, 'min_mandatory': 1}}
            for attrib, value in fu_field.attrib.items():
                if attrib == 'presence':
                    dict_fu['type'] = value
                elif attrib == 'inter_gene_max_space':
                    dict_fu['parameters']['max_separation'] = int(value)
            if len(dict_fu['type']) == 0:  # No type
                dict_fu['type'] = "neutral"
            # Write FU as Fam
            dict_fu["families"].append({'name': name_fu,
                                        'type': "mandatory",
                                        'relation': None,
                                        'parameters': {'max_separation': dict_fu['parameters']['max_separation']}
                                        if 'max_separation' in dict_fu['parameters'] else {}}
                                       )
            if len(list(fu_field)) > 0:
                families = list(fu_field)[0]
                for fam_field in families:
                    fam_dict = {'name': "".join(fam_field.get('name').split('__')[1:]),
                                'type': "mandatory",
                                'relation': None,
                                'parameters': {}}
                    if fam_field.tag in ['homologs', 'analogs']:
                        relation = fam_field.tag
                    else:
                        relation = None
                    fam_dict['relation'] = relation
                    for attribute, string in fam_field.attrib.items():
                        if attribute == 'presence':
                            dict_fam['type'] = string
                        if attribute == 'inter_gene_max_space':
                            dict_fam['parameters']['max_separation'] = int(string)
                    dict_fu['families'].append(fam_dict)
            data_json['func_units'].append(dict_fu)


def write(path_output: Path, name: str, dict_data: dict):
    """
    Create new json file and write json data dictionary

    :param path_output: path of output directory
    :param name: name of file to translate
    :param dict_data: json dictionary translated
    """
    with open(f"{path_output}/{name}.json", "w") as my_file:
        json.dump(dict_data, my_file, indent=2)


def launch(args: argparse.Namespace):
    """
    Launch functions with arguments

    :param args: arguments given
    """
    systems_path = Path(args.directory_input)
    translate_path = Path(args.directory_output)
    if args.format == "yaml":
        meta_df = pd.read_csv(filepath_or_buffer=args.meta.absolute().as_posix(), sep='\t', header=0)
        for file in systems_path.glob("*.yaml"):
            data_json = {"name": str(), 'func_units': dict(), 'parameters': dict()}
            canonical_sys = None
            name_file = file.stem
            if re.search("_other", name_file):
                basename = re.split("_other", name_file)[0]
                for canon_file in systems_path.glob("*.yaml"):
                    name_canon = canon_file.stem
                    if (re.search(f"{basename}_", name_canon) or re.search(f"{basename}$", name_canon)) \
                            and name_canon != name_file:
                        if canonical_sys is None:
                            canonical_sys = [name_canon]
                        else:
                            canonical_sys.append(name_canon)
            data = read_yaml(file)
            data_json['name'] = name_file
            translate_yaml(data, data_json, meta_df, canonical_sys)
            write(translate_path, name_file, data_json)
    elif args.format == "xml":
        for file in systems_path.glob("*.xml"):
            data_json = {"name": file.stem, 'func_units': [], 'parameters': dict()}
            root = read_xml(file)
            translate_xml(root, data_json)
            write(translate_path, file.stem, data_json)
    else:
        raise ValueError("Format must be yaml or xml")


def parse_annot(parser: argparse.ArgumentParser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument("-di", "--directory_input",
                          required=True,
                          type=str,
                          help="Name to rule directory to translate")
    required.add_argument("-do", "--directory_output",
                          required=True,
                          type=str,
                          help="Name to output rule directory")
    required.add_argument("--format",
                          required=True,
                          type=str,
                          help="Format to translate in Json")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional :")
    optional.add_argument("--meta", required=False, type=Path,
                          help='metadata')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="comparative pangenomic analyses toolsbox",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parse_annot(parser)
    args = parser.parse_args()
    launch(args)
    # data_json = {}
    # name_file = Path("padloc-db/sys/AbiE.yaml").stem
    # data = read_yaml("padloc-db/sys/AbiE.yaml")
    # data_json['name'] = name_file
    # translate_yaml(data, data_json)
    # print(json.dumps(data_json, indent=4))
