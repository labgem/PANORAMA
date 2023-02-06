#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import re
import tempfile
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

    def translate_gene(elem):
        try:
            name = "".join(elem.get('name').split('__')[1:])
        except KeyError:
            Exception(f"In {data_json['name']} of MacSyFinder, one gene doesn't have a name")
        else:
            dict_elem = {'name': name, 'type': 'neutral', 'parameters': {"max_separation": 0}}

            for attrib, value in elem.attrib.items():
                if attrib == 'presence':
                    dict_elem['type'] = value
                elif attrib == 'inter_gene_max_space':
                    dict_elem['parameters']['max_separation'] = int(value)

            for relation in elem:
                if relation.tag == "exchangeables":
                    dict_elem["exangeables"] = []
                    for gene in relation:
                        try:
                            name = "".join(elem.get('name').split('__')[1:])
                        except KeyError:
                            Exception(f"In {data_json['name']} of MacSyFinder, one gene doesn't have a name")
                        else:
                            dict_elem["exangeables"].append(name)
                else:
                    logging.getLogger().warning("Unexpected relation")
            return dict_elem

    def translate_fu(elem):
        try:
            name = elem.get('name')
        except KeyError:
            Exception(f"In {data_json['name']} of MacSyFinder, one gene doesn't have a name")
        else:
            dict_elem = {'name': name, 'type': 'neutral', "families": [],
                         'parameters': {"max_separation": 0, 'max_forbidden': 0,
                                        'min_mandatory': 1, "min_total": 1}}

            for attrib, value in elem.attrib.items():
                if attrib == 'presence':
                    dict_elem['type'] = value
                elif attrib == 'min_mandatory_genes_required':
                    dict_elem['parameters']['min_mandatory'] = int(value)
                elif attrib == 'min_genes_required':
                    dict_elem['parameters']['min_total'] = int(value)
                elif attrib == 'inter_gene_max_space':
                    dict_elem['parameters']['max_separation'] = int(value)

            for family in elem:
                dict_elem["families"].append(translate_gene(family))
            return dict_elem

    if root.attrib is not None:  # Read attributes root
        for parameter, number in root.attrib.items():
            if parameter == 'inter_gene_max_space':
                data_json['parameters']['max_separation'] = int(number)
            elif parameter == 'min_mandatory_genes_required':
                data_json['parameters']['min_mandatory'] = int(number)
            elif parameter == 'min_genes_required':
                data_json['parameters']['min_total'] = int(number)
    data_json['parameters']['max_forbidden'] = 0

    fu_list = []
    fam_list = []
    for elem in root:
        if elem.tag == "gene":
            fam_list.append(translate_gene(elem))

        elif elem.tag == "functional_unit":
            fu_list.append(translate_fu(elem))
    if len(fu_list) == 0:  # only genes
        fu_list.append({'name': data_json["name"], 'type': 'mandatory', "families": fam_list, 'parameters': {}})
    for fu in fu_list:
        data_json["func_units"].append(fu)


def add_fu_field(model_file: Path, tmp: Path):
    with open(model_file, "r") as model:
        lines = model.readlines()
    lines = list(filter(lambda element: element != "\n", lines))
    lines = list(filter(lambda element: element != " \n", lines))
    newfile_path = tmp.joinpath(model_file.name)
    with open(newfile_path, "w") as model_test:
        index = 0
        in_fu = False
        while index < len(lines):
            line = lines[index]
            if re.search("^#########################", line):
                index += 2
            elif re.search("^ *#+ ", line):
                in_fu = True
                fu_name = line.replace("#", "").replace(" ", "").replace("\n", "")
                presence = re.findall(r'presence="(\w+)"', lines[index + 1])[0]
                model_test.write(f'<functional_unit name="{fu_name}" presence="{presence}" '
                                 f'min_mandatory_genes_required="1" min_genes_required="1">\n')
            elif re.search('</gene>', line) and in_fu:
                model_test.write(line)
                if re.search("^#", lines[index + 1]):
                    model_test.write("</functional_unit>\n")
            elif re.search('</model>', line) and in_fu:
                model_test.write("</functional_unit>\n")
                model_test.write(line)
            # elif re.search("<gene.*/>", line) and not in_fu:
            #     fu_name = re.findall(r'name="(\w+)"', line)[0]
            #     presence = re.findall(r'presence="(\w+)"', line)[0]
            #     model_test.write(f"<functional_unit name={fu_name}, presence={presence}, "
            #                      f'min_mandatory_genes_required="1" min_genes_required="1">\n')
            else:
                model_test.write(line)
            index += 1
    return newfile_path


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
        for file in systems_path.rglob("*.xml"):
            new_file = add_fu_field(file, args.tmpdir)
            data_json = {"name": file.stem, 'func_units': [], 'parameters': dict()}
            root = read_xml(new_file)
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
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")


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
