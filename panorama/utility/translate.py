#!/usr/bin/env python3
# coding:utf-8
import logging
# default libraries
import re
import tempfile
from typing import List
import yaml
import json
from pathlib import Path
from lxml import etree as et

# installed libraries
import pandas as pd
from tqdm import tqdm

# local libraries


def read_yaml(model):
    """
    Open yaml file and return data file

    :param model: yaml file
    :return: data : data yaml file
    """
    try:
        model_file = open(model, "r")
    except IOError as ioerror:
        raise IOError(f"Problem to open {model}. The foolowing error is direct cause\n{ioerror}")
    except Exception as error:
        raise Exception(f"Unexepected Problem to read {model}. The foolowing error is direct cause\n{error}")
    else:
        data = yaml.safe_load(model_file)
        model_file.close()
        return data


def read_xml(model):
    """
    Open xml file and return the data root

    :param model: xml file
    :return: data root
    """
    try:
        my_tree = et.parse(model.absolute().as_posix())
    except IOError as ioerror:
        raise IOError(f"Problem to open {model}. The foolowing error is direct cause\n{ioerror}")
    except Exception as error:
        raise Exception(f"Unexepected Problem to read {model}. The foolowing error is direct cause\n{error}")
    else:
        my_root = my_tree.getroot()
        return my_root


def write(path_output: Path, dict_data: dict):
    """
    Create new json file and write json data dictionary

    :param path_output: path of output directory
    :param dict_data: json dictionary translated
    """
    with open(f"{path_output}/{dict_data['name']}.json", "w") as my_file:
        json.dump(dict_data, my_file, indent=2)


def parse_meta_padloc(meta: Path, output: Path) -> pd.DataFrame:
    meta_col_names = ["accession", "hmm_name", "protein_name", "secondary_name", "score_threshold",
                      "eval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
    df = pd.read_csv(meta, sep="\t", usecols=[0, 1, 2, 3, 4, 8, 9, 10], header=0)
    df[['protein_name', 'temp']] = df['protein.name'].str.split('|', expand=True)
    df['temp'].fillna('', inplace=True)
    df['secondary.name'].fillna('', inplace=True)
    df['secondary_name'] = df['temp'] + df['secondary.name']
    df = df.drop('protein.name', axis=1)
    df = df.drop('secondary.name', axis=1)
    df = df.drop('temp', axis=1)
    df = df.iloc[:, [0, 1, 6, 7, 3, 4, 5, 2]]
    df.insert(4, 'score_threshold', None)
    df.columns = meta_col_names
    df.to_csv(output / "padloc_meta.tsv", sep="\t", header=meta_col_names)
    return df


def translate_model_padloc(data_yaml: dict, model_name: str, meta: pd.DataFrame = None, canonical: list = None):
    """
    Translate the yaml data in json dictionary

    :param data_yaml: data yaml file to translate
    """
    assert canonical is not None and isinstance(canonical, list)
    assert meta is not None and isinstance(meta, pd.DataFrame)

    def add_family(families_list: List[str], secondary_names: List[str], max_separation: int, meta: pd.DataFrame, fam_type: str):
        family_list = []
        for fam_name in families_list:
            if fam_name != 'NA':
                if fam_name in secondary_names:
                    filter_df = meta.loc[meta['secondary_name'] == fam_name]
                    prime_names = filter_df['protein_name'].dropna().unique().tolist()
                    for prime in prime_names:
                        family_list.append({'name': prime,
                                            'type': fam_type,
                                            'relation': None,
                                            'parameters': {"max_separation": max_separation}})
                family_list.append({'name': fam_name,
                                    'type': fam_type,
                                    'relation': None,
                                    'parameters': {"max_separation": max_separation}})
        return family_list

    padloc_keys = ["maximum_separation", "minimum_core", "minimum_total", "core_genes", "optional_genes", "prohibited_genes"]
    if not all(key in padloc_keys for key in data_yaml.keys()):
        raise KeyError(f"Unexpected key in PADLOC model : {model_name}."
                       f"authorized keys are : {', '.join(padloc_keys)}")

    data_json = {"name": model_name, 'func_units': [],
                 'parameters': {"max_forbidden": 0, "max_separation": 1, "min_mandatory": 1, "min_total": 1}}
    if len(canonical) > 0:
        data_json['canonical'] = canonical
    secondary_names = meta['secondary_name'].dropna().unique().tolist()
    func_unit = {'name': model_name, 'type': 'mandatory', 'parameters': {'max_forbidden': 0}}
    func_unit['parameters']['max_separation'] = data_yaml["maximum_separation"] if "maximum_separation" in data_yaml else 0
    func_unit['parameters']['min_mandatory'] = data_yaml["minimum_core"] if "minimum_core" in data_yaml else 1
    func_unit['parameters']['min_total'] = data_yaml["min_total"] if "min_total" in data_yaml else 1

    family_list = list()
    family_list += add_family(families_list=data_yaml["core_genes"] if "core_genes" in data_yaml else [],
                              secondary_names=secondary_names, meta=meta, fam_type='mandatory',
                              max_separation=func_unit['parameters']['max_separation'])
    family_list += add_family(families_list=data_yaml["optional_genes"] if "optional_genes" in data_yaml else [],
                              secondary_names=secondary_names, meta=meta, fam_type='accessory',
                              max_separation=func_unit['parameters']['max_separation'])
    family_list += add_family(families_list=data_yaml["prohibited_genes"] if "prohibited_genes" in data_yaml else [],
                              secondary_names=secondary_names, meta=meta, fam_type='forbidden',
                              max_separation=func_unit['parameters']['max_separation'])
    func_unit["families"] = family_list
    data_json['func_units'].append(func_unit)
    return data_json


def search_canonical_padloc(model_name: str, models: Path) -> List[str]:
    canonical_sys = []
    if re.search("_other", model_name):
        basename = re.split("_other", model_name)[0]
        for canon_file in models.glob("*.yaml"):
            name_canon = canon_file.stem
            if (re.search(f"{basename}_", name_canon) or re.search(f"{basename}$", name_canon)) \
                    and name_canon != model_name:
                canonical_sys.append(name_canon)
    return canonical_sys


def translate_padloc(models: Path, meta: Path, output: Path, disable_bar: bool = False):
    list_data = []
    meta_df = parse_meta_padloc(meta, output)
    for model in tqdm(list(models.rglob("*.yaml")), unit="file", disable=disable_bar):
        canonical_sys = search_canonical_padloc(model.stem, models)
        data = read_yaml(model)
        list_data.append(translate_model_padloc(data, model.stem, meta_df, canonical_sys))
    return list_data


def translate_defense_finder_model(root, model_name: str):
    """
    Translate the xml data in json dictionary

    :param root: root of xml data
    :param model_name: name of the model
    """

    def translate_gene(elem):
        try:
            name = "".join(elem.get('name').split('__')[1:])
        except KeyError:
            raise KeyError(f"In {data_json['name']} of MacSyFinder, one gene doesn't have a name")
        else:
            dict_elem = {'name': name, 'type': 'neutral', 'parameters': {"max_separation": 0}}

            for attrib, value in elem.attrib.items():
                if attrib == 'presence':
                    dict_elem['type'] = value
                elif attrib == 'inter_gene_max_space':
                    dict_elem['parameters']['max_separation'] = int(value)

            for relation in elem:
                if relation.tag == "exchangeables":
                    dict_elem["exchangeables"] = []
                    for gene in relation:
                        try:
                            name = "".join(gene.get('name').split('__')[1:])
                        except KeyError:
                            raise KeyError(f"In {data_json['name']} of Defense Finder, one gene doesn't have a name")
                        else:
                            dict_elem["exchangeables"].append(name)
                else:
                    logging.getLogger().warning("Unexpected relation")
            return dict_elem

    def translate_fu(elem):
        try:
            name = elem.get('name')
        except KeyError:
            raise KeyError(f"In {data_json['name']} of MacSyFinder, one gene doesn't have a name")
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

    data_json = {"name": model_name, 'func_units': [], 'parameters': dict()}
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
        fu_list.append({'name': data_json["name"], 'type': 'mandatory', "families": fam_list,
                        'parameters': data_json["parameters"]})
        data_json["parameters"] = {"max_forbidden": 0, "max_separation": 1, "min_mandatory": 1, "min_total": 1}
    for fu in fu_list:
        data_json["func_units"].append(fu)
    return data_json


def parse_defense_finder(model_file: Path, tmpdir: Path):
    with open(model_file, "r") as model:
        lines = model.readlines()
    lines = list(filter(lambda element: not re.search("^ ?\n", element), lines))
    newfile_path = tmpdir.joinpath(model_file.name)
    with open(newfile_path, "w") as model_test:
        index = 0
        in_fu = False
        while index < len(lines):
            line = lines[index]
            if re.search("^###############+", line):
                index += 2
            elif re.search("^#+ ", line):
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
            else:
                model_test.write(line)
            index += 1
    return newfile_path


def translate_defense_finder(models: Path, tmpdir: Path = None, disable_bar: bool = False):
    assert tmpdir is not None and isinstance(tmpdir, Path)

    with tempfile.TemporaryDirectory(prefix="Dfinder_", dir=tmpdir) as tmp:
        tmp_path = Path(tmp)
        list_data = []
        for model in tqdm(list(models.rglob("*.xml")), unit='file', disable=disable_bar):
            try:
                parse_defense_finder(model, tmp_path)
            except Exception:
                raise Exception(f"Problem to parse {model.name}")
            root = read_xml(model)
            list_data.append(translate_defense_finder_model(root, model.stem))
    return list_data


def launch_translate(models: Path, source: str, output: Path, meta_data: Path = None, tmpdir: Path = None,
                     disable_bar: bool = False):
    if source == "padloc":
        logging.getLogger().info("Begin to translate padloc models...")
        list_data = translate_padloc(models=models, meta=meta_data, output=output, disable_bar=disable_bar)
    elif source == "defense-finder":
        logging.getLogger().info("Begin to translate defense finder models...")
        list_data = translate_defense_finder(models=models, tmpdir=tmpdir, disable_bar=disable_bar)
    elif source == "macsy-finder":
        raise NotImplementedError
    else:
        raise ValueError(f"The given source: {source} is not recognize. "
                         f"Please choose between padloc, defense-finder or macsy-finder")
    logging.getLogger().info("Write models for PANORAMA...")
    for data in tqdm(list_data, unit="model", disable=disable_bar):
        write(output, data)
