#!/usr/bin/env python3
# coding:utf-8
import logging
# default libraries
import re
from typing import Dict, List

import lxml.etree
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
    df = pd.read_csv(meta, sep="\t", usecols=[0, 1, 2, 3, 4, 6, 7, 8], header=0)
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

    def add_family(families_list: List[str], secondary_names: List[str], meta: pd.DataFrame,
                   fam_type: str) -> List[Dict[str, str]]:
        family_list = []
        for fam_name in families_list:
            if fam_name == "cas_adaptation":
                # TODO try to find an elegant way to manage this exception
                filter_df = meta.loc[meta['secondary_name'] == fam_name]
                for filter_fam in filter_df['protein_name'].dropna().unique().tolist():
                    fam_dict = {'name': filter_fam,
                                'presence': fam_type}
                    family_list.append(fam_dict)
            elif fam_name != 'NA':
                fam_dict = {'name': fam_name,
                            'presence': fam_type}
                if fam_name in secondary_names:
                    filter_df = meta.loc[meta['secondary_name'] == fam_name]
                    fam_dict.update({"exchangeable": filter_df['protein_name'].dropna().unique().tolist()})
                family_list.append(fam_dict)
        return family_list

    padloc_keys = ["maximum_separation", "minimum_core", "minimum_total", "force_strand",
                   "core_genes", "secondary_genes", "neutral_genes", "prohibited_genes"]
    if not all(key in padloc_keys for key in data_yaml.keys()):
        raise KeyError(f"Unexpected key in PADLOC model : {model_name}."
                       f"authorized keys are : {', '.join(padloc_keys)}")
    if "core_genes" not in data_yaml:
        raise KeyError("Core gene must be found in padloc keys")

    data_json = {"name": model_name, 'func_units': [],
                 'parameters': {"max_forbidden": 0, "max_separation": 0, "min_mandatory": 1, "min_total": 1,
                                "max_mandatory": -1, "max_total": -1}
                 }
    if len(canonical) > 0:
        data_json['canonical'] = canonical
    secondary_names = meta['secondary_name'].dropna().unique().tolist()
    func_unit = {'name': model_name, 'presence': 'mandatory', 'same_strand': data_yaml["force_strand"],
                 'parameters':
                     {
                         "max_separation": data_yaml["maximum_separation"] if "maximum_separation" in data_yaml else 0,
                         "min_mandatory": data_yaml["minimum_core"] if "minimum_core" in data_yaml else 1,
                         "min_total": data_yaml["minimum_total"] if "minimum_total" in data_yaml else 1
                     }
                 }

    family_list = list()
    family_list += add_family(families_list=data_yaml["core_genes"],
                              secondary_names=secondary_names, meta=meta, fam_type='mandatory')
    family_list += add_family(families_list=data_yaml["secondary_genes"] if "secondary_genes" in data_yaml else [],
                              secondary_names=secondary_names, meta=meta, fam_type='accessory')
    family_list += add_family(families_list=data_yaml["prohibited_genes"] if "prohibited_genes" in data_yaml else [],
                              secondary_names=secondary_names, meta=meta, fam_type='forbidden')
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

def translate_gene(elem: lxml.etree.Element, data: dict, hmm_dict: Dict[str, List[str]], splitter: str):
    try:
        name = "".join(elem.get('name').split(splitter)[1:])
    except KeyError:
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a name")
    else:
        dict_elem = {'name': name, 'presence': 'neutral', "parameters": {}}
        exchangeable_set = set()
        if name in hmm_dict:
            exchangeable_set |= set(hmm_dict[name])
        for attrib, value in elem.attrib.items():
            if attrib == 'presence':
                dict_elem['presence'] = value
            elif attrib == 'inter_gene_max_space':
                dict_elem['parameters']['max_separation'] = int(value)
            elif attrib == 'max_nb_genes':
                dict_elem['parameters']['max_total'] = int(value)
            elif attrib == 'multi_system' and (value == 1 or value == "True"):
                dict_elem['parameters']['multi_system'] = True
            elif attrib == "loner" and (value == 1 or value == "True"):
                dict_elem["parameters"]["max_separation"] = -1
            elif attrib == 'multi_model' and (value == 1 or value == "True"):
                dict_elem['parameters']['multi_model'] = True

        for relation in elem:
            if relation.tag == "exchangeables":
                for gene in relation:
                    try:
                        exchangeable_name = splitter.join(gene.get('name').split(splitter)[1:])
                    except KeyError:
                        raise KeyError(f"In {data['name']}, one gene doesn't have a name")
                    else:
                        exchangeable_set.add(exchangeable_name)
                        if exchangeable_name in hmm_dict:
                            exchangeable_set |= set(hmm_dict[exchangeable_name])
            else:
                logging.warning("Unexpected relation")
        if len(exchangeable_set) > 0:
            dict_elem["exchangeable"] = list(exchangeable_set.difference({name}))
        return dict_elem


def translate_fu(elem, data: dict, hmm_dict: Dict[str, List[str]], splitter: str):
    try:
        name = elem.get('name')
    except KeyError:
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a name")
    else:
        dict_elem = {'name': name, 'presence': 'neutral', "families": [],
                     'parameters': {}}

        for attrib, value in elem.attrib.items():
            if attrib == 'presence':
                dict_elem['presence'] = value
            elif attrib == 'min_mandatory_genes_required':
                dict_elem['parameters']['min_mandatory'] = int(value)
            elif attrib == 'min_genes_required':
                dict_elem['parameters']['min_total'] = int(value)
            elif attrib == 'inter_gene_max_space':
                dict_elem['parameters']['max_separation'] = int(value)
            elif attrib == 'max_nb_genes':
                dict_elem['parameters']['max_total'] = -1
            elif attrib == 'multi_system' and value == 1:
                dict_elem['parameters']['multi_system'] = True
            elif attrib == "loner" and (value == 1 or value is True):
                dict_elem["parameters"]["max_separation"] = -1
            elif attrib == 'multi_model' and value == 1:
                dict_elem['parameters']['multi_model'] = True

        for family in elem:
            dict_elem["families"].append(translate_gene(family, data, hmm_dict, splitter))
        return dict_elem


def translate_macsyfinder_model(root, model_name: str, hmm_dict: Dict[str, List[str]], canonical: List[str], splitter: str):
    """
    Translate the xml data in json dictionary

    :param root: root of xml data
    :param model_name: name of the model
    """

    data_json = {"name": model_name, 'func_units': [], 'parameters': dict(), "canonical": canonical}
    if root.attrib is not None:  # Read attributes root
        for parameter, value in root.attrib.items():
            if parameter == 'inter_gene_max_space':
                data_json['parameters']['max_separation'] = int(value)
            elif parameter == 'min_mandatory_genes_required':
                data_json['parameters']['min_mandatory'] = int(value)
            elif parameter == 'min_genes_required':
                data_json['parameters']['min_total'] = int(value)
            elif parameter == 'max_nb_genes':
                data_json['parameters']['max_total'] = -1
            elif parameter == 'multi_system' and (value == 1 or value is True):
                data_json['parameters']['multi_system'] = True
            elif parameter == "loner" and (value == 1 or value is True):
                data_json["parameters"]["max_separation"] = -1
            elif parameter == 'multi_model' and (value == 1 or value is True):
                data_json['parameters']['multi_model'] = True
    data_json['parameters']['max_forbidden'] = 0

    fu_list = []
    fam_list = []
    for elem in root:
        if elem.tag == "gene":
            fam_list.append(translate_gene(elem, data_json, hmm_dict, splitter))
        elif elem.tag == "functional_unit":
            fu_list.append(translate_fu(elem, data_json, hmm_dict, splitter))
    if len(fu_list) == 0:  # only genes
        data_json["func_units"].append({'name': data_json["name"], 'presence': 'mandatory',
                                        "families": fam_list, 'parameters': data_json["parameters"]})
        data_json["parameters"] = {"max_forbidden": 0, "max_separation": 1, "min_mandatory": -1, "min_total": -1,
                                   "max_mandatory": -1, "max_total": -1}
    else:
        for fu in fu_list:
            data_json["func_units"].append(fu)
    return data_json


def parse_dfinder_hmm(hmms_path: Path):
    hmm_dict = {}
    for hmm_path in tqdm(list(hmms_path.glob('*.hmm')), unit="HMM"):
        hmm_names = []
        with open(hmm_path, 'r') as hmm_file:
            for line in hmm_file.readlines():
                if re.search('NAME', line):
                    hmm_names.append(re.split(' +', line.replace("\n", ""))[1])
        if len(hmm_names) > 1:
            hmm_dict[hmm_path.stem.split('__')[-1]] = hmm_names
        else:
            if hmm_names[0] != hmm_path.stem.split('__')[-1]:
                hmm_dict[hmm_path.stem.split('__')[-1]] = hmm_names
    return hmm_dict


def search_canonical_dfinder(model_name: str, models: Path) -> List[str]:
    canonical_sys = []
    if re.search("-Type-", model_name):
        basename = re.split("-Type-", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        for canon_file in models.glob("*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}-Subtype-{base_type}-", name_canon) and name_canon != model_name:
                canonical_sys.append(name_canon)
    elif model_name == "CAS_Cluster" or model_name == "CBASS":
        for canon_file in models.glob("*.xml"):
            if canon_file.stem != model_name:
                canonical_sys.append(canon_file.stem)
    return canonical_sys


def translate_defense_finder(models: Path, hmms_path: Path, tmpdir: Path, disable_bar: bool = False):
    assert tmpdir is not None and isinstance(tmpdir, Path)
    hmm_dict = parse_dfinder_hmm(hmms_path)
    list_data = []
    for model in tqdm(list(models.rglob("*.xml")), unit='file', disable=disable_bar):
        canonical_sys = search_canonical_dfinder(model.stem, model.parent)
        # try:
        #     parse_defense_finder(model, tmp_path)
        # except Exception:
        #     raise Exception(f"Problem to parse {model.name}")
        root = read_xml(model)
        list_data.append(translate_macsyfinder_model(root, model.stem, hmm_dict, canonical_sys, splitter="__"))
    return list_data


def search_canonical_macsyfinder(model_name: str, models: Path) -> List[str]:
    canonical_sys = []
    if re.search("-Type-", model_name):
        basename = re.split("-Type-", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        for canon_file in models.glob("*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}-Subtype-{base_type}-", name_canon) and name_canon != model_name:
                canonical_sys.append(name_canon)
    elif model_name == "CAS_Cluster" or model_name == "CBASS":
        for canon_file in models.glob("*.xml"):
            if canon_file.stem != model_name:
                canonical_sys.append(canon_file.stem)
    return canonical_sys


def parse_macsyfinder_hmm(hmms_path: Path):
    hmm_dict = {}
    for hmm_path in tqdm(list(hmms_path.glob('*.hmm')), unit="HMM"):
        hmm_names = []
        with open(hmm_path, 'r') as hmm_file:
            for line in hmm_file.readlines():
                if re.search('NAME', line):
                    hmm_names.append(re.split(' +', line.replace("\n", ""))[1])
        if len(hmm_names) > 1:
            hmm_dict[hmm_path.stem.split('_')[-1]] = hmm_names
        else:
            split_name = hmm_path.stem.split('_')
            if len(split_name) == 2:
                if hmm_names[0] != hmm_path.stem.split('_')[-1]:
                    hmm_dict[hmm_path.stem.split('_')[-1]] = hmm_names
            elif len(split_name) > 2:
                if hmm_names[0] != "_".join(split_name[1:]):
                    hmm_dict["_".join(split_name[1:])] = hmm_names
    return hmm_dict


def translate_macsyfinder(models: Path, hmms_path: Path, tmpdir: Path, disable_bar: bool = False):
    assert tmpdir is not None and isinstance(tmpdir, Path)
    hmm_dict = parse_macsyfinder_hmm(hmms_path)
    list_data = []
    for model in tqdm(list(models.rglob("*.xml")), unit='file', disable=disable_bar):
        # canonical_sys = search_canonical_dfinder(model.stem, model.parent)
        # try:
        #     parse_defense_finder(model, tmp_path)
        # except Exception:
        #     raise Exception(f"Problem to parse {model.name}")
        root = read_xml(model)
        list_data.append(translate_macsyfinder_model(root, model.stem, hmm_dict, [], splitter="_"))
    return list_data


def launch_translate(models: Path, source: str, output: Path, meta_data: Path = None, hmms_path: Path = None,
                     tmpdir: Path = None, disable_bar: bool = False):
    if source == "padloc":
        logging.info("Begin to translate padloc models...")
        list_data = translate_padloc(models=models, meta=meta_data, output=output, disable_bar=disable_bar)
    elif source == "defense-finder":
        logging.info("Begin to translate defense finder models...")
        list_data = translate_defense_finder(models=models, hmms_path=hmms_path, tmpdir=tmpdir, disable_bar=disable_bar)
    elif source == "macsy-finder":
        logging.info("Begin to translate macsy finder models...")
        list_data = translate_macsyfinder(models=models, hmms_path=hmms_path, tmpdir=tmpdir, disable_bar=disable_bar)
    else:
        raise ValueError(f"The given source: {source} is not recognize. "
                         f"Please choose between padloc, defense-finder or macsy-finder")
    logging.info("Write models for PANORAMA...")
    for data in tqdm(list_data, unit="model", disable=disable_bar):
        write(output, data)
