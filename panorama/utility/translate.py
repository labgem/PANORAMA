#!/usr/bin/env python3
# coding:utf-8

"""
Model Translation Module for PANORAMA

This module provides comprehensive functionality to translate models from different
bioinformatics databases (PADLOC, DefenseFinder, MacSyFinder-based tools) into
PANORAMA-compatible formats. It handles HMM processing, metadata parsing and
model structure conversion while maintaining compatibility across different
annotation frameworks.

Supported Sources:
- PADLOC: Prokaryotic Antiviral Defense Location predictor
- DefenseFinder: Antiviral defense systems identification
- CONJScan: Conjugation system detection
- TXSScan: Type secretion system detection
- TFFScan: Type IV-A pilus detection
"""

# default libraries
import re
from random import choice
from string import digits
from typing import Callable, Dict, List, Set, Union
import logging
from pathlib import Path
from lxml import etree as et
import lxml.etree
import yaml
import json

# installed libraries
from tqdm import tqdm
import pandas as pd
from numpy import nan
from pyhmmer.plan7 import HMM

# local libraries
from panorama.utils import mkdir
from panorama.utility.genInput import create_hmm_list_file, read_hmm, write_hmm, gen_acc


# Constants
KNOWN_SOURCES = ["padloc", "defense-finder", "CONJScan", "TXSScan", "TFFscan"]
PADLOC_MODEL_KEYS = [
    "maximum_separation", "minimum_core", "minimum_total", "force_strand",
    "core_genes", "secondary_genes", "neutral_genes", "prohibited_genes"
]
DEFAULT_EVAL_THRESHOLD = 0.1
DEFAULT_IEVAL_THRESHOLD = 0.001
PANORAMA_ACC_PREFIX = "PAN"
PANORAMA_ACC_SUFFIX_LENGTH = 6


class ModelTranslationError(Exception):
    """Custom exception for model translation errors."""
    pass


class HMMProcessingError(Exception):
    """Custom exception for HMM processing errors."""
    pass


def read_yaml(model_path: Path) -> Dict[str, Union[List[str], int, bool]]:
    """
    Read and parse a YAML file safely.

    Args:
        model_path (Path): Path to the YAML file to be read

    Returns:
        Dict: The contents of the YAML file as a Python dictionary

    Raises:
        IOError: If there is a problem opening the file
        yaml.YAMLError: If there is a problem parsing the YAML content
        ModelTranslationError: For any unexpected errors during file processing
    """
    try:
        with open(model_path, "r", encoding="utf-8") as file:
            data = yaml.safe_load(file)
            if data is None:
                logging.getLogger("PANORAMA").warning(f"Empty YAML file: {model_path}")
                return {}
            return data
    except IOError as e:
        raise IOError(f"Problem opening {model_path}: {e}") from e
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Problem parsing YAML file {model_path}: {e}") from e
    except Exception as e:
        raise ModelTranslationError(f"Unexpected error reading {model_path}: {e}") from e


def read_xml(model_path: Path) -> et.Element:
    """
    Read and parse an XML file with security considerations.

    Args:
        model_path (Path): Path to the XML file to be read

    Returns:
        et.Element: The root element of the parsed XML document

    Raises:
        IOError: If there is a problem opening the file
        et.XMLSyntaxError: If there is a problem parsing the XML content
        ModelTranslationError: For any unexpected errors during file processing
    """
    try:
        # Configure parser for security: disable external entities and comments
        parser = et.XMLParser(
            remove_comments=True,
            resolve_entities=False,
            no_network=True,
            huge_tree=False
        )
        tree = et.parse(str(model_path.absolute()), parser=parser)
        return tree.getroot()
    except IOError as e:
        raise IOError(f"Problem opening {model_path}: {e}") from e
    except et.XMLSyntaxError as e:
        raise et.XMLSyntaxError(f"Problem parsing XML file {model_path}: {e}") from e
    except Exception as e:
        raise ModelTranslationError(f"Unexpected error reading {model_path}: {e}") from e


def write_model(output_path: Path, model_data: Dict[str, Union[str, List[Dict], Dict[str, int], List[str]]]) -> Path:
    """
    Write a translated model to a JSON file with proper formatting.

    Args:
        output_path (Path): Path to the output directory
        model_data (Dict): Dictionary containing the model data to write

    Returns:
        Path: Path to the written JSON file

    Raises:
        IOError: If there is a problem writing the file
        KeyError: If the model_data doesn't contain the required 'name' field
    """
    if 'name' not in model_data:
        raise KeyError("Model data must contain 'name' field")

    output_file = output_path / f"{model_data['name']}.json"

    try:
        with open(output_file, "w", encoding="utf-8") as file:
            json.dump(model_data, file, indent=2, ensure_ascii=False)
        return output_file
    except IOError as e:
        raise IOError(f"Problem writing model {model_data['name']} to {output_file}: {e}") from e


def parse_meta_padloc(meta_path: Path) -> pd.DataFrame:
    """
    Parse the PADLOC metadata file and return a structured DataFrame.

    The function processes the PADLOC hmm_meta.txt file which contains HMM metadata
    including accession numbers, names, thresholds and descriptions. It handles
    protein name parsing and secondary name merging.

    Args:
        meta_path (Path): Path to the PADLOC metadata file (hmm_meta.txt)

    Returns:
        pd.DataFrame: DataFrame with parsed metadata indexed by accession number,
                     containing columns: name, protein_name, secondary_name,
                     score_threshold, eval_threshold, ieval_threshold,
                     hmm_cov_threshold, target_cov_threshold, description

    Raises:
        IOError: If the metadata file cannot be read
        ValueError: If the metadata format is unexpected
    """
    def _merge_secondary_names(row: pd.Series) -> str:
        """Helper function to merge secondary name columns."""
        temp_val = row.get('temp', '')
        secondary_val = row.get('secondary.name', '')

        if pd.isna(temp_val) and pd.isna(secondary_val):
            return ''
        elif pd.isna(temp_val):
            return str(secondary_val)
        elif pd.isna(secondary_val):
            return str(temp_val)
        else:
            return f"{temp_val},{secondary_val}"

    # Define expected column structure
    meta_columns = [
        "accession", "name", "protein_name", "secondary_name", "score_threshold",
        "eval_threshold", "ieval_threshold", "hmm_cov_threshold",
        "target_cov_threshold", "description"
    ]

    # Read the TSV file with specific columns
    df = pd.read_csv(meta_path, sep="\t", usecols=[0, 1, 2, 3, 4, 6, 7, 8], header=0)

    # Parse protein names (format: main_name|secondary_name)
    df[['protein_name', 'temp']] = df['protein.name'].str.split('|', expand=True)

    # Merge secondary names
    df['secondary_name'] = df.apply(_merge_secondary_names, axis=1)

    # Clean up temporary columns
    df = df.drop(columns=['protein.name', 'secondary.name', 'temp'], errors='ignore')

    # Reorder columns and add missing ones
    df = df.iloc[:, [0, 1, 6, 7, 3, 4, 5, 2]]  # Reorder existing columns

    # Insert missing threshold columns
    df.insert(4, 'score_threshold', None)
    df.insert(5, 'eval_threshold', None)

    # Set proper column names
    df.columns = meta_columns

    # Set index and handle missing descriptions
    df = df.set_index('accession')
    df['description'] = df['description'].fillna('unknown')

    return df


def _add_families_to_functional_unit(
        families_list: List[str],
        secondary_names: List[str],
        family_type: str,
        metadata_df: pd.DataFrame,
        seen_families: Set[str]
) -> List[Dict[str, Union[str, List[str]]]]:
    """
    Add families to a functional unit with proper metadata integration.

    This helper function processes a list of family names and creates properly
    formatted family dictionaries for PANORAMA models. It handles special cases
    like cas_adaptation and manages exchangeable proteins.

    Args:
        families_list (List[str]): List of family names to process
        secondary_names (List[str]): List of known secondary names for exchangeable lookup
        family_type (str): Type of family ('mandatory', 'accessory', 'forbidden', 'neutral')
        metadata_df (pd.Dataframe): DataFrame containing family metadata
        seen_families (Set[str]): Set to track already processed families (modified in place)

    Returns:
        List of family dictionaries with structure:
        - name: str (family name)
        - presence: str (family type)
        - exchangeable: List[str] (optional, list of exchangeable proteins)
    """
    family_list = []

    for family_name in families_list:
        if family_name in seen_families or family_name in ['NA', 'cas_accessory']:
            continue

        # Special handling for cas_adaptation families
        if family_name == "cas_adaptation":
            filtered_df = metadata_df.loc[metadata_df['secondary_name'] == family_name]
            for protein_name in filtered_df['protein_name'].dropna().unique():
                if protein_name not in seen_families:
                    family_dict = {
                        'name': protein_name,
                        'presence': family_type
                    }
                    seen_families.add(protein_name)

                    # Add exchangeable proteins if available
                    if protein_name in secondary_names:
                        exchangeable_proteins = filtered_df['protein_name'].dropna().unique().tolist()
                        if len(exchangeable_proteins) > 1:
                            family_dict["exchangeable"] = exchangeable_proteins

                    family_list.append(family_dict)
        else:
            # Standard family processing
            family_dict = {
                'name': family_name,
                'presence': family_type
            }
            seen_families.add(family_name)

            # Add exchangeable proteins if this family has secondary names
            if family_name in secondary_names:
                filtered_df = metadata_df.loc[metadata_df['secondary_name'] == family_name]
                exchangeable_proteins = filtered_df['protein_name'].dropna().unique().tolist()
                if exchangeable_proteins:
                    family_dict["exchangeable"] = exchangeable_proteins

            family_list.append(family_dict)

    return family_list


def translate_model_padloc(
        data_yaml: Dict[str, Union[List[str], int, bool]],
        model_name: str,
        metadata_df: pd.DataFrame,
        canonical_models: List[str] = None
) -> Dict[str, Union[str, List[Dict], Dict[str, int], List[str]]]:
    """
    Translate a PADLOC model from YAML format to PANORAMA JSON format.

    This function converts PADLOC defense system models into the standardized
    PANORAMA format, handling gene categories, parameters and canonical relationships.

    Args:
        data_yaml: PADLOC model data loaded from YAML file
        model_name: Name identifier for the model
        metadata_df: DataFrame containing HMM metadata for gene information
        canonical_models: List of canonical model names (optional)

    Returns:
        Dict: Translated model in PANORAMA format with structure:
        - name: str (model name)
        - func_units: List[Dict] (functional units with families)
        - parameters: Dict (global model parameters)
        - canonical: List[str] (optional, canonical model references)

    Raises:
        KeyError: If required keys are missing from the PADLOC model
        AssertionError: If input parameters are invalid
        ModelTranslationError: For translation-specific errors
    """
    assert canonical_models is not None and isinstance(canonical_models, list), "canonical_models must be a list"
    assert isinstance(metadata_df, pd.DataFrame), "metadata_df must be a pandas DataFrame"

    # Validate inputs
    if canonical_models is None:
        canonical_models = []

    # Validate PADLOC model structure
    unexpected_keys = set(data_yaml.keys()) - set(PADLOC_MODEL_KEYS)
    if unexpected_keys:
        raise KeyError(
            f"Unexpected keys in PADLOC model '{model_name}': {unexpected_keys}. "
            f"Expected keys: {PADLOC_MODEL_KEYS}"
        )

    if "core_genes" not in data_yaml:
        raise KeyError(f"Missing required 'core_genes' key in PADLOC model '{model_name}'")

    # Initialize PANORAMA model structure
    panorama_model = {
        "name": model_name,
        'func_units': [],
        'parameters': {
            "transitivity": 0,
            "window": 1,
            "min_mandatory": 1,
            "min_total": 1
        }
    }

    if canonical_models:
        panorama_model['canonical'] = canonical_models

    # Extract secondary names for exchangeable protein lookup
    secondary_names = metadata_df['secondary_name'].dropna().unique().tolist()

    # Create a functional unit
    functional_unit = {
        'name': model_name,
        'presence': 'mandatory',
        'same_strand': data_yaml.get("force_strand", False),
        'parameters': {
            "transitivity": data_yaml.get("maximum_separation", 0),
            "min_mandatory": data_yaml.get("minimum_core", 1),
            "min_total": data_yaml.get("minimum_total", 1)
        }
    }

    # Process gene families by category
    family_list = []
    seen_families = set()

    # Core genes (mandatory)
    family_list.extend(_add_families_to_functional_unit(
        data_yaml["core_genes"], secondary_names, 'mandatory', metadata_df, seen_families
    ))

    # Secondary genes (accessory)
    if "secondary_genes" in data_yaml:
        family_list.extend(_add_families_to_functional_unit(
            data_yaml["secondary_genes"], secondary_names, 'accessory', metadata_df, seen_families
        ))

    # Prohibited genes (forbidden)
    if "prohibited_genes" in data_yaml:
        family_list.extend(_add_families_to_functional_unit(
            data_yaml["prohibited_genes"], secondary_names, 'forbidden', metadata_df, seen_families
        ))

    # Neutral genes
    if "neutral_genes" in data_yaml:
        family_list.extend(_add_families_to_functional_unit(
            data_yaml["neutral_genes"], secondary_names, 'neutral', metadata_df, seen_families
        ))

    functional_unit["families"] = family_list
    panorama_model['func_units'].append(functional_unit)

    return panorama_model


def search_canonical_padloc(model_name: str, models_dir: Path) -> List[str]:
    """
    Search for canonical models related to a PADLOC model.

    PADLOC uses a naming convention where models ending with '_other' are
    variants of base models. This function identifies the canonical (base)
    models for such variants.

    Args:
        model_name: Name of the current model being processed
        models_dir: Directory containing all PADLOC model files

    Returns:
        List[str]: Names of canonical models related to the input model
    """
    # Extract base name (everything before '_other')
    base_name = re.split("_other", model_name)[0]

    # Create regex pattern to match canonical models
    pattern = re.compile(f"^{re.escape(base_name)}(_.*)?$")

    # Filter files using both glob and regex
    return [
        model_file.stem
        for model_file in models_dir.rglob(f"{base_name}*.yaml")
        if model_file.stem != model_name and pattern.match(model_file.stem)
    ]


def translate_padloc(
        padloc_db: Path,
        output: Path,
        binary_hmm: bool = False,
        hmm_coverage: float = None,
        target_coverage: float = None,
        force: bool = False,
        disable_bar: bool = False
) -> List[Dict]:
    """
    Translate all PADLOC models to PANORAMA format.

    This function orchestrates the complete translation process for the PADLOC database:
    1. Parses HMM metadata
    2. Creates HMM list file for annotation
    3. Translates all model files from YAML to JSON format
    4. Handles canonical model relationships

    Args:
        padloc_db: Path to the PADLOC database directory
        output: Path to output directory for translated files
        binary_hmm: Whether to output HMMs in binary format
        hmm_coverage: Global HMM coverage threshold (optional)
        target_coverage: Global target coverage threshold (optional)
        force: Whether to overwrite existing output files
        disable_bar: Whether to disable progress bars

    Returns:
        List[Dict]: List of translated PANORAMA models

    Raises:
        FileNotFoundError: If required PADLOC database files are missing
        ModelTranslationError: If translation fails for any model
    """
    # Validate database structure
    required_files = ["hmm_meta.txt", "sys", "hmm"]
    for required_file in required_files:
        if not (padloc_db / required_file).exists():
            raise FileNotFoundError(f"Required PADLOC file/directory missing: {padloc_db / required_file}")

    # Parse metadata
    logging.getLogger("PANORAMA").info("Parsing PADLOC metadata...")
    metadata_df = parse_meta_padloc(padloc_db / "hmm_meta.txt")

    # Create the HMM list file
    logging.getLogger("PANORAMA").info("Creating HMM list file...")
    create_hmm_list_file(
        hmm_path=[padloc_db / "hmm"],
        output=output,
        metadata_df=metadata_df,
        hmm_coverage=hmm_coverage,
        target_coverage=target_coverage,
        binary_hmm=binary_hmm,
        force=force,
        disable_bar=disable_bar
    )

    # Translate models
    logging.getLogger("PANORAMA").info("Translating PADLOC models...")
    translated_models = []
    models_dir = padloc_db / "sys"

    model_files = list(models_dir.rglob("*.yaml"))

    for model_file in tqdm(model_files, unit="file", desc="Translating models", disable=disable_bar):
        try:
            # Find canonical models
            canonical_models = []
            if re.search("_other", model_file.stem):  # Only process '_other' variant models
                canonical_models = search_canonical_padloc(model_file.stem, models_dir)

            # Load and translate model
            model_data = read_yaml(model_file)
            translated_model = translate_model_padloc(
                model_data, model_file.stem, metadata_df, canonical_models
            )
            translated_models.append(translated_model)

        except Exception as e:
            logging.getLogger("PANORAMA").error(f"Failed to translate PADLOC model {model_file.stem}: {e}")
            raise ModelTranslationError(f"Translation failed for {model_file.stem}") from e

    logging.getLogger("PANORAMA").info(f"Successfully translated {len(translated_models)} PADLOC models")
    return translated_models


def translate_gene(elem: lxml.etree.Element, data: dict, hmm_df: pd.DataFrame, transitivity_mut: Callable) -> dict:
    """ Translate a gene from Defense-finder or MacSyFinder models into PANORAMA family models

    Args:
        elem: XML element corresponding to the gene to be translated
        data: Dictionary containing PANORAMA model information
        hmm_df: HMM information useful to translate gene into PANORAMA family models

    Returns:
        dict: Dictionary containing PANORAMA family model information
    """
    try:
        hmm_info = hmm_df.loc[elem.get('name')]
    except KeyError:
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a name")
    else:
        if isinstance(hmm_info, pd.Series):
            dict_elem = {'name': hmm_info["protein_name"], 'presence': 'neutral', "parameters": {}}
            exchangeable_set = {hmm_info.name, hmm_info['secondary_name']}
        else:  # pd.Dataframe
            prot_name_list = hmm_info["protein_name"].mode().to_list()
            secondary_name_list = hmm_info["secondary_name"].mode().to_list()
            dict_elem = {'name': prot_name_list[0], 'presence': 'neutral', "parameters": {}}
            exchangeable_set = set(prot_name_list) | {elem.get('name')} | set(secondary_name_list)
        for attrib, value in elem.attrib.items():
            if attrib == 'presence':
                dict_elem['presence'] = value
            elif attrib == 'inter_gene_max_space':
                dict_elem['parameters']['transitivity'] = transitivity_mut(int(value))
            elif attrib == 'multi_system' and (value == 1 or value == "True"):
                dict_elem['parameters']['multi_system'] = True
            elif attrib == "loner" and (value == 1 or value == "True"):
                dict_elem["parameters"]["transitivity"] = -1
            elif attrib == 'multi_model' and (value == 1 or value == "True"):
                dict_elem['parameters']['multi_model'] = True

        for relation in elem:
            if relation.tag == "exchangeables":
                for gene in relation:
                    try:
                        exchangeable_info = hmm_df.loc[gene.get('name')]
                    except KeyError:
                        raise KeyError(f"In {data['name']}, one gene doesn't have a name")
                    else:
                        if isinstance(exchangeable_info, pd.Series):
                            exchangeable_set |= {exchangeable_info["protein_name"], exchangeable_info.name,
                                                 exchangeable_info['secondary_name']}
                        else:  # pd.Dataframe
                            prot_name_set = set(exchangeable_info["protein_name"].to_list())
                            secondary_name_set = set(exchangeable_info["secondary_name"].to_list())
                            name_set = set(exchangeable_info.index.to_list())
                            exchangeable_set |= prot_name_set | name_set | secondary_name_set
            else:
                logging.getLogger("PANORAMA").warning("Unexpected relation")
        if len(exchangeable_set) > 0:
            if isinstance(hmm_info, pd.Series):
                dict_elem["exchangeable"] = list(exchangeable_set.difference({hmm_info["protein_name"], ""}))
            else:  # pd.Dataframe
                prot_name_list = hmm_info["protein_name"].mode().to_list()
                dict_elem["exchangeable"] = list(exchangeable_set.difference({prot_name_list[0], ""}))
        return dict_elem


def translate_fu(elem: lxml.etree.Element, data: dict, hmm_df: pd.DataFrame, transitivity_mut: Callable):
    """ Translate a gene from Defense-finder or MacSyFinder models into PANORAMA functional unit models

    Args:
        elem: XML element corresponding to the gene to be translated
        data: Dictionary containing PANORAMA model information
        hmm_df: HMM information useful to translate gene into PANORAMA family models

    Returns:
        dict: Dictionary containing PANORAMA functional unit model information
    """
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
                dict_elem['parameters']['transitivity'] = transitivity_mut(int(value))
            elif attrib == 'multi_system' and value == 1:
                dict_elem['parameters']['multi_system'] = True
            elif attrib == "loner" and (value == 1 or value is True):
                dict_elem["parameters"]["transitivity"] = -1
            elif attrib == 'multi_model' and value == 1:
                dict_elem['parameters']['multi_model'] = True

        for family in elem:
            dict_elem["families"].append(translate_gene(family, data, hmm_df, transitivity_mut))
        return dict_elem


def translate_macsyfinder_model(root: et.Element, model_name: str, hmm_df: pd.DataFrame,
                                canonical: List[str], transitivity_mut: Callable) -> dict:
    """
    Translate macsyfinder model into PANORAMA model

    Args:
        root: Model root information
        model_name: name of the model
        hmm_df: HMM information to translate models into PANORAMA models
        canonical: List of canonical models

    Returns:
        dict: PANORAMA model information
    """
    data_json = {"name": model_name, "func_units": [], "parameters": {}}
    if len(canonical) > 0:
        data_json["canonical"] = canonical
    if root.attrib is not None:  # Read attributes root
        for parameter, value in root.attrib.items():
            if parameter == 'inter_gene_max_space':
                data_json['parameters']['transitivity'] = transitivity_mut(int(value))
            elif parameter == 'min_mandatory_genes_required':
                data_json['parameters']['min_mandatory'] = int(value)
            elif parameter == 'min_genes_required':
                data_json['parameters']['min_total'] = int(value)
            elif parameter == 'multi_system' and (value == 1 or value is True):
                data_json['parameters']['multi_system'] = True
            elif parameter == "loner" and (value == 1 or value is True):
                data_json["parameters"]["transitivity"] = -1
            elif parameter == 'multi_model' and (value == 1 or value is True):
                data_json['parameters']['multi_model'] = True

    fu_list = []
    fam_list = []
    for elem in root:
        if elem.tag == "gene":
            fam_list.append(translate_gene(elem, data_json, hmm_df, transitivity_mut))
        elif elem.tag == "functional_unit":
            fu_list.append(translate_fu(elem, data_json, hmm_df, transitivity_mut))
    if len(fu_list) == 0:  # only genes
        data_json["func_units"].append({'name': data_json["name"], 'presence': 'mandatory',
                                        "families": fam_list, 'parameters': data_json["parameters"]})
        data_json["parameters"] = {"transitivity": 0, "window": 1,  "min_mandatory": 1, "min_total": 1}
    else:
        for fu in fu_list:
            data_json["func_units"].append(fu)
    return data_json


def parse_macsyfinder_hmm(hmm: HMM, hmm_file: Path, panorama_acc: Set[str]) -> Dict[str, Union[str, int, float]]:
    """
    Read DefenseFinder HMM and get information for PANORAMA annotation step

    Args:
        hmm: HMM object
        hmm_file: Path to the DefenseFinder HMM
        panorama_acc: Set of PANORAMA accession ID for HMM

    Returns:
        List[Dict[str, Union[str, int, float]]]: List of dictionaries with all necessary information to write HMM metadata
    """
    hmm_dict = {"name": "", 'accession': "", "length": len(hmm.consensus), 'protein_name': "", 'secondary_name': "",
                "score_threshold": nan, "eval_threshold": 0.1, "ieval_threshold": 0.001, "hmm_cov_threshold": nan,
                "target_cov_threshold": nan, "description": ""}

    name = hmm.name.decode('UTF-8').split()
    if len(name) > 1:
        logging.getLogger("PANORAMA").error(f"HMM with problem to translate: {hmm_file.absolute().as_posix()}")
        raise IOError("It's not possible to get the name of the HMM. Please report an issue on GitHub.")
    else:
        if name[0] != hmm_file.stem.split('__')[-1]:
            hmm_dict["name"] = hmm_file.stem
            hmm_dict["protein_name"] = hmm_file.stem.split('__')[-1]
        else:
            hmm_dict["name"] = hmm_file.stem
            hmm_dict["protein_name"] = name[0]
    hmm.name = hmm_dict["protein_name"].encode('UTF-8')

    if hmm.accession is None or hmm.accession == "".encode("UTF-8"):
        hmm_dict["accession"] = gen_acc("PAN" + ''.join(choice(digits) for _ in range(6)), panorama_acc)
        hmm.accession = hmm_dict["accession"].encode("UTF-8")
    else:
        hmm_dict["accession"] = hmm.accession.decode("UTF-8")

    if hmm.description is not None:
        hmm_dict["description"] = hmm.description.decode("UTF-8")

    return hmm_dict


def create_macsyfinder_hmm_list(hmms_path: Path, output: Path, binary_hmm: bool = False, hmm_coverage: float = None,
                                target_coverage: float = None, force: bool = False,
                                disable_bar: bool = False) -> pd.DataFrame:
    """
    Read and parse all DefenseFinder HMM files and write a HMM list file for PANORAMA annotation step

    Args:
        hmms_path: Path to the HMM directory
        output: Path to the output directory where HMM list file will be written
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        force: Flag to overwrite the output directory
        disable_bar: Flag to disable the progress bar

    Returns:
        pd.Dataframe: Dataframe containing the HMM information needed for model translation
    """
    logging.getLogger("PANORAMA").info("Begin to translate HMM files from DefenseFinder...")
    hmm_dir = mkdir(output / 'hmm', force, erase=True)
    hmm_info_list = []
    panorama_acc = set()
    for hmm_file in tqdm(list(hmms_path.glob('*.hmm')), unit="HMM", desc='Translate HMM', disable=disable_bar):
        hmm_list = read_hmm(hmm_path=hmm_file)
        for hmm in hmm_list:
            hmm_dict = parse_macsyfinder_hmm(hmm, hmm_file, panorama_acc)
            hmm_dict["path"] = write_hmm(hmm, hmm_dir, binary_hmm).absolute().as_posix()
            hmm_info_list.append(hmm_dict)

    hmm_df = pd.DataFrame(hmm_info_list)
    hmm_df = hmm_df.sort_values(by=["name", "accession", "protein_name"], ascending=[True, True, True])
    hmm_df = hmm_df[['name', 'accession', 'path', 'length', 'protein_name', 'secondary_name', 'score_threshold',
                     'eval_threshold', "ieval_threshold", 'hmm_cov_threshold', 'target_cov_threshold', 'description']]
    if hmm_coverage is not None:
        logging.getLogger("PANORAMA").warning("HMM coverage threshold will be overwritten")
        hmm_df["hmm_cov_threshold"] = hmm_coverage
    if target_coverage is not None:
        logging.getLogger("PANORAMA").warning("target coverage threshold will be overwritten")
        hmm_df["target_cov_threshold"] = target_coverage
    hmm_df.to_csv(output / "hmm_list.tsv", sep="\t", index=False)
    hmm_df.set_index('name', inplace=True)
    logging.getLogger("PANORAMA").info("HMM list file created.")
    return hmm_df


def search_canonical_macsyfinder(model_name: str, models: Path) -> List[str]:
    """
    Search Canonical models for MacSyFinder models

    Args:
        model_name: Name of the model
        models: Path to all other models

    Returns:
        List[str]: List with name of canonical model
    """
    canonical_sys = []
    if re.search("-Type-", model_name):
        basename = re.split("-Type-", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        for canon_file in models.rglob(f"{base_class}*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}-Subtype-{base_type}-", name_canon) and name_canon != model_name:
                logging.getLogger("PANORAMA").debug("")
                canonical_sys.append(name_canon)
    elif re.search("_Type_", model_name):
        basename = re.split("_Type_", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        if len(base_type.split('_')) == 2:
            base_type, subtype = base_type.split('_')
        else:
            subtype = ""
        for canon_file in models.rglob("*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}_Type_{base_type}_{subtype}", name_canon) and name_canon != model_name:
                logging.getLogger("PANORAMA").debug("")
                canonical_sys.append(name_canon)
    elif model_name in ["CAS_Cluster", "CBASS", "Wadjet", "BREX"]:
        for canon_file in models.rglob(f"{model_name}*.xml"):
            if canon_file.stem != model_name:
                canonical_sys.append(canon_file.stem)
    return canonical_sys


def get_models_path(models: Path, source: str) -> Dict[str, Path]:
    """
    Associate a unique name for the model to its path in function of the source

    Args:
        models: Models database path
        source: Name of the source.

    Returns:
        Dict[str, Path]: A dictionary with model name for PANORAMA as key and the model path as value
    """
    model2path = {}

    for model in models.rglob("*.xml"):
        if source == "defense-finder":
            model2path[model.stem] = model
        elif source == "CONJScan":
            model2path[f"{model.parent.stem}_{model.stem}"] = model
        elif source == "TXSScan":
            model2path[model.stem] = model
        elif source == "TFFscan":
            model2path[model.stem] = model
        else:
            raise ValueError(f"Unknown source: {source}")

    return model2path


def translate_macsyfinder(macsy_db: Path, output: Path, binary_hmm: bool = False, hmm_coverage: float = None,
                          target_coverage: float = None, source: str = "", force: bool = False,
                          disable_bar: bool = False) -> List[dict]:
    """
    Translate MacSyFinder models into PANORAMA models and write all necessary file for PANORAMA steps

    Args:
        macsy_db: Models database path
        output: Path to output directory for PANORAMA files
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        source: Name of the source. Defaults to ""
        force: Flag to force overwrite files
        disable_bar: Flag to disable progress bar

    Returns:
         List[dict]: List of dictionaries containing translated models
    """
    def default_transitivity_mut(x: int):
        return x

    def CONJScan_transitivity_mut(x: int):
        return int(x/2) if int(x/2) >= x/2 else int(x/2) + 1

    if hmm_coverage is None:
        if source == "defense-finder":
            hmm_coverage = 0.4
        else:
            hmm_coverage = 0.5

    hmm_df = create_macsyfinder_hmm_list(macsy_db / "profiles", output, binary_hmm, hmm_coverage, target_coverage,
                                         force, disable_bar)
    list_data = []
    logging.getLogger('PANORAMA').info(f"Begin to translate {source} models")
    model_path = macsy_db / "definitions"
    if source == "CONJScan":
        transitivity_mut = CONJScan_transitivity_mut
    else:
        transitivity_mut = default_transitivity_mut
    for name, path in tqdm(get_models_path(model_path, source).items(), unit='file',
                           desc='Translate models', disable=disable_bar):
        try:
            canonical_sys = search_canonical_macsyfinder(name, model_path)
            root = read_xml(path)
            list_data.append(translate_macsyfinder_model(root, name, hmm_df, canonical_sys, transitivity_mut))
        except Exception as error:
            logging.getLogger('PANORAMA').error(f"Error translating {path.as_posix()}")
            raise Exception(f"Error translating {name} due to error: {error}")
    return list_data


def launch_translate(db: Path, source: str, output: Path,  binary_hmm: bool = False, hmm_coverage: float = None,
                     target_coverage: float = None, force: bool = False, disable_bar: bool = False):
    """
    Launch models translation process and write results for PANORAMA

    Args:
        db: Path to the models database that need to be translated
        source: Name of the source model. PADLOC, DefenseFinder and MacSyFinder are supported for now
        output: Path to the output directory to write all files needed for PANORAMA
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
        force: Flag to force overwrite existing files
        disable_bar: Flag to disable progress bar
    """

    if source == "padloc":
        list_data = translate_padloc(padloc_db=db, output=output, binary_hmm=binary_hmm, hmm_coverage=hmm_coverage,
                                     target_coverage=target_coverage, force=force, disable_bar=disable_bar)
    elif source in ["defense-finder", "CONJScan", "TXSScan", "TFFscan"]:
        list_data = translate_macsyfinder(macsy_db=db, output=output, binary_hmm=binary_hmm, hmm_coverage=hmm_coverage,
                                          target_coverage=target_coverage, source=source,
                                          force=force, disable_bar=disable_bar)
    else:
        raise ValueError(f"The given source: {source} is not recognize. "
                         f"Please choose between padloc, defense-finder or macsy-finder")
    logging.getLogger("PANORAMA").info("Write models for PANORAMA...")
    model_list = []
    model_dir = mkdir(output / 'models', force)
    for data in tqdm(list_data, unit="model", desc='Write models', disable=disable_bar):
        model_path = write_model(model_dir, data)
        model_list.append([data['name'], model_path.resolve()])
    model_df = pd.DataFrame(model_list, columns=['name', 'path'])
    model_df = model_df.sort_values('name')
    if source == "CONJScan":
        plasmid_mask = model_df['name'].str.contains('Plasmids_')
        chromosome_mask = model_df['name'].str.contains('Chromosome_')
        model_df[plasmid_mask].to_csv(output / 'models_plasmids_list.tsv', sep="\t", header=False, index=False)
        model_df[chromosome_mask].to_csv(output / 'models_chromosome_list.tsv', sep="\t", header=False, index=False)
    model_df.to_csv(output / "models_list.tsv", sep="\t", header=False, index=False)
