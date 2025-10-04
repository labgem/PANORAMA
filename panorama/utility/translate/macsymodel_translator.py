#!/usr/bin/env python3
# coding:utf-8

"""
Model Translation Module for PANORAMA

This module provides comprehensive functionality to translate models from MacSyFinder-based tools into
PANORAMA-compatible formats. It handles HMM processing, metadata parsing and
model structure conversion.

Supported Sources:
- DefenseFinder: Antiviral defense systems identification
- CasFinder: CRISPR-Cas Antiviral defense systems identification
- CONJScan: Conjugation system detection
- TXSScan: Type secretion system detection
- TFFScan: Type IV-A pilus detection
"""

# default libraries
import logging
import re
from pathlib import Path
from typing import Callable, Dict, List, Set, Union

# installed libraries
import lxml.etree
import pandas as pd
from lxml import etree as et
from numpy import nan
from pyhmmer.plan7 import HMM
from tqdm import tqdm

# local libraries
from panorama.utility.genInput import (
    process_hmm_accession,
    process_hmm_name,
    read_hmm,
    write_hmm,
)
from panorama.utility.translate import ModelTranslationError, read_xml
from panorama.utils import is_true_value, mkdir

DEFAULT_EVAL_THRESHOLD = 0.1
DEFAULT_IEVAL_THRESHOLD = 0.001
# Source-specific coverage thresholds
COVERAGE_THRESHOLDS = {"defense-finder": 0.4, "default": 0.5}


def process_attributes(
    elem: lxml.etree.Element,
    dict_elem: Dict[str, Union[str, Dict[str, int]]],
    transitivity_mut: Callable[[int], int],
) -> Dict[str, Union[str, Dict[str, int], List]]:
    """
    Process XML attributes for an element (Family, Functional Unit or Model).

    Args:
        elem (lxml.etree.Element): XML element with attributes to process
        dict_elem (Dict): Dictionary to update
        transitivity_mut (Callable): Function to transform transitivity values

    Returns:
        Dict: Updated gene dictionary
    """
    for attrib, value in elem.attrib.items():
        if attrib == "presence":
            dict_elem["presence"] = value
        elif attrib == "min_mandatory_genes_required":
            dict_elem["parameters"]["min_mandatory"] = int(value)
        elif attrib == "min_genes_required":
            dict_elem["parameters"]["min_total"] = int(value)
        elif attrib == "inter_gene_max_space":
            dict_elem["parameters"]["transitivity"] = transitivity_mut(int(value))
        elif attrib == "multi_system" and is_true_value(value):
            dict_elem["parameters"]["multi_system"] = True
        elif attrib == "loner" and is_true_value(value):
            dict_elem["parameters"]["transitivity"] = -1  # Special case for loner genes
        elif attrib == "multi_model" and is_true_value(value):
            dict_elem["parameters"]["multi_model"] = True
    return dict_elem


def process_exchangeable(
    elem: lxml.etree.Element,
    data: Dict[str, Union[str, List, Dict[str, int]]],
    hmm_df: pd.DataFrame,
    exchangeable_set: Set[str],
) -> Set[str]:
    """
    Process exchangeable genes from XML relations.

    Args:
        elem (lxml.etree.Element): XML element containing relation
        data (Dict): Model data dictionary
        hmm_df (pd.DataFrame): HMM information DataFrame
        exchangeable_set (Set): Set of exchangeable gene names to update

    Returns:
        Set: The updated set of exchangeable gene names
    """
    for relation in elem:
        if relation.tag == "exchangeables":
            for gene in relation:
                gene_name = gene.get("name")
                if not gene_name:
                    raise KeyError(f"In {data['name']}, one exchangeable gene doesn't have a name")

                try:
                    exchangeable_info = hmm_df.loc[gene_name]
                except KeyError as error:
                    raise KeyError(f"Exchangeable gene '{gene_name}' not found in HMM dataframe") from error

                # Add exchangeable gene information to set
                if isinstance(exchangeable_info, pd.Series):
                    exchangeable_set.update(
                        {
                            exchangeable_info["protein_name"],
                            exchangeable_info.name,
                            exchangeable_info["secondary_name"],
                        }
                    )
                else:  # DataFrame case
                    prot_name_set = set(exchangeable_info["protein_name"].to_list())
                    secondary_name_set = set(exchangeable_info["secondary_name"].to_list())
                    name_set = set(exchangeable_info.index.to_list())
                    exchangeable_set.update(prot_name_set | name_set | secondary_name_set)
        else:
            logging.getLogger("PANORAMA").warning(f"Unexpected relation tag: {relation.tag}")

    return exchangeable_set


def translate_gene(
    elem: lxml.etree.Element,
    data: Dict[str, Union[str, List, Dict[str, int]]],
    hmm_df: pd.DataFrame,
    transitivity_mut: Callable[[int], int],
) -> Dict[str, Union[str, Dict[str, int]]]:
    """
    Translate a gene element from MacSyFinder models (or like) into PANORAMA family models.

    This function processes XML gene elements and converts them to PANORAMA Family,
    handling various attributes like presence, transitivity, and exchangeable genes.

    Args:
        elem (lxml.etree.Element): XML element corresponding to the gene to be translated
        data (Dict[str, Any]): Dictionary containing PANORAMA model information
        hmm_df (pd.DataFrame): HMM information DataFrame indexed by gene name for translation
        transitivity_mut (Callable[[int], int]): Function to mutate transitivity values

    Returns:
        Dict[str, Any]: Dictionary containing PANORAMA family model information with keys:
            - name: protein name
            - presence: gene presence requirement ('neutral', 'mandatory', etc.)
            - parameters: dict with transitivity, multi_system, multi_model flags
            - exchangeable: list of exchangeable gene names (if any)

    Raises:
        KeyError: If the gene name is missing or not found in HMM dataframe
        ModelTranslationError: For unexpected translation errors
    """
    gene_name = elem.get("name")
    if not gene_name:
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a name")

    try:
        hmm_info = hmm_df.loc[gene_name]
    except KeyError as error:
        raise KeyError(f"Gene '{gene_name}' not found in HMM dataframe for model {data['name']}") from error

    # Initialize a gene dictionary with default values
    dict_elem = {
        "name": "",
        "presence": "",  # Allow having the presence field written before parameters' field
        "parameters": {},
    }

    if not elem.get("presence"):
        raise KeyError(f"In {data['name']} of MacSyFinder, one gene doesn't have a presence requirement")
    # Handle both Series (single match) and DataFrame (multiple matches)
    if isinstance(hmm_info, pd.Series):
        dict_elem["name"] = hmm_info["protein_name"]
        exchangeable_set = {hmm_info.name, hmm_info["secondary_name"]}
    else:  # DataFrame case - multiple matches
        # Use mode to get the most common protein name
        prot_name_list = hmm_info["protein_name"].mode().to_list()
        secondary_name_list = hmm_info["secondary_name"].mode().to_list()
        dict_elem["name"] = prot_name_list[0]
        exchangeable_set = set(prot_name_list) | {gene_name} | set(secondary_name_list)

    if not elem.get("presence"):
        raise KeyError(f"In {data['name']}, one gene doesn't have a presence requirement")
    # Process gene attributes
    dict_elem = process_attributes(elem, dict_elem, transitivity_mut)

    # Process exchangeable genes from child elements
    exchangeable_set = process_exchangeable(elem, data, hmm_df, exchangeable_set)

    # Set the exchangeable genes list (excluding self and empty strings)
    if exchangeable_set:
        primary_name = dict_elem["name"]
        dict_elem["exchangeable"] = list(exchangeable_set.difference({primary_name, ""}))

    return dict_elem


def translate_functional_unit(
    elem: lxml.etree.Element,
    data: Dict[str, Union[str, List, Dict[str, int]]],
    hmm_df: pd.DataFrame,
    transitivity_mut: Callable[[int], int],
) -> Dict[str, Union[str, List, Dict[str, int]]]:
    """
    Translate a functional unit from MacSyFinder models (or like) into PANORAMA format.

    Functional units represent collections of genes/families that work together as a unit.

    Args:
        elem (lxml.etree.Element): XML element corresponding to the functional unit
        data (Dict[str, Any]): Dictionary containing PANORAMA model information
        hmm_df (pd.DataFrame): HMM information DataFrame for gene translation
        transitivity_mut (Callable[[int], int]): Function to transform transitivity values

    Returns:
        Dict[str, Any]: Dictionary containing PANORAMA functional unit model information with keys:
            - name: functional unit name
            - presence: presence requirement
            - families: list of gene families in the unit
            - parameters: dict with various thresholds and flags

    Raises:
        KeyError: If the functional unit name is missing
    """
    name = elem.get("name")
    if not name:
        raise KeyError(f"In {data['name']} of MacSyFinder, one functional unit doesn't have a name")

    # Initialize functional unit dictionary
    dict_elem = {"name": name, "families": [], "parameters": {}}

    if not elem.get("presence"):
        raise KeyError(f"In {data['name']}, one functional unit doesn't have a presence requirement")
    # Process functional unit attributes
    dict_elem = process_attributes(elem, dict_elem, transitivity_mut)

    # Translate all gene families within this functional unit
    for family in elem:
        dict_elem["families"].append(translate_gene(family, data, hmm_df, transitivity_mut))

    return dict_elem


def translate_macsyfinder_model(
    root: et.Element,
    model_name: str,
    hmm_df: pd.DataFrame,
    canonical: List[str],
    transitivity_mut: Callable[[int], int],
) -> Dict[str, Union[str, List, Dict[str, int]]]:
    """
    Translate a complete MacSyFinder model (or like) into PANORAMA format.

    This function processes the root XML element of a MacSyFinder model and converts
    it to PANORAMA format, handling both gene-only models and models with functional units.

    Args:
        root (et.Element): Root XML element of the model
        model_name (str): Name of the model being translated
        hmm_df (pd.DataFrame): HMM information DataFrame for gene translation
        canonical (List[str]): List of canonical model names related to this model
        transitivity_mut (Callable[[int], int]): Function to transform transitivity values

    Returns:
        Dict[str, Any]: Complete PANORAMA model information with keys:
            - name: model name
            - func_units: list of functional units
            - parameters: model-level parameters
            - canonical: list of canonical models (if any)

    Raises:
        ModelTranslationError: For errors during model translation
    """
    # Initialize model dictionary with default structure
    data_json = {"name": model_name, "func_units": [], "parameters": {}}

    # Add canonical models if specified
    if canonical:
        data_json["canonical"] = canonical

    # Process root-level attributes
    if root.attrib:
        data_json = process_attributes(root, data_json, transitivity_mut)

    # Separate genes and functional units
    fu_list = []
    fam_list = []

    for elem in root:
        if elem.tag == "gene":
            fam_list.append(translate_gene(elem, data_json, hmm_df, transitivity_mut))
        elif elem.tag == "functional_unit":
            fu_list.append(translate_functional_unit(elem, data_json, hmm_df, transitivity_mut))
        else:
            logging.getLogger("PANORAMA").warning(f"Unknown element tag: {elem.tag}")

    # Handle model structure based on content
    if not fu_list:  # Gene-only model
        # Create a single functional unit containing all genes
        data_json["func_units"].append(
            {
                "name": data_json["name"],
                "presence": "mandatory",
                "families": fam_list,
                "parameters": data_json["parameters"],
            }
        )
        # Reset model parameters to defaults for gene-only models
        data_json["parameters"] = {
            "transitivity": 0,
            "window": 1,
            "min_mandatory": 1,
            "min_total": 1,
        }
    else:
        # Model with explicit functional units
        data_json["func_units"].extend(fu_list)

    return data_json


def parse_macsyfinder_hmm(hmm: HMM, hmm_file: Path, panorama_acc: Set[str]) -> Dict[str, Union[str, int, float]]:
    """
    Parse a MacSyFinder HMM and extract information for PANORAMA annotation.

    This function processes HMM objects and extracts metadata needed for PANORAMA,
    including generating accession numbers and handling naming conventions.

    Args:
        hmm (HMM): HMM object from the HMM file
        hmm_file (Path): Path to the HMM file being processed
        panorama_acc (Set[str]): Set of existing PANORAMA accession IDs to avoid duplicates

    Returns:
        Dict[str, Union[str, int, float]]: Dictionary with HMM metadata including:
            - name: HMM name
            - accession: unique accession ID
            - length: consensus sequence length
            - protein_name: protein name
            - secondary_name: alternative name
            - score_threshold: score threshold (NaN if not set)
            - eval_threshold: E-value threshold
            - ieval_threshold: independent E-value threshold
            - hmm_cov_threshold: HMM coverage threshold (NaN if not set)
            - target_cov_threshold: target coverage threshold (NaN if not set)
            - description: HMM description

    Raises:
        IOError: If HMM name cannot be determined
    """
    # Initialize HMM dictionary with default values
    hmm_dict = {
        "name": "",
        "accession": "",
        "length": len(hmm.consensus),
        "protein_name": "",
        "secondary_name": "",
        "score_threshold": nan,
        "eval_threshold": DEFAULT_EVAL_THRESHOLD,
        "ieval_threshold": DEFAULT_IEVAL_THRESHOLD,
        "hmm_cov_threshold": nan,
        "target_cov_threshold": nan,
        "description": "",
    }

    # Process HMM name
    hmm_dict = process_hmm_name(hmm, hmm_file, hmm_dict)

    # Generate or extract accession ID
    hmm_dict = process_hmm_accession(hmm, panorama_acc, hmm_dict)

    # Extract description if available
    if hmm.description is not None:
        hmm_dict["description"] = hmm.description.decode("UTF-8")

    return hmm_dict


def create_macsyfinder_hmm_list(
    hmms_path: Path,
    output: Path,
    binary_hmm: bool = False,
    hmm_coverage: float = None,
    target_coverage: float = None,
    force: bool = False,
    disable_bar: bool = False,
) -> pd.DataFrame:
    """
    Create a comprehensive HMM list file for PANORAMA annotation from MacSyFinder HMM files.

    This function processes all HMM files in a directory, extracts the metadata, and creates
    a standardized HMM list file for use in PANORAMA annotation steps.

    Args:
        hmms_path (Path): Path to the directory containing HMM files
        output (Path): Path to the output directory for generated files
        binary_hmm (bool, optional): Whether to write HMMs in binary format. Defaults to False.
        hmm_coverage (float, optional): Global HMM coverage threshold override. Defaults to None.
        target_coverage (float, optional): Global target coverage threshold override. Defaults to None.
        force (bool, optional): Whether to overwrite the existing output directory. Defaults to False.
        disable_bar (bool, optional): Whether to disable the progress bar. Defaults to False.

    Returns:
        pd.DataFrame: DataFrame containing HMM information indexed by name, with columns:
            - accession: unique accession ID
            - path: path to processed HMM file
            - length: HMM consensus length
            - protein_name: protein name
            - secondary_name: alternative name
            - score_threshold: score threshold
            - eval_threshold: E-value threshold
            - ieval_threshold: independent E-value threshold
            - hmm_cov_threshold: HMM coverage threshold
            - target_cov_threshold: target coverage threshold
            - description: HMM description

    Raises:
        IOError: If HMM files cannot be read or processed
    """
    logging.getLogger("PANORAMA").info("Starting HMM file translation from MacSyFinder...")

    # Create HMM output directory
    hmm_dir = mkdir(output / "hmm", force, erase=True)

    hmm_info_list = []
    panorama_acc = set()

    # Process all HMM files
    hmm_files = list(hmms_path.glob("*.hmm"))
    for hmm_file in tqdm(hmm_files, unit="HMM", desc="Translating HMMs", disable=disable_bar):
        try:
            hmm_list = read_hmm(hmm_file)
            for hmm in hmm_list:
                hmm_dict = parse_macsyfinder_hmm(hmm, hmm_file, panorama_acc)
                hmm_dict["path"] = write_hmm(hmm, hmm_dir, binary_hmm).absolute().as_posix()
                hmm_info_list.append(hmm_dict)
        except Exception as e:
            logging.getLogger("PANORAMA").error(f"Error processing HMM file {hmm_file}: {e}")
            raise

    # Create and format DataFrame
    hmm_df = pd.DataFrame(hmm_info_list)
    hmm_df = hmm_df.sort_values(by=["name", "accession", "protein_name"], ascending=True)

    # Reorder columns for consistency
    column_order = [
        "name",
        "accession",
        "path",
        "length",
        "protein_name",
        "secondary_name",
        "score_threshold",
        "eval_threshold",
        "ieval_threshold",
        "hmm_cov_threshold",
        "target_cov_threshold",
        "description",
    ]
    hmm_df = hmm_df[column_order]

    # Apply global coverage thresholds if specified
    if hmm_coverage is not None:
        logging.getLogger("PANORAMA").warning("Overriding HMM coverage thresholds globally")
        hmm_df["hmm_cov_threshold"] = hmm_coverage
    if target_coverage is not None:
        logging.getLogger("PANORAMA").warning("Overriding target coverage thresholds globally")
        hmm_df["target_cov_threshold"] = target_coverage

    # Write a HMM list file and prepare for return
    hmm_df.to_csv(output / "hmm_list.tsv", sep="\t", index=False)
    hmm_df.set_index("name", inplace=True)

    logging.getLogger("PANORAMA").info("HMM list file created successfully.")
    return hmm_df


def find_type_subtype_canonical(model_name: str, models: Path) -> List[str]:
    """
    Find canonical models for Type-Subtype naming patterns.

    Args:
        model_name: Name of the model
        models: Path to models directory

    Returns:
        List of canonical model names

    TODO:
        Try to have a function more general and refactor
    """
    canonical_systems = []
    if re.search("-Type-", model_name):
        basename = re.split("-Type-", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        for canon_file in models.rglob(f"{base_class}*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}-Subtype-{base_type}-", name_canon) and name_canon != model_name:
                logging.getLogger("PANORAMA").debug("")
                canonical_systems.append(name_canon)

    elif re.search("_Type_", model_name):
        basename = re.split("_Type_", model_name)
        base_class = basename[0]
        base_type = basename[-1]
        if len(base_type.split("_")) == 2:
            base_type, subtype = base_type.split("_")
        else:
            subtype = ""
        for canon_file in models.rglob("*.xml"):
            name_canon = canon_file.stem
            if re.search(f"{base_class}_Type_{base_type}_{subtype}", name_canon) and name_canon != model_name:
                logging.getLogger("PANORAMA").debug("")
                canonical_systems.append(name_canon)

    return canonical_systems


def find_cluster_canonical(model_name: str, models: Path) -> List[str]:
    """
    Find canonical models for cluster-based systems.

    Args:
        model_name: Name of the cluster model
        models: Path to models directory

    Returns:
        List of canonical model names
    """
    canonical_systems = []

    for canon_file in models.rglob(f"{model_name}*.xml"):
        if canon_file.stem != model_name:
            canonical_systems.append(canon_file.stem)
            logging.getLogger("PANORAMA").debug(f"Found cluster canonical: {canon_file.stem}")

    return canonical_systems


def search_canonical_macsyfinder(model_name: str, models: Path) -> List[str]:
    """
    Search for canonical models related to a given MacSyFinder model.

    This function identifies canonical (parent/related) models based on naming patterns
    specific to different bioinformatics tools and their model hierarchies.

    Args:
        model_name (str): Name of the model to find canonical models for
        models (Path): Path to the directory containing all model files

    Returns:
        List[str]: List of canonical model names related to the input model.
                  Empty list if no canonical models are found.

    Note:
        Canonical models are identified based on naming patterns:
        - Type/Subtype relationships (Defense systems)
        - Cluster relationships (CAS, CBASS, etc.)
        - Family relationships with underscore notation
    """
    canonical_systems = []

    # Handle Type-Subtype relationships (hyphen-separated)
    if "-Type-" in model_name:
        canonical_systems.extend(find_type_subtype_canonical(model_name, models))

    # Handle Type_Subtype relationships (underscore-separated)
    elif "_Type_" in model_name:
        canonical_systems.extend(find_type_subtype_canonical(model_name, models))

    # Handle specific system clusters
    elif model_name in ["CAS_Cluster", "CBASS", "Wadjet", "BREX"]:
        canonical_systems.extend(find_cluster_canonical(model_name, models))

    return canonical_systems


def get_models_path(models: Path, source: str) -> Dict[str, Path]:
    """
    Create a mapping of model names to their file paths based on the source tool.

    Different bioinformatics tools have different naming conventions and directory
    structures. This function standardizes the model name extraction process.

    Args:
        models (Path): Path to models database directory
        source (str): Name of the source tool ("defense-finder", "CONJScan", "TXSScan", "TFFscan")

    Returns:
        Dict[str, Path]: Dictionary mapping standardized model names to their file paths

    Raises:
        ValueError: If the source tool is not supported
    """
    model2path = {}

    # Define source-specific naming strategies
    naming_strategies = {
        "defense-finder": lambda model: model.stem,
        "CONJScan": lambda model: f"{model.parent.stem}_{model.stem}",
        "TXSScan": lambda model: model.stem,
        "TFFscan": lambda model: model.stem,
    }

    if source not in naming_strategies:
        raise ValueError(f"Unsupported source: {source}. Supported sources: {list(naming_strategies.keys())}")

    naming_func = naming_strategies[source]

    # Process all XML model files
    for model_file in models.rglob("*.xml"):
        try:
            model_name = naming_func(model_file)
            model2path[model_name] = model_file
        except Exception as e:
            logging.getLogger("PANORAMA").warning(f"Could not process model {model_file}: {e}")

    logging.getLogger("PANORAMA").info(f"Found {len(model2path)} models for source: {source}")
    return model2path


def translate_macsyfinder(
    macsy_db: Path,
    output: Path,
    binary_hmm: bool = False,
    hmm_coverage: float = None,
    target_coverage: float = None,
    source: str = "",
    force: bool = False,
    disable_bar: bool = False,
) -> List[Dict[str, Union[str, List, Dict[str, int]]]]:
    """
    Translate MacSyFinder models into PANORAMA format with comprehensive error handling.

    This is the main translation function that orchestrates the entire process of
    converting MacSyFinder-compatible models into PANORAMA format, including HMM
    processing, model translation, and file generation.

    Args:
        macsy_db (Path): Path to the MacSyFinder models database directory
        output (Path): Path to output directory for PANORAMA files
        binary_hmm (bool, optional): Whether to write HMMs in binary format. Defaults to False.
        hmm_coverage (float, optional): Global HMM coverage threshold. Defaults to source-specific values.
        target_coverage (float, optional): Global target coverage threshold. Defaults to None.
        source (str, optional): Name of the source tool. Defaults to "".
        force (bool, optional): Whether to overwrite existing files. Defaults to False.
        disable_bar (bool, optional): Whether to disable progress bars. Defaults to False.

    Returns:
        List[Dict[str, Any]]: List of dictionaries containing translated PANORAMA models

    Raises:
        ValueError: If an unsupported source is provided
        ModelTranslationError: If model translation fails
        IOError: If file operations fail
    """

    # Define transitivity mutation functions for different sources
    def default_transitivity_mut(x: int) -> int:
        """Default transitivity transformation - no change."""
        return x

    def conjscan_transitivity_mut(x: int) -> int:
        """CONJScan-specific transitivity transformation - divide by 2 and round up."""
        return int(x / 2) if int(x / 2) >= x / 2 else int(x / 2) + 1

    # Set the default coverage threshold based on the source
    if hmm_coverage is None:
        hmm_coverage = COVERAGE_THRESHOLDS.get(source, COVERAGE_THRESHOLDS["default"])

    # Process HMM files and create HMM DataFrame
    logging.getLogger("PANORAMA").info(f"Processing HMM files for {source}")
    hmm_df = create_macsyfinder_hmm_list(
        macsy_db / "profiles",
        output,
        binary_hmm,
        hmm_coverage,
        target_coverage,
        force,
        disable_bar,
    )

    # Select the appropriate transitivity mutation function
    transitivity_mut = conjscan_transitivity_mut if source == "CONJScan" else default_transitivity_mut

    # Translate all models
    list_data = []
    logging.getLogger("PANORAMA").info(f"Starting translation of {source} models")

    model_path = macsy_db / "definitions"
    model_mapping = get_models_path(model_path, source)

    for name, path in tqdm(
        model_mapping.items(),
        unit="model",
        desc="Translating models",
        disable=disable_bar,
    ):
        try:
            # Find canonical models for this model
            canonical_sys = search_canonical_macsyfinder(name, model_path)

            # Parse XML and translate model
            root = read_xml(path)
            translated_model = translate_macsyfinder_model(root, name, hmm_df, canonical_sys, transitivity_mut)
            list_data.append(translated_model)

        except Exception as error:
            logging.getLogger("PANORAMA").error(f"Failed to translate model {name} from {path.as_posix()}")
            raise ModelTranslationError(f"Error translating {name}: {error}") from error

    logging.getLogger("PANORAMA").info(f"Successfully translated {len(list_data)} models")
    return list_data
