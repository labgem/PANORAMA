#!/usr/bin/env python3
# coding:utf-8

"""
Padloc Model Translation Module for PANORAMA

This module provides functionality to translate models from PADLOC into
PANORAMA-compatible formats. It handles HMM processing, metadata parsing and
model structure conversion.
"""

# default libraries
import logging
import re
from pathlib import Path
from typing import Dict, List, Set, Union

# installed libraries
import pandas as pd
from tqdm import tqdm

# local libraries
from panorama.utility.genInput import create_hmm_list_file
from panorama.utility.translate import ModelTranslationError, read_yaml

PADLOC_MODEL_KEYS = [
    "maximum_separation",
    "minimum_core",
    "minimum_total",
    "force_strand",
    "core_genes",
    "secondary_genes",
    "neutral_genes",
    "prohibited_genes",
]


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
        temp_val = row.get("temp", "")
        secondary_val = row.get("secondary.name", "")

        if pd.isna(temp_val) and pd.isna(secondary_val):
            return ""
        elif pd.isna(temp_val):
            return str(secondary_val)
        elif pd.isna(secondary_val):
            return str(temp_val)
        else:
            return f"{temp_val},{secondary_val}"

    # Define expected column structure
    meta_columns = [
        "accession",
        "name",
        "protein_name",
        "secondary_name",
        "score_threshold",
        "eval_threshold",
        "ieval_threshold",
        "hmm_cov_threshold",
        "target_cov_threshold",
        "description",
    ]

    # Read the TSV file with specific columns
    df = pd.read_csv(meta_path, sep="\t", usecols=[0, 1, 2, 3, 4, 6, 7, 8], header=0)

    # Parse protein names (format: main_name|secondary_name)
    df[["protein_name", "temp"]] = df["protein.name"].str.split("|", expand=True)

    # Merge secondary names
    df["secondary_name"] = df.apply(_merge_secondary_names, axis=1)

    # Clean up temporary columns
    df = df.drop(columns=["protein.name", "secondary.name", "temp"], errors="ignore")

    # Reorder columns and add missing ones
    df = df.iloc[:, [0, 1, 6, 7, 3, 4, 5, 2]]  # Reorder existing columns

    # Insert missing threshold columns
    df.insert(4, "score_threshold", None)
    df.insert(5, "eval_threshold", None)

    # Set proper column names
    df.columns = meta_columns

    # Set index and handle missing descriptions
    df = df.set_index("accession")
    df["description"] = df["description"].fillna("unknown")

    return df


def _add_families_to_functional_unit(
    families_list: List[str],
    secondary_names: List[str],
    family_type: str,
    metadata_df: pd.DataFrame,
    seen_families: Set[str],
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
        if family_name in seen_families or family_name in ["NA", "cas_accessory"]:
            continue

        # Special handling for cas_adaptation families
        if family_name == "cas_adaptation":
            filtered_df = metadata_df.loc[metadata_df["secondary_name"] == family_name]
            for protein_name in filtered_df["protein_name"].dropna().unique():
                if protein_name not in seen_families:
                    family_dict = {"name": protein_name, "presence": family_type}
                    seen_families.add(protein_name)

                    # Add exchangeable proteins if available
                    if protein_name in secondary_names:
                        exchangeable_proteins = filtered_df["protein_name"].dropna().unique().tolist()
                        if len(exchangeable_proteins) > 1:
                            family_dict["exchangeable"] = exchangeable_proteins

                    family_list.append(family_dict)
        else:
            # Standard family processing
            family_dict = {"name": family_name, "presence": family_type}
            seen_families.add(family_name)

            # Add exchangeable proteins if this family has secondary names
            if family_name in secondary_names:
                filtered_df = metadata_df.loc[metadata_df["secondary_name"] == family_name]
                exchangeable_proteins = filtered_df["protein_name"].dropna().unique().tolist()
                if exchangeable_proteins:
                    family_dict["exchangeable"] = exchangeable_proteins

            family_list.append(family_dict)

    return family_list


def translate_model_padloc(
    data_yaml: Dict[str, Union[List[str], int, bool]],
    model_name: str,
    metadata_df: pd.DataFrame,
    canonical_models: List[str] = None,
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
            f"Unexpected keys in PADLOC model '{model_name}': {unexpected_keys}. Expected keys: {PADLOC_MODEL_KEYS}"
        )

    if "core_genes" not in data_yaml:
        raise KeyError(f"Missing required 'core_genes' key in PADLOC model '{model_name}'")

    # Initialize PANORAMA model structure
    panorama_model = {
        "name": model_name,
        "func_units": [],
        "parameters": {
            "transitivity": 0,
            "window": 1,
            "min_mandatory": 1,
            "min_total": 1,
        },
    }

    if canonical_models:
        panorama_model["canonical"] = canonical_models

    # Extract secondary names for exchangeable protein lookup
    secondary_names = metadata_df["secondary_name"].dropna().unique().tolist()

    # Create a functional unit
    functional_unit = {
        "name": model_name,
        "presence": "mandatory",
        "same_strand": data_yaml.get("force_strand", False),
        "parameters": {
            "transitivity": data_yaml.get("maximum_separation", 0),
            "min_mandatory": data_yaml.get("minimum_core", 1),
            "min_total": data_yaml.get("minimum_total", 1),
        },
    }

    # Process gene families by category
    family_list = []
    seen_families = set()

    # Core genes (mandatory)
    family_list.extend(
        _add_families_to_functional_unit(
            data_yaml["core_genes"],
            secondary_names,
            "mandatory",
            metadata_df,
            seen_families,
        )
    )

    # Secondary genes (accessory)
    if "secondary_genes" in data_yaml:
        family_list.extend(
            _add_families_to_functional_unit(
                data_yaml["secondary_genes"],
                secondary_names,
                "accessory",
                metadata_df,
                seen_families,
            )
        )

    # Prohibited genes (forbidden)
    if "prohibited_genes" in data_yaml:
        family_list.extend(
            _add_families_to_functional_unit(
                data_yaml["prohibited_genes"],
                secondary_names,
                "forbidden",
                metadata_df,
                seen_families,
            )
        )

    # Neutral genes
    if "neutral_genes" in data_yaml:
        family_list.extend(
            _add_families_to_functional_unit(
                data_yaml["neutral_genes"],
                secondary_names,
                "neutral",
                metadata_df,
                seen_families,
            )
        )

    functional_unit["families"] = family_list
    panorama_model["func_units"].append(functional_unit)

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
    disable_bar: bool = False,
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
        disable_bar=disable_bar,
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
            translated_model = translate_model_padloc(model_data, model_file.stem, metadata_df, canonical_models)
            translated_models.append(translated_model)

        except Exception as e:
            logging.getLogger("PANORAMA").error(f"Failed to translate PADLOC model {model_file.stem}: {e}")
            raise ModelTranslationError(f"Translation failed for {model_file.stem}") from e

    logging.getLogger("PANORAMA").info(f"Successfully translated {len(translated_models)} PADLOC models")
    return translated_models
