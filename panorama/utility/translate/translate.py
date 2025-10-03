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
- CasFinder: CRISPR-Cas Antiviral defense systems identification
- CONJScan: Conjugation system detection
- TXSScan: Type secretion system detection
- TFFScan: Type IV-A pilus detection
"""

# default libraries
from typing import Dict, List, Union
import logging
from pathlib import Path
from lxml import etree as et
import yaml
import json

# installed libraries
from tqdm import tqdm
import pandas as pd

# local libraries
from panorama.utils import mkdir


# Constants
KNOWN_SOURCES = ["padloc", "defense-finder", "CONJScan", "TXSScan", "TFFscan"]


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
        raise ModelTranslationError(
            f"Unexpected error reading {model_path}: {e}"
        ) from e


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
            huge_tree=False,
        )
        tree = et.parse(str(model_path.absolute()), parser=parser)
        return tree.getroot()
    except IOError as e:
        raise IOError(f"Problem opening {model_path}: {e}") from e
    except et.XMLSyntaxError as e:
        raise et.XMLSyntaxError(f"Problem parsing XML file {model_path}: {e}") from e
    except Exception as e:
        raise ModelTranslationError(
            f"Unexpected error reading {model_path}: {e}"
        ) from e


def write_model(
    output_path: Path,
    model_data: Dict[str, Union[str, List[Dict], Dict[str, int], List[str]]],
) -> Path:
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
    if "name" not in model_data:
        raise KeyError("Model data must contain 'name' field")

    output_file = output_path / f"{model_data['name']}.json"

    try:
        with open(output_file, "w", encoding="utf-8") as file:
            json.dump(model_data, file, indent=2, ensure_ascii=False)
        return output_file
    except IOError as e:
        raise IOError(
            f"Problem writing model {model_data['name']} to {output_file}: {e}"
        ) from e


def launch_translate(
    database_path: Path,
    source: str,
    output: Path,
    binary_hmm: bool = False,
    hmm_coverage: float = None,
    target_coverage: float = None,
    force: bool = False,
    disable_bar: bool = False,
):
    """
    Launch the complete model translation process for any supported source.

    This is the main entry point for translating models from different
    bioinformatics databases to PANORAMA format. It handles the entire
    workflow from parsing to writing output files.

    Args:
        database_path: Path to the source database directory
        source: Source identifier (padloc, defense-finder, CONJScan, TXSScan, TFFscan)
        output: Path to output directory for translated files
        binary_hmm: Whether to output HMMs in binary format (default: False)
        hmm_coverage: Global HMM coverage threshold (optional, defaults vary by source)
        target_coverage: Global target coverage threshold (optional)
        force: Whether to overwrite existing output files (default: False)
        disable_bar: Whether to disable progress bars (default: False)

    Raises:
        ValueError: If the source is not recognized
        FileNotFoundError: If the database path doesn't exist
        ModelTranslationError: If the translation process fails
    """

    if source == "padloc":
        from panorama.utility.translate.padloc_translator import translate_padloc

        list_data = translate_padloc(
            padloc_db=database_path,
            output=output,
            binary_hmm=binary_hmm,
            hmm_coverage=hmm_coverage,
            target_coverage=target_coverage,
            force=force,
            disable_bar=disable_bar,
        )
    elif source in ["defense-finder", "CONJScan", "TXSScan", "TFFscan"]:
        from panorama.utility.translate.macsymodel_translator import (
            translate_macsyfinder,
        )

        list_data = translate_macsyfinder(
            macsy_db=database_path,
            output=output,
            binary_hmm=binary_hmm,
            hmm_coverage=hmm_coverage,
            target_coverage=target_coverage,
            source=source,
            force=force,
            disable_bar=disable_bar,
        )
    else:
        raise ValueError(
            f"The given source: {source} is not recognize. "
            f"Please choose between padloc, defense-finder or macsy-finder"
        )
    logging.getLogger("PANORAMA").info("Write models for PANORAMA...")
    model_list = []
    model_dir = mkdir(output / "models", force)
    for data in tqdm(list_data, unit="model", desc="Write models", disable=disable_bar):
        model_path = write_model(model_dir, data)
        model_list.append([data["name"], model_path.resolve()])
    model_df = pd.DataFrame(model_list, columns=["name", "path"])
    model_df = model_df.sort_values("name")
    if source == "CONJScan":
        plasmid_mask = model_df["name"].str.contains("Plasmids_")
        chromosome_mask = model_df["name"].str.contains("Chromosome_")
        model_df[plasmid_mask].to_csv(
            output / "models_plasmids_list.tsv", sep="\t", header=False, index=False
        )
        model_df[chromosome_mask].to_csv(
            output / "models_chromosome_list.tsv", sep="\t", header=False, index=False
        )
    model_df.to_csv(output / "models_list.tsv", sep="\t", header=False, index=False)
