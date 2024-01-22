#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import Dict, List, Set, Union
from pathlib import Path
from random import choice
from string import digits

# install libraries
from tqdm import tqdm
import pandas as pd
from numpy import nan


def read_metadata(metadata: Path) -> pd.DataFrame:
    """ Read metadata associate with HMM

    Args:
        metadata (Path): path to the metadata file

    Raises:
        FileNotFoundError: If metadata path is not found.
        IOError: If the metadata path is not a file
        ValueError: If the number of field is unexpected
        NameError: If the column names use in metadata are not allowed

    Returns:
        str: metadata dataframe with hmm information
    """
    logging.debug("Reading HMM metadata...")
    authorize_names = ["accession", "name", "protein_name", "secondary_name", "score_threshold",
                       "eval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
    dtype = {"accession": "string", "name": "string", "protein_name": "string", "secondary_name": "string",
             "score_threshold": "float", "eval_threshold": "float", "hmm_cov_threshold": "float",
             "target_cov_threshold": "float", "description": "string"}
    if not metadata.exists():
        raise FileNotFoundError(f"Metadata file does not exist at the given path: {metadata}")
    if not metadata.is_file():
        raise IOError(f"Metadata path is not a file: {metadata}")
    metadata_df = pd.read_csv(metadata, delimiter="\t", header=0)
    if metadata_df.shape[1] == 1 or metadata_df.shape[1] > len(authorize_names):
        raise ValueError("The number of field is unexpected. Please check that tabulation is used as separator and "
                         "that you give not more than the expected columns")
    if any(name not in authorize_names for name in metadata_df.columns):
        logging.getLogger("PANORAMA").error(f"Authorized keys: {authorize_names}")
        logging.getLogger("PANORAMA").debug(f"metadata_df.columns.names: {metadata_df.columns}")
        raise NameError("The column names use in metadata are not allowed")
    metadata_df.astype(dtype)
    metadata_df = metadata_df.set_index('accession')
    metadata_df['description'] = metadata_df["description"].fillna('unknown')
    return metadata_df


def gen_acc(acc: str, panorama_acc: Set[str]) -> str:
    """
    Generates a unique accession number for the given HMM.

    Args:
        acc (str): The accession number to check.
        panorama_acc (Set[str]): The set of existing accession numbers.

    Returns:
        str: A unique accession number.
    """
    if acc in panorama_acc:
        return gen_acc("PAN" + ''.join(choice(digits) for _ in range(6)), panorama_acc)
    else:
        panorama_acc.add(acc)
        return acc


def read_hmm(hmm_file: Path, metadata: pd.DataFrame = None) -> Dict[str, Union[str, int, float]]:
    """
    Reads the given HMM file and returns a dictionary containing information about the HMM.

    Args:
        hmm_file (Path): The path to the HMM file.
        metadata (pd.DataFrame, optional): The metadata dataframe. Defaults to None.

    Returns:
        Dict[str, Union[str, int, float]]: A dictionary containing information about the HMM.
    """
    hmm_dict = {"name": "", 'accession': "", 'path': hmm_file, "length": nan, "description": ""}
    stop = False
    with open(hmm_file, "r") as hmm:
        hmm.readline()  # Skip first line
        while not stop:
            line = hmm.readline()
            if line.startswith("HMM"):
                stop = True
            else:
                if line.startswith("NAME"):
                    hmm_dict["name"] = "_".join(line.split()[1:])
                elif line.startswith("ACC"):
                    hmm_dict["accession"] = line.split()[1]
                elif line.startswith("DESC"):
                    hmm_dict["description"] = "_".join(line.split()[1:])
                elif line.startswith("LENG"):
                    hmm_dict["length"] = int(line.split()[1])
    if metadata is not None and hmm_dict["accession"] in metadata.index:
        hmm_info = metadata.loc[hmm_dict["accession"]]
        hmm_dict.update(hmm_info.to_dict())
    else:
        hmm_dict.update({'protein_name': "", 'secondary_name': "", "score_threshold": nan, "eval_threshold": nan,
                         "hmm_cov_threshold": nan, "target_cov_threshold": nan})
    return hmm_dict


def create_hmm_list_file(hmm_path: List[Path], output: Path, metadata_df: pd.DataFrame = None,
                         recursive: bool = False, disable_bar: bool = False) -> None:
    """
    Creates a TSV file containing information about the given HMM files.

    Args:
        hmm_path (List[Path]): The paths to the HMM files.
        output (Path): The path to the output directory.
        metadata_df (pd.DataFrame, optional): The metadata dataframe. Defaults to None.
        recursive (bool, optional): Whether to search for HMM files recursively in the given directory. Defaults to False.
        disable_bar (bool, optional): Whether to disable the progress bar. Defaults to False.

    Raises:
        FileNotFoundError: If any of the given paths are not found.
        Exception: If an unexpected error occurs.

    Returns:
        None
    """
    logging.getLogger("PANORAMA").info("Begin to create hmm list file...")
    hmm_list = []
    hmm_path_list = []
    panorama_acc = set()
    for path in hmm_path:
        if path.is_file():
            hmm_path_list.append(path)
        elif path.is_dir():
            for hmm_file in path.rglob("*.hmm") if recursive else path.glob("*.hmm"):
                hmm_path_list.append(hmm_file)
        else:
            if not path.exists():
                raise FileNotFoundError(f"The given path is not found: {path}")
            else:
                raise Exception("Unexpected error")
    for hmm in tqdm(hmm_path_list, unit="HMM", disable=disable_bar):
        hmm_dict = read_hmm(hmm, metadata_df)
        if hmm_dict["accession"] == "":
            hmm_dict["accession"] = gen_acc("PAN" + ''.join(choice(digits) for _ in range(6)), panorama_acc)
        if hmm_dict["description"] == "":
            hmm_dict["description"] = 'unknown'
        hmm_list.append(hmm_dict)
    pd.DataFrame(hmm_list).to_csv(output / "hmm_list.tsv", sep="\t", index=False)
    logging.getLogger("PANORAMA").info("HMM list file created.")
