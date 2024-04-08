#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import Dict, List, Set, Tuple, Union
from pathlib import Path
from random import choice
from string import digits
import functools

# install libraries
from tqdm import tqdm
import pandas as pd
from numpy import nan
from pyhmmer.easel import Alphabet
from pyhmmer.plan7 import HMM, HMMFile

# local libraries
from panorama.utils import mkdir


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
    logging.getLogger("PANORAMA").debug("Reading HMM metadata...")
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


def read_hmm(hmm_path: Path) -> HMM:
    try:
        hmm_file = HMMFile(hmm_path)
    except Exception as error:
        logging.getLogger("PANORAMA").error(f"Problem reading HMM: {hmm_path}")
        raise error
    end = False
    while not end:
        try:
            hmm = next(hmm_file)
        except StopIteration:
            end = True
        except Exception as error:
            raise Exception(f'Unexpected error on HMM file {hmm_path}, caused by {error}')
        else:
            return hmm


def parse_hmm_info(hmm: HMM, panorama_acc: Set[str], metadata: pd.DataFrame = None):
    hmm_dict = {"name": hmm.name.decode('UTF-8'), 'accession': "", "length": len(hmm.consensus)}

    hmm.name = hmm_dict["name"].encode("UTF-8")

    if hmm.accession is None or hmm.accession == "".encode("UTF-8"):
        hmm_dict["accession"] = gen_acc("PAN" + ''.join(choice(digits) for _ in range(6)), panorama_acc)
        hmm.accession = hmm_dict["accession"].encode("UTF-8")
    else:
        hmm_dict["accession"] = hmm.accession.decode("UTF-8")

    if hmm.description is not None:
        hmm_dict["description"] = hmm.description.decode("UTF-8")

    if metadata is not None and hmm_dict["accession"] in metadata.index:
        hmm_info = metadata.loc[hmm_dict["accession"]]
        hmm_dict.update(hmm_info.to_dict())
    else:
        hmm_dict.update({'protein_name': "", 'secondary_name': "", "score_threshold": nan, "eval_threshold": nan,
                         "hmm_cov_threshold": nan, "target_cov_threshold": nan})
    return hmm_dict


def write_hmm(hmm: HMM, output: Path, binary: bool = False):
    with open(output/f"{hmm.name.decode('UTF-8')}.{'h3m' if binary else 'hmm'}", "wb") as file:
        hmm.write(file, binary)


def create_hmm_list_file(hmm_path: List[Path], output: Path, metadata_df: pd.DataFrame = None,
                         hmm_coverage: float = None, target_coverage: float = None, binary_hmm: bool = False,
                         recursive: bool = False, force: bool = False, disable_bar: bool = False) -> None:
    """
    Creates a TSV file containing information about the given HMM files.

    Args:
        hmm_path (List[Path]): The paths to the HMM files.
        output (Path): The path to the output directory.
        metadata_df (pd.DataFrame, optional): The metadata dataframe. Defaults to None.
        hmm_coverage: Set a global value of HMM coverage threshold for all HMM. Defaults to None
        target_coverage: Set a global value of target coverage threshold for all target. Defaults to None
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

    hmm_out = mkdir(output/"hmm", force, True) if binary_hmm else None
    for hmm_file in tqdm(hmm_path_list, unit="HMM", disable=disable_bar):
        hmm = read_hmm(hmm_path=hmm_file)
        hmm_dict = parse_hmm_info(hmm, panorama_acc, metadata_df)
        hmm_list.append(hmm_dict)
        if binary_hmm:
            write_hmm(hmm, hmm_out, True)

    hmm_df = pd.DataFrame(hmm_list)
    if hmm_coverage is not None:
        hmm_df["hmm_cov_threshold"] = hmm_coverage
    if target_coverage is not None:
        hmm_df["target_cov_threshold"] = target_coverage
    hmm_df.to_csv(output / "hmm_list.tsv", sep="\t", index=False)
    logging.getLogger("PANORAMA").info("HMM list file created.")
