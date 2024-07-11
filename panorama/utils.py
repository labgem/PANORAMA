#!/usr/bin/env python3
# coding:utf-8

"""
This module contains functions for managing files and directories, and checking the sanity of a TSV file.
"""

# default libraries
import logging
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, Union
from multiprocessing import Manager, Lock


# File managing system
def mkdir(output: Path, force: bool = False, erase: bool = False) -> Path:
    """
    Create a directory at the given path.

    Args:
        output (Path): The path to the output directory.
        force (bool, optional): Whether to raise an exception if the directory already exists. Defaults to False.
        erase (bool, optional): Whether to erase the directory if it already exists and force is True. Defaults to False.

    Returns:
        Path: The path to the output directory.

    Raises:
        FileExistsError: If the directory already exists and force is False.
        Exception: If an unexpected error occurs.
    """
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError:
        if not force:
            raise FileExistsError(f"{output} already exists."
                                  f"Use --force if you want to overwrite the files in the directory")
        else:
            if erase:
                logging.getLogger("PANORAMA").warning(f"Erasing the directory: {output}")
                try:
                    shutil.rmtree(output)
                except Exception:
                    raise Exception(f"It's not possible to remove {output}. Could be due to read-only files.")
                else:
                    return mkdir(output, force=force, erase=erase)
            else:
                logging.getLogger("PANORAMA").warning(f"{output.as_posix()} already exist and file could be overwrite by the new generated")
                return Path(output)
    except Exception:
        raise Exception("An unexpected error happened. Please report on our GitHub")
    else:
        return Path(output)


def check_tsv_sanity(tsv_path: Path) -> Dict[str, Dict[str, Union[int, str, Path]]]:
    """
    Check if the given TSV file is readable for the next PANORAMA step.

    Args:
        tsv_path (Path): The path to the TSV file with the list of pangenomes.

    Returns:
        Dict[str, Dict[str, Union[int, str, Path]]]: A dictionary with pangenome name as key and a dictionary with path and taxid as values.

    Raises:
        SyntaxError: If the TSV file has less than 2 columns.
        ValueError: If there is a line with no value in pangenome name or if the pangenome names contain spaces.
        FileNotFoundError: If unable to locate one or more pangenomes in the TSV file.
    """
    tsv = pd.read_csv(tsv_path, sep='\t', header=None)
    if tsv.shape[1] < 2:
        raise SyntaxError("Format not readable. You need at least 2 columns (name and path to pangenome)")
    else:
        col_names = ['sp', 'path', 'taxid']
        for col_idx in range(tsv.shape[1], len(col_names)):
            tsv[col_names[col_idx]] = None
    tsv.columns = col_names
    if tsv['sp'].isnull().values.any():
        err_df = pd.DataFrame([tsv['sp'], tsv['sp'].isnull()], ["pangenome_name", "error"]).T
        logging.getLogger("PANORAMA").error("\n" + err_df.to_string().replace('\n', '\n\t'))
        raise ValueError("There is a line with no value in pangenome name (first column)")
    else:
        if tsv['sp'].str.count(" ").any():
            err_df = pd.DataFrame([tsv['sp'], np.where(tsv['sp'].str.count(" ") >= 1, True, False)],
                                  ["pangenome_name", "error"]).T
            logging.getLogger("PANORAMA").error("\n" + err_df.to_string().replace('\n', '\n\t'))
            raise ValueError("Your pangenome names contain spaces. "
                             "To ensure compatibility with all of the dependencies this is not allowed. "
                             "Please remove spaces from your pangenome names.")
    if tsv['path'].isnull().values.any():
        err_df = pd.DataFrame([tsv['path'], tsv['path'].isnull()], ["pangenome_path", "error"]).T
        logging.getLogger("PANORAMA").error("\n" + err_df.to_string().replace('\n', '\n\t'))
        raise ValueError("There is a line with no path (second column)")
    else:
        if not tsv["path"].map(lambda x: Path(x).exists()).all():
            err_df = pd.DataFrame([tsv['path'], ~tsv["path"].map(lambda x: Path(x).exists())],
                                  ["pangenome_path", "error"]).T
            logging.getLogger("PANORAMA").error("\n" + err_df.to_string().replace('\n', '\n\t'))
            raise FileNotFoundError("Unable to locate one or more pangenome in your file.}")
        tsv["path"] = tsv["path"].map(lambda x: Path(x).resolve().absolute())

        return tsv.set_index('sp').to_dict('index')


def init_lock(lock: Lock = None):
    """
    Initialize the loading lock.

    Args:
        lock (Lock, optional): The lock object to be assigned to `loading_lock`. Defaults to None.

    Returns:
        Lock: The lock object assigned to `loading_lock`.
    """
    if lock is None:
        manager = Manager()
        return manager.Lock()
