#!/usr/bin/env python3
# coding:utf-8

"""
This module contains functions for managing files and directories, and checking the sanity of a TSV file.
"""

# default libraries
import sys
import os
import argparse
import logging
from typing import Dict, TextIO, Set, Union
import shutil
from pathlib import Path
from multiprocessing import Manager, Lock
from importlib.metadata import distribution
from pathlib import Path

# installed libraries
import numpy as np
import pandas as pd


# Create a custom formatter that combines both
class RawTextArgumentDefaultsHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def check_log(name: str) -> TextIO:
    """Check if the output log is writable

    :param name: Path to the log output

    :return: file object to write log
    """
    if name == "stdout":
        return sys.stdout
    elif name == "stderr":
        return sys.stderr
    else:
        log_file = Path(name)
        if log_file.exists():
            logging.getLogger("PANORAMA").warning(
                f"{log_file} exists, and will be overwriting"
            )
        return log_file


def pop_specific_action_grp(
    sub: argparse.ArgumentParser, title: str
) -> argparse._SubParsersAction:
    existing_titles = []

    for action_group in sub._action_groups:
        existing_titles.append(action_group.title)

        if action_group.title == title:
            sub._action_groups.remove(action_group)
            return action_group

    raise KeyError(
        f"{title} is not found in the provided subparser. Subparser contains {existing_titles}"
    )


def add_common_arguments(subparser: argparse.ArgumentParser) -> None:
    """
    Add common argument to the input subparser.

    :param subparser: A subparser object from any subcommand.
    """
    try:
        common = pop_specific_action_grp(
            subparser, "Optional arguments"
        )  # get the 'optional arguments' action group
    except KeyError:
        common = subparser.add_argument_group(title="Optional arguments")
    common.add_argument(
        "--verbose",
        required=False,
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)",
    )
    common.add_argument(
        "--log",
        required=False,
        type=check_log,
        default="stdout",
        help="log output file",
    )
    common.add_argument(
        "-d",
        "--disable_prog_bar",
        required=False,
        action="store_true",
        help="disables the progress bars",
    )
    common.add_argument(
        "--force",
        action="store_true",
        help="Force writing in output directory and in pangenome output file.",
    )
    subparser._action_groups.append(common)


def set_verbosity_level(args: argparse.Namespace) -> None:
    """Set the verbosity level

    :param args: argument pass by command line
    """
    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

        if (
            args.log != sys.stdout and not args.disable_prog_bar
        ):  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True
        str_format = "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s"
        datefmt = "%Y-%m-%d %H:%M:%S"
        if args.log in [sys.stdout, sys.stderr]:
            # use stream
            logging.basicConfig(
                stream=args.log, level=level, format=str_format, datefmt=datefmt
            )
        else:
            # log is written in a files. basic condif uses filename
            logging.basicConfig(
                filename=args.log, level=level, format=str_format, datefmt=datefmt
            )
        logging.getLogger("PANORAMA").debug(
            "Command: " + " ".join([arg for arg in sys.argv])
        )
        logging.getLogger("PANORAMA").debug(
            f"PANORAMA version: {distribution('panorama').version}"
        )


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
    if not isinstance(output, Path):
        raise Exception(f"Invalid directory type: {type(output)}")
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError:
        if not force:
            raise FileExistsError(
                f"{output} already exists."
                f"Use --force if you want to overwrite the files in the directory"
            )
        else:
            if erase:
                if any(output.iterdir()):
                    logging.getLogger("PANORAMA").warning(
                        f"Erasing the non empty directory: {output}"
                    )
                try:
                    shutil.rmtree(output)
                except Exception:
                    raise Exception(
                        f"It's not possible to remove {output}. Could be due to read-only files."
                    )
                else:
                    return mkdir(output, force=force, erase=erase)
            else:
                logging.getLogger("PANORAMA").warning(
                    f"{output.as_posix()} already exist and file could be overwrite by the new generated"
                )
                return Path(output)
    except Exception:
        raise Exception("An unexpected error happened. Please report on our GitHub")
    else:
        return Path(output)


def is_empty(filepath):
    """
    Checks if a file is empty.

    Args:
        filepath (str): The path to the file to check.

    Returns:
        bool: True if the file is empty, False otherwise.
    """
    return os.path.getsize(filepath) == 0


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
    tsv = pd.read_csv(tsv_path, sep="\t", header=None)
    if tsv.shape[1] < 2:
        raise SyntaxError(
            "Format not readable. You need at least 2 columns (name and path to pangenome)"
        )
    else:
        col_names = ["sp", "path", "taxid"]
        for col_idx in range(tsv.shape[1], len(col_names)):
            tsv[col_names[col_idx]] = None
    tsv.columns = col_names
    if tsv["sp"].isnull().values.any():
        err_df = pd.DataFrame(
            [tsv["sp"], tsv["sp"].isnull()], ["pangenome_name", "error"]
        ).T
        logging.getLogger("PANORAMA").error(
            "\n" + err_df.to_string().replace("\n", "\n\t")
        )
        raise ValueError(
            "There is a line with no value in pangenome name (first column)"
        )
    else:
        if tsv["sp"].str.count(" ").any():
            err_df = pd.DataFrame(
                [tsv["sp"], np.where(tsv["sp"].str.count(" ") >= 1, True, False)],
                ["pangenome_name", "error"],
            ).T
            logging.getLogger("PANORAMA").error(
                "\n" + err_df.to_string().replace("\n", "\n\t")
            )
            raise ValueError(
                "Your pangenome names contain spaces. "
                "To ensure compatibility with all of the dependencies this is not allowed. "
                "Please remove spaces from your pangenome names."
            )
    if tsv["path"].isnull().values.any():
        err_df = pd.DataFrame(
            [tsv["path"], tsv["path"].isnull()], ["pangenome_path", "error"]
        ).T
        logging.getLogger("PANORAMA").error(
            "\n" + err_df.to_string().replace("\n", "\n\t")
        )
        raise ValueError("There is a line with no path (second column)")
    else:
        if not tsv["path"].map(lambda x: Path(x).exists()).all():
            err_df = pd.DataFrame(
                [tsv["path"], ~tsv["path"].map(lambda x: Path(x).exists())],
                ["pangenome_path", "error"],
            ).T
            logging.getLogger("PANORAMA").error(
                "\n" + err_df.to_string().replace("\n", "\n\t")
            )
            raise FileNotFoundError(
                "Unable to locate one or more pangenome in your file.}"
            )
        tsv["path"] = tsv["path"].map(lambda x: Path(x).resolve().absolute())

        return tsv.set_index("sp").to_dict("index")


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


def is_true_value(value: Union[str, int, bool]) -> bool:
    """
    Check if a value represents a true condition.

    True conditions are: 1, "True" or True

    Args:
        value: Value to check

    Returns:
        True if value represents a true condition
    """
    return value == 1 or value == "True" or value is True
