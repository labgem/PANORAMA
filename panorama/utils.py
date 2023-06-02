#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import os
import sys
import logging
from pathlib import Path
import csv
from typing import TextIO, Dict, Union
import pkg_resources
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Lock
from tqdm import tqdm

# installed libraries
from ppanggolin.formats import check_pangenome_info

# local libraries
from panorama.pangenomes import Pangenome



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
        return open(name, "w")


def set_verbosity_level(args: argparse.Namespace):
    """Set the verbosity level

    :param args: argument pass by command line
    """
    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

        if args.log != sys.stdout and not args.disable_prog_bar:  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True

        logging.basicConfig(stream=args.log, level=level,
                            format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.getLogger().info("Command: " + " ".join([arg for arg in sys.argv]))
        logging.getLogger().info("Panorama version: " + pkg_resources.get_distribution("panorama").version)


# File managing system
def mkdir(output: str, force: bool = False) -> Path:
    """Create a directory at the given path

    :param output: Path to output directory
    :param force: Pass exception if directory already exist

    :raise FileExistError: If the given directory already exist and no force is used
    :raise Exception: Handle all others exception

    :return: Path object to output directory
    """
    try:
        os.makedirs(output)
    except OSError:
        if not force:
            raise FileExistsError(f"{output} already exists."
                                  f"Use --force if you want to overwrite the files in the directory")
        else:
            logging.getLogger().warning(f"The {output} directory already exists. The file it contains may be potentially overwritten by the newly generated results.")
            return Path(output)
    except Exception:
        raise Exception("An unexpected error happened. Please report on our GitHub")
    else:
        return Path(output)


def path_exist(path: Path) -> Path:
    try:
        abs_path = path.resolve()
        if not abs_path.exists():
            raise FileNotFoundError
    except FileNotFoundError:
        raise FileNotFoundError(f"{path.resolve()} not exist")
    except Exception:
        raise Exception(f"An unexpected error happened. Please report to our github.")
    else:
        return abs_path


def path_is_dir(path: Path) -> bool:
    abs_path = path_exist(path)
    if abs_path.is_dir():
        return True
    else:
        return False


def path_is_file(path: Path) -> bool:
    abs_path = path_exist(path)
    if abs_path.is_file():
        return True
    else:
        return False


def check_tsv_sanity(tsv_path: Path) -> Dict[str, Dict[str, Union[int, str]]]:
    """ Check if the given tsv is readable for the next PANORAMA step

    :param tsv_path: Path to tsv file with list of pangenome

    :raise IOError: If tsv or a pangenome not exist raise IOError
    :raise Exception: Handle all others exception

    :return: Dictionary with pangenome name as key and path to hdf5 file as value
    """
    pan_to_path = {}
    try:
        p_file = open(tsv_path.absolute(), 'r')
        tsv = csv.reader(p_file, delimiter="\t")
    except IOError as ios_error:
        raise IOError(ios_error)
    except Exception as exception_error:
        raise Exception(f"The following unexpected error happened when opening the list of pangenomes : "
                        f"{exception_error}")
    else:
        for line in tsv:
            if len(line) < 2:
                raise SyntaxError("Format not readable. You need at least 2 columns (name and path to pangenome)")
            if " " in line[0]:
                raise ValueError(f"Your pangenome names contain spaces (The first encountered pangenome name that had "
                                f"this string: '{line[0]}'). To ensure compatibility with all of the dependencies of "
                                f"PPanGGOLiN this is not allowed. Please remove spaces from your pangenome names.")
            try:
                abs_path = path_exist(Path(line[1]))
            except FileNotFoundError:
                try:
                    abs_path = path_exist(tsv_path.parent/Path(line[1]))
                except FileNotFoundError as file_error:
                    raise FileNotFoundError(f"{file_error}")
                else:
                    pan_to_path[line[0]] = {"path": f"{abs_path.as_posix()}",
                                            "taxid": line[2] if len(line) > 2 else None}
            except Exception:
                raise Exception("Unexpected error")
            else:
                pan_to_path[line[0]] = {"path": f"{abs_path.as_posix()}",
                                        "taxid": line[2] if len(line) > 2 else None}
        p_file.close()
        return pan_to_path


loading_lock = None


def init_lock(lock):
    """
    Initialize the loading lock.

    This function initializes the `loading_lock` variable as a global variable, assigning it the value of the `lock` parameter.
    If the `loading_lock` is already initialized, the function does nothing.

    :param lock: The lock object to be assigned to `loading_lock`.
    """
    global loading_lock
    if loading_lock is None:
        loading_lock = lock


def load_pangenome(name: str, path: str, taxid: int, need_info: dict) -> bool:
    """
    Load a pangenome from a given path and check the required information.

    This function loads a pangenome from the specified `path` and assigns it the provided `name` and `taxid`.
    The pangenome file is added to the pangenome object. The function then checks that the required information
    are present in the pangenome and if they are, it loads them.

    :param name: The name of the pangenome.
    :param path: The path to the pangenome file.
    :param taxid: The taxonomic ID associated with the pangenome.
    :param need_info: A dictionary containing information required to load in the Pangenome object.
    :return: True if the pangenome is loaded successfully and the required information is present.
    """
    pangenome = Pangenome(name=name, taxid=taxid)
    pangenome.add_file(path)

    check_pangenome_info(pangenome, disable_bar=True, **need_info)

    return True


def load_multiple_pangenomes(pan_name_to_path: Dict[str, Dict[str, Union[str, int]]], 
                             max_workers: int, disable_bar: bool, lock: Lock, need_info: bool):
    """
    Load multiple pangenomes in parallel using a process pool executor.

    This function loads multiple pangenomes in parallel using a process pool executor. It takes a dictionary 
    `pan_name_to_path` containing the mapping of pangenome names to their corresponding paths and other
    information. The pangenomes are loaded using the `load_pangenome` function. The loading progress is
    displayed using a tqdm progress bar.

    :param pan_name_to_path: A dictionary mapping pangenome names to their corresponding paths and information.
    :param max_workers: The maximum number of worker processes to use in the process pool executor.
    :param disable_bar: A flag indicating whether to disable the tqdm progress bar.
    :param lock: A multiprocessing lock used for synchronization.
    :param need_info: A flag indicating what information is needed during pangenome loading.
    """

    with ProcessPoolExecutor(max_workers=max_workers, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pan_name_to_path), unit='pangenome', disable=disable_bar) as progress:

            futures = []
            
            for pangenome_name, pangenome_path_info in pan_name_to_path.items():
                
                future = executor.submit(load_pangenome, pangenome_name, pangenome_path_info["path"], 
                                         pangenome_path_info["taxid"], need_info)

                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            
            for future in futures:
                results = future.result()
