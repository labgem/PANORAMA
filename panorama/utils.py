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
    global loading_lock
    if loading_lock is None:
        loading_lock = lock