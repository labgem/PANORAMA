#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import sys
import logging
from pathlib import Path
import csv
from typing import TextIO, Dict, Union
from multiprocessing import Manager, Lock
from importlib.metadata import distribution


# installed libraries

# local libraries


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


def add_common_arguments(subparser: argparse.ArgumentParser) -> None:
    """
    Add common argument to the input subparser.

    :param subparser: A subparser object from any subcommand.
    """
    common = subparser._action_groups.pop(1)  # get the 'optional arguments' action group.
    common.title = "Common arguments"
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
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

        if args.log != sys.stdout and not args.disable_prog_bar:  # if output is not to stdout we remove progress bars.
            args.disable_prog_bar = True
        str_format = "%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s"
        datefmt = '%Y-%m-%d %H:%M:%S'
        if args.log in [sys.stdout, sys.stderr]:
            # use stream
            logging.basicConfig(stream=args.log, level=level,
                                format=str_format,
                                datefmt=datefmt)
        else:
            # log is written in a files. basic condif uses filename
            logging.basicConfig(filename=args.log, level=level,
                                format=str_format,
                                datefmt=datefmt)
        logging.getLogger("PANORAMA").info("Command: " + " ".join([arg for arg in sys.argv]))
        logging.getLogger("PANORAMA").info(f"PPanGGOLiN version: {distribution('ppanggolin').version}")


# File managing system
def mkdir(output: Path, force: bool = False) -> Path:
    """Create a directory at the given path

    :param output: Path to output directory
    :param force: Pass exception if directory already exist

    :raise FileExistError: If the given directory already exist and no force is used
    :raise Exception: Handle all others exception

    :return: Path object to output directory
    """
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError:
        if not force:
            raise FileExistsError(f"{output} already exists."
                                  f"Use --force if you want to overwrite the files in the directory")
        else:
            logging.warning(f"{output.as_posix()} already exist and file could be overwrite by the new generated")
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
        raise Exception("An unexpected error happened. Please report to our github.")
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
                    abs_path = path_exist(tsv_path.parent / Path(line[1]))
                except FileNotFoundError as file_error:
                    raise FileNotFoundError(f"{file_error}")
                else:
                    pan_to_path[line[0]] = {"path": abs_path,
                                            "taxid": line[2] if len(line) > 2 else None}
            except Exception:
                raise Exception("Unexpected error")
            else:
                pan_to_path[line[0]] = {"path": abs_path,
                                        "taxid": line[2] if len(line) > 2 else None}
        p_file.close()
        return pan_to_path


def init_lock(lock: Lock = None):
    """
    Initialize the loading lock.

    This function initializes the `loading_lock` variable as a global variable,
    assigning it the value of the `lock` parameter.
    If the `loading_lock` is already initialized, the function does nothing.

    :param lock: The lock object to be assigned to `loading_lock`.
    """
    if lock is None:
        manager = Manager()
        return manager.Lock()
