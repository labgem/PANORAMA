#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from pathlib import Path
from typing import Dict
from multiprocessing import Manager, Lock
from concurrent.futures import ThreadPoolExecutor

# installed libraries
from tqdm import tqdm

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import init_lock, mkdir
from panorama.format.read_binaries import load_multiple_pangenomes
from panorama.info.module import get_module_info, export_modules
from panorama.info.content import get_content_info, export_content


def check_former_info(args: argparse.Namespace, manager: Manager = None) -> (Dict[str, Path], Dict[str, dict]):
    """ Check pangenomes path and needed information

    :param args: Argument from command line

    :return:
    """
    manager = Manager() if manager is None else manager
    need_info = manager.dict()
    info_dict = manager.dict()
    if not any(arg for arg in [args.modules, args.content]):
        raise argparse.ArgumentError(argument=None, message="You did not indicate which information you want.")
    if args.modules:
        need_info['need_families'] = True
        need_info['need_modules'] = True
        info_dict["modules"] = manager.dict()
    if args.content:
        need_info['need_content'] = True
        info_dict["content"] = manager.dict()

    return info_dict, need_info


def get_info_pangenome(pangenome: Pangenome, info_dict: dict, manager: Manager = None):
    """ Compute in multiprocessing the genome fluidity of multiple pangenome

    :param pangenome: Pangenome


    :return: The computed pangenome
    """
    manager = Manager() if manager is None else manager
    if "content" in info_dict:
        info_dict["content"][pangenome.name] = get_content_info(pangenome)
    if "modules" in info_dict:
        info_dict["modules"][pangenome.name] = get_module_info(pangenome)


def get_info(pangenomes: Pangenomes, info_dict: dict, threads: int = 1, lock: Lock = None, disable_bar: bool = False) -> Pangenome:
    """Allow to get information with multiprocessing

    :param args: pan_name: str, pan_file: Path, manager_dict: dict, disable_bar: bool = False

    :return: The computed pangenome
    """
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            logging.info("Get pangenome information...")
            for pangenome in pangenomes:
                future = executor.submit(get_info_pangenome, pangenome, info_dict)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def export_info(info_dict: dict, output: Path):
    """ Export information to HTML file

    :param manager_dict: Dictionary with information readable in multiprocessing
    :param output: Path to output directory
    :param cpu: Number of available CPU
    """
    if "content" in info_dict:
        logging.debug("Write content info...")
        export_content(info_dict["content"], output)
    if "modules" in info_dict:
        logging.debug("Write modules info...")
        export_modules(info_dict["modules"], output)


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    logging.debug("launch info command")
    manager = Manager()
    lock = manager.Lock()
    outdir = mkdir(args.output, args.force)
    info_dict, need_info = check_former_info(args, manager)
    pangenomes = load_multiple_pangenomes(pangenome_list=args.pangenomes, need_info=need_info, lock=lock,
                                          max_workers=args.threads, force=args.force, disable_bar=args.disable_prog_bar)
    get_info(pangenomes=pangenomes, info_dict=info_dict, threads=args.threads,
             lock=lock, disable_bar=args.disable_prog_bar)
    export_info(info_dict, outdir)

    logging.info("Done")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_info(parser)
    return parser


def parser_info(parser):
    """
    Parser for specific argument of info command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help="A list of pangenome .h5 files")
    required.add_argument('-o', '--output', required=True, type=Path, nargs='?')
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument("--content", required=False, action="store_true",
                        help="Create a detailed information TSV file about pangenomes content")
    onereq.add_argument('--modules', required=False, action="store_true",
                        help="Create a detailed information TSV file about pangenomes modules")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
        formatter_class=argparse.RawTextHelpFormatter)

    parser_info(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
