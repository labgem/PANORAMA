#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from pathlib import Path
from typing import Dict

from tqdm import tqdm

# installed libraries
from ppanggolin.formats.readBinaries import check_pangenome_info
from multiprocessing import get_context, Manager

# local libraries
from panorama.pangenomes import Pangenome
from panorama.utils import check_tsv_sanity, mkdir
from panorama.info.module import get_module_info, export_modules
from panorama.info.content import get_content_info, export_content

need_info = {"need_annotations": False,
             "need_families": False,
             "need_graph": False,
             "need_partitions": False,
             "need_rgp": False,
             "need_spots": False,
             "need_gene_sequences": False,
             "need_modules": False,
             "need_content": False}


def check_former_info(args: argparse.Namespace) -> (Dict[str, Path], Dict[str, dict]):
    """ Check pangenomes path and needed information

    :param args: Argument from command line

    :return:
    """
    global need_info
    if not any(arg for arg in [args.modules, args.content]):
        raise argparse.ArgumentError(argument=None, message="You did not indicate which information you want.")
    manager = Manager()
    manager_dict = {}
    if args.modules is not None:
        need_info['need_families'] = True
        need_info['need_modules'] = True
        manager_dict["modules"] = manager.dict()
    if args.content:
        need_info['need_content'] = True
        manager_dict["content"] = manager.dict()

    return check_tsv_sanity(tsv_path=args.pangenomes), manager_dict


def mp_info(pan_name: str, pan_file: Path, manager_dict: dict, disable_bar: bool = False) -> Pangenome:
    """ Compute in multiprocessing the genome fluidity of multiple pangenome

    :param pan_name: Pangenome name
    :param pan_file: Path to pangenome
    :param manager_dict: Dictionary to save pangenome information in multiprocessing
    :param disable_bar: Allow to disable progress bar

    :return: The computed pangenome
    """
    global need_info
    pangenome = Pangenome(pan_name)
    pangenome.add_file(pan_file.absolute().as_posix())
    check_pangenome_info(pangenome=pangenome, **{k: v for k, v in need_info.items() if k != "need_content"},
                         disable_bar=disable_bar)
    if need_info['need_content']:
        manager_dict["content"][pangenome.name] = get_content_info(pangenome)
    if need_info['need_modules']:
        manager_dict["modules"][pangenome.name] = get_module_info(pangenome)
    return pangenome


def run_info(args: list) -> Pangenome:
    """Allow to get information with multiprocessing

    :param args: pan_name: str, pan_file: Path, manager_dict: dict, disable_bar: bool = False

    :return: The computed pangenome
    """
    return mp_info(*args)


def export_info(manager_dict: dict, output: Path, cpu: int = 1):
    """ Export information to HTML file

    :param manager_dict: Dictionary with information readable in multiprocessing
    :param output: Path to output directory
    :param cpu: Number of available CPU
    """
    processes = []
    with get_context('fork').Pool(processes=cpu) as p:
        if need_info['need_modules']:
            processes.append(p.apply_async(func=export_modules, args=(manager_dict["modules"], output)))
        if need_info['need_content']:
            processes.append(p.apply_async(func=export_content, args=(manager_dict["content"], output)))
        for process in processes:
            process.get()


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    global need_info
    logging.getLogger().debug("launch info command")
    outdir = mkdir(args.output, args.force)
    path_2_pang, manager_dict = check_former_info(args)
    mp_args = [(name, pan, manager_dict, args.disable_prog_bar) for name, pan in path_2_pang.items()]
    with get_context('fork').Pool(args.cpu) as p:
        for pangenome in tqdm(p.imap_unordered(run_info, mp_args), unit='pangenome',
                              total=len(path_2_pang), disable=args.disable_prog_bar):
            logging.getLogger().debug(f"{pangenome.name} Done")
    export_info(manager_dict, outdir, args.cpu)

    logging.getLogger().info("Done")


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
    required.add_argument('-o', '--output', required=True, type=str, nargs='?')
    onereq = parser.add_argument_group(title="Input file", description="One of the following argument is required :")
    onereq.add_argument("--content", required=False, action="store_true", default=False,
                        help="Create a detailed information TSV file about pangenomes content")
    onereq.add_argument('--modules', required=False, action="store_true", default=False,
                        help="Create a detailed information TSV file about pangenomes modules")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")


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
