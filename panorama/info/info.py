#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from pathlib import Path

from tqdm import tqdm

# installed libraries
from ppanggolin.formats.readBinaries import check_pangenome_info
from multiprocessing import get_context, Manager

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import check_tsv_sanity, mkdir
from panorama.info.module import get_module_info, export_modules

# def mp_info(pangenome_file: str):
#     pass

modules_dic = {}
need_info = {}


def check_info(args) -> dict:
    check_dict = {"need_annotations": False,
                  "need_families": False,
                  "need_graph": False,
                  "need_partitions": False,
                  "need_rgp": False,
                  "need_spots": False,
                  "need_gene_sequences": False,
                  "need_modules": False
                  }
    if args.modules is not None:
        check_dict['need_families'] = True
        check_dict['need_modules'] = True

    return check_dict


def check_former_info(args: argparse.Namespace):
    """ Check pangenomes path and needed information

    :param args:

    :return:
    """
    if not any(arg for arg in [args.modules, args.content]):
        raise Exception("You did not indicate which information you want.")
    return check_info(args), check_tsv_sanity(tsv_path=args.pangenomes)


def mp_info(pan_name: str, pan_file: Path) -> tuple:
    """ Compute in multiprocessing the genome fluidity of multiple pangenome

    :param pan_name: Pangenome name
    :param pan_file: Path to pangenome

    :return: Name of the computed pangenome
    """
    # TODO Allow to show progress bar of each process
    global need_info
    modules = None
    pangenome = Pangenome(pan_name)
    print(pangenome.name)
    pangenome.add_file(pan_file.absolute().as_posix())
    check_pangenome_info(pangenome=pangenome, **need_info)
    if need_info['need_modules']:
        modules = get_module_info(pangenome)
    return pangenome, modules


def run_info(args: list):
    return mp_info(*args)

def export_info(output: Path):
    if need_info['need_modules']:
        export_modules(modules_dic, output)


def launch(args: argparse.Namespace):
    """
    Command launcher
    :param args: All arguments provide by user
    """
    global need_info
    global modules_dic
    logging.getLogger().debug("launch info command")
    outdir = mkdir(args.output, args.force)
    need_info, path_2_pang = check_former_info(args)
    pangenomes = Pangenomes()
    mp_args = [(name, pan) for name, pan in path_2_pang.items()]
    with get_context('fork').Pool(args.cpu) as p:
        for pangenome, modules_info in tqdm(p.imap_unordered(run_info, mp_args),
                                            unit='pangenome', total=len(path_2_pang), disable=args.disable_prog_bar):
            pangenomes.add_pangenome(pangenome)
            modules_dic[pangenome.name] = modules_info
            logging.getLogger().debug(f"{pangenome} Done")
    export_info(outdir)

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
    import argparse

    from panorama.main import check_log

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

    launch(main_parser.parse_args())
