#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from pathlib import Path
from typing import Dict, Tuple, Union
from multiprocessing import Manager, Lock
from concurrent.futures import ThreadPoolExecutor

# installed libraries
from tqdm import tqdm
import pandas as pd
import tables
from bokeh.plotting import save
from bokeh.models import ColumnDataSource, DataTable, TableColumn

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import init_lock, mkdir
from panorama.format.read_binaries import load_multiple_pangenomes


def check_former_info(args: argparse.Namespace, manager: Manager = None) -> Tuple[Dict[str, bool], Dict[str, dict]]:
    """ Check pangenomes path and needed information

    :param args: Argument from command line
    :param manager: Manager instance to use to create dictionaries

    :return: One dictionary to load pangenome info and one to store pangenome content
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

    return need_info, info_dict


def get_content_info(pangenome: Pangenome) -> Dict[str, Union[int, str]]:
    """ Get content information from pangenome file

    :param pangenome: Pangenome to get content information

    :return: Content information
    """
    h5f = tables.open_file(pangenome.file, "r")
    if "/info" in h5f:
        info_group = h5f.root.info
        req_info_list = ['numberOfGenes', 'numberOfOrganisms', 'numberOfClusters', 'numberOfEdges', 'numberOfCloud',
                         'numberOfPersistent', 'persistentStats', 'numberOfShell', 'shellStats', 'cloudStats',
                         'numberOfPartitions', 'numberOfSubpartitions']
        if all(info in info_group._v_attrs._f_list() for info in req_info_list):
            content_dic = {'Number of Genes': info_group._v_attrs['numberOfGenes'],
                           'Number of Genomes': info_group._v_attrs['numberOfOrganisms'],
                           'Number of Families': info_group._v_attrs['numberOfClusters'],
                           'Number of Edges': info_group._v_attrs['numberOfEdges'],
                           'Number of Persistent': info_group._v_attrs['numberOfPersistent'],
                           'Minimum of Persistent': info_group._v_attrs['persistentStats']['min'],
                           'Maximum of Persistent': info_group._v_attrs['persistentStats']['min'],
                           'Mean of Persistent': info_group._v_attrs['persistentStats']['mean'],
                           'Sd of Persistent': info_group._v_attrs['persistentStats']['sd'],
                           'Number of Shell': info_group._v_attrs['numberOfShell'],
                           'Minimum of Shell': info_group._v_attrs['shellStats']['min'],
                           'Maximum of Shell': info_group._v_attrs['shellStats']['min'],
                           'Mean of Shell': info_group._v_attrs['shellStats']['mean'],
                           'Sd of Shell': info_group._v_attrs['shellStats']['sd'],
                           'Number of Cloud': info_group._v_attrs['numberOfCloud'],
                           'Minimum of Cloud': info_group._v_attrs['cloudStats']['min'],
                           'Maximum of Cloud': info_group._v_attrs['cloudStats']['min'],
                           'Mean of Cloud': info_group._v_attrs['cloudStats']['mean'],
                           'Sd of Cloud': info_group._v_attrs['cloudStats']['sd']
                           }
        else:
            miss_list = []
            for info in req_info_list:
                if info not in info_group._v_attrs._f_list():
                    miss_list.append(info)
            raise Exception(f"No information about {','.join(miss_list)} find in {pangenome.file}"
                            f"Please use ppanggolin module -p {pangenome.file} to compute module and "
                            f"ppanggolin metrics --info_modules {pangenome.file} to compute information about module")
        if 'genome_fluidity' in info_group._v_attrs._f_list():
            for part, value in info_group._v_attrs['genomes_fluidity'].items():
                content_dic[f'Genome Fluidity {part}'] = value
        if 'family_fluidity' in info_group._v_attrs._f_list():
            for part, value in info_group._v_attrs['genomes_fluidity'].items():
                content_dic[f'Family Fluidity {part}'] = value
        if 'numberOfRGP' in info_group._v_attrs._f_list():
            content_dic["Number of RGPs"] = info_group._v_attrs['numberOfRGP']
        if 'numberOfSpots' in info_group._v_attrs._f_list():
            content_dic["Number of Spots"] = info_group._v_attrs['numberOfSpots']
        if 'numberOfModules' in info_group._v_attrs._f_list():
            content_dic["Number of Modules"] = info_group._v_attrs['numberOfModules']
            content_dic["Families in Modules"] = info_group._v_attrs['numberOfFamiliesInModules']
        h5f.close()
        return content_dic
    else:
        raise Exception(f"No information find in {pangenome.file}")


def get_module_info(pangenome: Pangenome) -> Dict[str, int]:
    """ Get module information from pangenome file

    :param pangenome: Pangenome with module computed

    :return: Dictionary with module information
    """
    h5f = tables.open_file(pangenome.file, "r")
    if "/info" in h5f:
        info_group = h5f.root.info
        module_info_list = ['numberOfModules', 'numberOfFamiliesInModules', 'PersistentSpecInModules',
                            'ShellSpecInModules', 'CloudSpecInModules', 'StatOfFamiliesInModules']
        if all(info in info_group._v_attrs._f_list() for info in module_info_list):
            module_dic = {'Number of Module': info_group._v_attrs['numberOfModules'],
                          'Number of families in modules': info_group._v_attrs['numberOfFamiliesInModules'],
                          "Percent of persistent families": info_group._v_attrs['PersistentSpecInModules']['percent'],
                          "Percent of shell families": info_group._v_attrs['ShellSpecInModules']['percent'],
                          "Percent of cloud families": info_group._v_attrs['CloudSpecInModules']['percent'],
                          "Minimum Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['min'],
                          "Maximum Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['max'],
                          "Sd Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['sd'],
                          "Mean Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['mean'],
                          }
            h5f.close()
            return module_dic
        else:
            raise Exception(f"No information about modules find in {pangenome.file} "
                            f"Please use ppanggolin module -p {pangenome.file} to compute module and "
                            f"ppanggolin metrics --info_modules {pangenome.file} to compute information about module")
    else:
        raise Exception(f"No information find in {pangenome.file}")


def get_info_pangenome(pangenome: Pangenome, info_dict: Dict[str, dict]):
    """ Compute in multiprocessing the genome fluidity of multiple pangenome

    :param pangenome: Pangenome
    :param info_dict: Dictionary that will contain information

    :return: The computed pangenome
    """
    if "content" in info_dict:
        info_dict["content"][pangenome.name] = get_content_info(pangenome)
    if "modules" in info_dict:
        info_dict["modules"][pangenome.name] = get_module_info(pangenome)


def get_info(pangenomes: Pangenomes, info_dict: Dict[str, dict], threads: int = 1,
             lock: Lock = None, disable_bar: bool = False) -> None:
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

    :param info_dict: Dictionary with information readable
    :param output: Path to output directory
    """
    def export_dict_to_tsv_and_html(info: dict, out_type: str):
        df = pd.DataFrame.from_dict(info, orient='index').T
        df.to_csv(path_or_buf=output / f'{out_type}_info.tsv', sep='\t')
        df = df.reset_index().rename(columns={"index": "Information"})
        source = ColumnDataSource(df)
        columns = [TableColumn(field=col, title=col) for col in df.columns]
        dt = DataTable(source=source, columns=columns, index_position=None,
                       sizing_mode='stretch_both')
        save(dt, filename=output / f'{out_type}_info.html', title=f"Pangenome {out_type} Information",
             resources="cdn")

    if "content" in info_dict:
        logging.debug("Write content info...")
        export_dict_to_tsv_and_html(info_dict["content"], "content")
    if "modules" in info_dict:
        logging.debug("Write modules info...")
        export_dict_to_tsv_and_html(info_dict["modules"], "module")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    logging.debug("launch info command")
    manager = Manager()
    lock = manager.Lock()
    outdir = mkdir(args.output, args.force)
    need_info, info_dict = check_former_info(args, manager)
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
