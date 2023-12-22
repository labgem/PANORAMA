#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from tqdm import tqdm
from typing import Any, Callable, Dict, List
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock

# installed libraries
import tables
from ppanggolin.formats import read_chunks
from ppanggolin.formats import check_pangenome_info as check_pp
from ppanggolin.formats import get_status as super_get_status

# local libraries
from panorama.system import System
from panorama.models import Models
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import check_tsv_sanity, init_lock


def get_status(pangenome, pangenome_file: str):
    """
        Checks which elements are already present in the file.
    """

    super_get_status(pangenome, pangenome_file)
    h5f = tables.open_file(pangenome_file, "r")
    status_group = h5f.root.status
    if hasattr(status_group._v_attrs, "systems") and status_group._v_attrs.systems:
        pangenome.status["systems"] = "inFile"
        pangenome.status["systems_sources"] = status_group._v_attrs.systems_sources
    h5f.close()


def read_systems_by_source(pangenome: Pangenome, system_table, models: Models, disable_bar: bool = False):
    """read sytem for one source

    :param pangenome: pangenome containing systems
    :param system_table: Table of systems corresponding to one source
    :param models: Models associated with systems
    :param disable_bar: Allow to disable progress bar
    """
    systems = {}
    for row in tqdm(read_chunks(system_table), total=system_table.nrows, unit="line", disable=disable_bar):
        curr_sys = systems.get(row["ID"].decode())
        model = models.get_model(row["name"].decode())
        if curr_sys is None:
            curr_sys = System(system_id=row["ID"].decode(), model=model,
                              source=system_table.name)
            systems[row["ID"].decode()] = curr_sys
        curr_sys.add_family(pangenome.get_gene_family(row["geneFam"].decode()))
    logging.info(f"Add system from {system_table.name} to pangenome...")
    for system in tqdm(systems.values(), unit="system", disable=disable_bar):
        pangenome.add_system(system)


def read_systems(pangenome: Pangenome, h5f: tables.File, models_path: List[Path], sources: List[str],
                 disable_bar: bool = False):
    """Read information about systems in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param models_path: list of path to models for each source
    :param sources: list of different source
    :param disable_bar: Disable the progress bar
    """
    systems_group = h5f.root.systems
    if sources is None:
        sources = pangenome.status["systems_sources"]

    for index, source in enumerate(sources):
        models = Models()
        models.read(models_path[index], disable_bar)
        systems_table = h5f.get_node(systems_group, source)
        logging.info(f"Read system from {source}...")
        read_systems_by_source(pangenome, systems_table, models, disable_bar)
        logging.debug(f"{source} has been read and added")
    pangenome.status["systems"] = "Loaded"


def check_pangenome_info(pangenome: Pangenome, sources: List[str] = None,
                         need_systems: bool = False, models: List[Path] = None,
                         disable_bar: bool = False, **kwargs):
    """
    Defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`

    :param pangenome: Pangenome object without some information
    :param need_systems:
    :param sources:
    :param models:
    :param disable_bar: Allow to disable the progress bar
    """

    need_annotations = kwargs["need_annotations"] if "need_annotations" in kwargs else False
    need_gene_sequences = kwargs["need_gene_sequences"] if "need_gene_sequences" in kwargs else False
    need_families = kwargs["need_families"] if "need_families" in kwargs else False
    need_partitions = kwargs["need_partitions"] if "need_partitions" in kwargs else False
    need_graph = kwargs["need_graph"] if "need_graph" in kwargs else False
    need_rgp = kwargs["need_rgp"] if "need_rgp" in kwargs else False
    need_spots = kwargs["need_spots"] if "need_spots" in kwargs else False
    need_modules = kwargs["need_modules"] if "need_modules" in kwargs else False
    need_metadata = kwargs["need_metadata"] if "need_metadata" in kwargs else False

    check_pp(pangenome, need_annotations=need_annotations, need_families=need_families, need_graph=need_graph,
             need_partitions=need_partitions, need_rgp=need_rgp, need_spots=need_spots,
             need_gene_sequences=need_gene_sequences, need_modules=need_modules,
             need_metadata=need_metadata, metatypes={"families"}, sources=sources, disable_bar=disable_bar)

    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")

    if need_systems:
        assert models is not None and sources is not None
        read_systems(pangenome, h5f, models, sources, disable_bar)

    h5f.close()


def load_pangenome(name: str, path: str, taxid: int, need_info: Dict[str, bool],
                   check_function: Callable[[Pangenome, Any], None] = None,
                   disable_bar: bool = False, **kwargs) -> Pangenome:
    """
    Load a pangenome from a given path and check the required information.

    This function loads a pangenome from the specified `path` and assigns it the provided `name` and `taxid`.
    The pangenome file is added to the pangenome object. The function then checks that the required information
    are present in the pangenome and if they are, it loads them.

    :param name: The name of the pangenome.
    :param path: The path to the pangenome file.
    :param taxid: The taxonomic ID associated with the pangenome.
    :param need_info: A dictionary containing information required to load in the Pangenome object.
    :param check_function:
    :param disable_bar:
    :return: The pangenome object with the loaded information.
    """
    pangenome = Pangenome(name=name, taxid=taxid)
    pangenome.add_file(path)

    if check_function is not None:
        try:
            check_function(pangenome, **kwargs)
        except Exception as error:
            logging.error(f"Pangenome {pangenome.name} reading return the below error")
            raise error
    check_pangenome_info(pangenome, disable_bar=disable_bar, **need_info)
    return pangenome


def load_pangenomes(pangenome_list: Path, need_info: Dict[str, bool],
                             check_function: Callable[[Pangenome, Any], None] = None,
                             max_workers: int = 1, lock: Lock = None,
                             disable_bar: bool = False, **kwargs) -> Pangenomes:
    """
    Load multiple pangenomes in parallel using a process pool executor.

    This function loads multiple pangenomes in parallel using a process pool executor. It takes a dictionary
    `pan_name_to_path` containing the mapping of pangenome names to their corresponding paths and other
    information. The pangenomes are loaded using the `load_pangenome` function. The loading progress is
    displayed using a tqdm progress bar.

    :param pangenome_list: Path to the pangenomes list files.
    :param need_info: A flag indicating what information is needed during pangenome loading.
    :param check_function: Function which check the pangenome before to load information
    :param max_workers: The maximum number of worker processes to use in the process pool executor.
    :param lock: A multiprocessing lock used for synchronization.

    :param disable_bar: A flag indicating whether to disable the tqdm progress bar.
    :return pangenomes: List of loaded pangenomes with required information
    """
    pangenomes = Pangenomes()
    pan_to_path = check_tsv_sanity(pangenome_list)
    with ThreadPoolExecutor(max_workers=max_workers, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pan_to_path), unit='pangenome', disable=disable_bar) as progress:
            futures = []

            for pangenome_name, pangenome_path_info in pan_to_path.items():
                future = executor.submit(load_pangenome, pangenome_name, pangenome_path_info["path"],
                                         pangenome_path_info["taxid"], need_info, check_function,
                                         disable_bar, **kwargs)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                with lock:
                    pangenomes.add_pangenome(future.result())
    return pangenomes
