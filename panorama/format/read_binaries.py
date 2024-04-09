#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging

from tqdm import tqdm
from typing import Callable, Dict, List
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock

# installed libraries
import tables
from ppanggolin.formats import (read_chunks, read_annotation, read_graph, read_rgp, read_modules,
                                read_gene_sequences, read_metadata, get_need_info)
from ppanggolin.formats import get_status as super_get_status
from ppanggolin.geneFamily import Gene

# local libraries
from panorama.systems.system import System
from panorama.systems.models import Models
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.region import Spot
from panorama.utils import check_tsv_sanity, init_lock


def get_status(pangenome, pangenome_file: Path):
    """
        Checks which elements are already present in the file.
    """

    super_get_status(pangenome, pangenome_file)
    h5f = tables.open_file(pangenome_file.absolute().as_posix(), "r")
    status_group = h5f.root.status
    if hasattr(status_group._v_attrs, "systems") and status_group._v_attrs.systems:
        pangenome.status["systems"] = "inFile"
        pangenome.status["systems_sources"] = status_group._v_attrs.systems_sources
    h5f.close()


def read_systems_by_source(pangenome: Pangenome, system_table, models: Models, disable_bar: bool = False):
    """read system for one source

    :param pangenome: pangenome containing systems
    :param system_table: Table of systems corresponding to one source
    :param models: Models associated with systems
    :param disable_bar: Allow to disable progress bar
    """
    systems = {}
    for row in tqdm(read_chunks(system_table), total=system_table.nrows, unit="line", disable=disable_bar):
        sys_id = row["ID"].decode()
        if sys_id not in systems:
            model = models.get_model(row["name"].decode())
            system = System(system_id=row["ID"].decode(), model=model,
                            source=system_table.name)
            systems[row["ID"].decode()] = system
        else:
            system = systems[sys_id]
        system.add_family(pangenome.get_gene_family(row["geneFam"].decode()))
    logging.getLogger("PANORAMA").info(f"Add system from {system_table.name} to pangenome...")
    logging.getLogger("PANORAMA").debug(f"Number of systems found: {len(systems)}")
    for system in tqdm(sorted(systems.values(), key=lambda x: len(x), reverse=True), unit="system",
                       disable=False if logging.getLogger().level == logging.DEBUG else True):
        pangenome.add_system(system)


def read_systems(pangenome: Pangenome, h5f: tables.File, models: List[Models], sources: List[str],
                 disable_bar: bool = False):
    """Read information about systems in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param models: list of models for each source
    :param sources: list of different source
    :param disable_bar: Disable the progress bar
    """
    systems_group = h5f.root.systems

    for index, source in enumerate(sources):
        systems_table = h5f.get_node(systems_group, source)
        logging.getLogger("PANORAMA").info(f"Read system from {source}...")
        read_systems_by_source(pangenome, systems_table, models[index], disable_bar)
        logging.getLogger("PANORAMA").debug(f"{source} has been read and added")
    pangenome.status["systems"] = "Loaded"


def read_gene_families_info(pangenome: Pangenome, h5f: tables.File, information: bool = False,
                            sequences: bool = False, disable_bar: bool = False):
    """
    Read information about gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.geneFamiliesInfo
    description = "Reading gene families "
    if information:
        description += f"information {'and sequences' if sequences else ''}"
    else:
        description += f"{'sequences' if sequences else ''}"

    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene family",
                    desc=description, disable=disable_bar):
        fam = pangenome.get_gene_family(row["name"].decode())
        if information:
            fam.partition = row["partition"].decode()
        if sequences:
            fam.add_sequence(row["protein"].decode())

    if information and h5f.root.status._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "Loaded"
    if sequences and h5f.root.status._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "Loaded"


def read_gene_families(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.geneFamilies

    link = True if pangenome.status["genomesAnnotated"] in ["Computed", "Loaded"] else False

    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="gene",
                    desc="Associate gene to gene families", disable=disable_bar):
        try:
            fam = pangenome.get_gene_family(name=row["geneFam"].decode())
        except KeyError:
            fam = GeneFamily(family_id=pangenome.max_fam_id, name=row["geneFam"].decode())
            pangenome.add_gene_family(fam)
        if link:  # linking if we have loaded the annotations
            gene_obj = pangenome.get_gene(row["gene"].decode())
        else:  # else, no
            gene_obj = Gene(row["gene"].decode())
        fam.add(gene_obj)
    pangenome.status["genesClustered"] = "Loaded"


def read_spots(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read hotspot in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without spot
    :param h5f: Pangenome HDF5 file with spot computed
    :param disable_bar: Disable the progress bar
    """
    table = h5f.root.spots
    spots = {}
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="spot", disable=disable_bar):
        curr_spot = spots.get(int(row["spot"]))
        if curr_spot is None:
            curr_spot = Spot(int(row["spot"]))
            spots[row["spot"]] = curr_spot
        region = pangenome.get_region(row["RGP"].decode())
        curr_spot.add(region)
        curr_spot.spot_2_families()
    for spot in spots.values():
        pangenome.add_spot(spot)
    pangenome.status["spots"] = "Loaded"


def read_pangenome(pangenome: Pangenome, annotation: bool = False, gene_families: bool = False, graph: bool = False,
                   rgp: bool = False, spots: bool = False, gene_sequences: bool = False, modules: bool = False,
                   metadata: bool = False, systems: bool = False, disable_bar: bool = False, **kwargs):
    """
    Reads a previously written pangenome, with all of its parts, depending on what is asked,
    with regard to what is filled in the 'status' field of the hdf5 file.

    :param pangenome: Pangenome object without some information
    :param annotation: get annotation
    :param gene_families: get gene families
    :param graph: get graph
    :param rgp: get RGP
    :param spots: get hotspot
    :param gene_sequences: get gene sequences
    :param modules: get modules
    :param metadata: get metadata
    :param systems: get systems
    :param disable_bar: Allow to disable the progress bar
    :param kwargs: others parameters to get attributes
    """
    if pangenome.file is None:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")

    h5f = tables.open_file(pangenome.file, "r")

    if annotation:  # I place annotation here, to link gene to gene families if organism are not loaded
        if h5f.root.status._v_attrs.genomesAnnotated:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome annotations...")
            read_annotation(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise ValueError(f"The pangenome in file '{pangenome.file}' has not been annotated, "
                             "or has been improperly filled")

    if gene_sequences:
        if h5f.root.status._v_attrs.geneSequences:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome gene dna sequences...")
            read_gene_sequences(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise ValueError(f"The pangenome in file '{pangenome.file}' does not have gene sequences, "
                             "or has been improperly filled")

    if gene_families:
        if h5f.root.status._v_attrs.genesClustered:
            logging.getLogger("PPanGGOLiN").info("Reading pangenome gene families...")
            read_gene_families(pangenome, h5f, disable_bar=disable_bar)
            if kwargs["gene_families_info"] or kwargs["gene_families_sequences"]:
                debug_msg = "Reading pangenome gene families "
                if kwargs["gene_families_info"]:
                    debug_msg += f"info{' and sequences...' if kwargs['gene_families_sequences'] else '...'}"
                elif kwargs["gene_families_sequences"]:
                    debug_msg += "sequences..."
                logging.getLogger("PPanGGOLiN").debug(debug_msg)
                read_gene_families_info(pangenome, h5f, kwargs["gene_families_info"],
                                        kwargs["gene_families_sequences"], disable_bar)
        else:
            raise ValueError(f"The pangenome in file '{pangenome.file}' does not have gene families, "
                             "or has been improperly filled")

    if graph:
        if h5f.root.status._v_attrs.NeighborsGraph:
            logging.getLogger("PPanGGOLiN").info("Reading the neighbors graph edges...")
            read_graph(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise AttributeError(f"The pangenome in file '{pangenome.file}' does not have graph information, "
                                 f"or has been improperly filled")

    if rgp:
        if h5f.root.status._v_attrs.predictedRGP:
            logging.getLogger("PPanGGOLiN").info("Reading the RGP...")
            read_rgp(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise AttributeError(f"The pangenome in file '{pangenome.file}' does not have RGP information, "
                                 f"or has been improperly filled")

    if spots:
        if h5f.root.status._v_attrs.spots:
            logging.getLogger("PPanGGOLiN").info("Reading the spots...")
            read_spots(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise AttributeError(f"The pangenome in file '{pangenome.file}' does not have spots information, "
                                 f"or has been improperly filled")

    if modules:
        if h5f.root.status._v_attrs.modules:
            logging.getLogger("PPanGGOLiN").info("Reading the modules...")
            read_modules(pangenome, h5f, disable_bar=disable_bar)
        else:
            raise AttributeError(f"The pangenome in file '{pangenome.file}' does not have modules information, "
                                 f"or has been improperly filled")

    if metadata:
        for metatype in kwargs["metatypes"]:

            if h5f.root.status._v_attrs.metadata:
                metastatus = h5f.root.status._f_get_child("metastatus")
                metasources = h5f.root.status._f_get_child("metasources")

                metatype_sources = set(metasources._v_attrs[metatype]) & set(kwargs["meta_sources"])
                if metastatus._v_attrs[metatype] and len(metatype_sources) > 0:
                    logging.getLogger("PPanGGOLiN").info(
                        f"Reading the {metatype} metadata from sources {metatype_sources}...")
                    read_metadata(pangenome, h5f, metatype, metatype_sources, disable_bar=disable_bar)
            else:
                raise KeyError(
                    f"The pangenome in file '{pangenome.file}' does not have metadata associated to {metatype}, ")

    if systems:
        read_systems(pangenome, h5f, kwargs["models"], kwargs["systems_sources"], disable_bar)
    h5f.close()


def check_pangenome_info(pangenome, need_families_info: bool = False, need_families_sequences: bool = False,
                         need_systems: bool = False, models: List[Models] = None,
                         systems_sources: List[str] = None, disable_bar: bool = False, **kwargs):
    """
    Defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`

    :param pangenome: Pangenome object without some information
    :param disable_bar: Allow to disable the progress bar
    """

    need_info = get_need_info(pangenome, **kwargs)
    need_info["meta_sources"] = need_info.pop("sources")

    if need_families_info or need_families_sequences:
        if kwargs.get('gene_families'):
            raise AssertionError("gene families need to be loaded to load either information or sequences.")
    need_info["gene_families_info"] = need_families_info
    need_info["gene_families_sequences"] = need_families_sequences

    if need_systems:
        assert models is not None and systems_sources is not None
        need_info["systems"] = True
        need_info["models"] = models
        need_info["systems_sources"] = systems_sources

    logging.getLogger("PANORAMA").debug(f"need_info: {need_info}")

    if any(need_info.values()):
        # if no flag is true, then nothing is needed.
        read_pangenome(pangenome, disable_bar=disable_bar, **need_info)


def load_pangenome(name: str, path: Path, taxid: int, need_info: Dict[str, bool],
                   check_function: Callable[[Pangenome, ...], None] = None,
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
            logging.getLogger("PANORAMA").error(f"Pangenome {pangenome.name} reading return the below error")
            raise error
    check_pangenome_info(pangenome, disable_bar=disable_bar, **need_info)
    return pangenome


def load_pangenomes(pangenome_list: Path, need_info: dict, check_function: callable = None,
                    max_workers: int = 1, lock: Lock = None, disable_bar: bool = False, **kwargs: object) -> Pangenomes:
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
                    pangenomes.add(future.result())
    return pangenomes
