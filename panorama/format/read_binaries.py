#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to read and load pangenome data from HDF5 files.
"""

# default libraries
import logging

from tqdm import tqdm
from typing import Callable, Dict, List
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock

# installed libraries
import tables
from ppanggolin.formats import (
    read_chunks, read_annotation, read_graph, read_rgp,
    read_gene_sequences, read_metadata, get_need_info
)
from ppanggolin.formats import get_status as super_get_status
from ppanggolin.geneFamily import Gene

# local libraries
from panorama.systems.system import System, SystemUnit
from panorama.systems.models import Models
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.region import Spot, Module
from panorama.utils import check_tsv_sanity, init_lock


def get_status(pangenome, pangenome_file: Path):
    """
    Check which elements are already present in the file and update the pangenome status.

    Args:
        pangenome (Pangenome): The pangenome object to update.
        pangenome_file (Path): The path to the pangenome HDF5 file.
    """
    super_get_status(pangenome, pangenome_file)
    h5f = tables.open_file(pangenome_file.absolute().as_posix(), "r")
    status_group = h5f.root.status
    if hasattr(status_group._v_attrs, "systems") and status_group._v_attrs.systems:
        pangenome.status["systems"] = "inFile"
        if hasattr(status_group._v_attrs, 'systems_sources'):
            pangenome.status["systems_sources"] = status_group._v_attrs.systems_sources
        else:
            pangenome.status["systems_sources"] = set()
    h5f.close()


def read_systems_by_source(pangenome: Pangenome, source_group: tables.Group, models: Models,
                           read_canonical: bool = True, disable_bar: bool = False):
    """
    Read systems for one source and add them to the pangenome.

    Args:
        pangenome (Pangenome): Pangenome containing systems.
        source_group: Source group with 3 tables to read systems.
        models (Models): Models associated with systems.
        read_canonical (bool, optional): Read canonical systems (default True)
        disable_bar (bool, optional): Whether to disable the progress bar (default False).
    """

    def read_system_unit(unit_row: tables.Table.row, units_dict: Dict[str, SystemUnit]):
        """
        Global function to read a line in unit table

        Args:
            unit_row: the row with unit information to read
            units_dict: Dictionary of read units, with identifier as key and unit as value
        """
        unit_id = unit_row["ID"]
        if unit_id not in units_dict:
            fu_name = unit_row["name"].decode()
            model = models.get_model(fu2model[fu_name])
            unit = SystemUnit(functional_unit=model.get(fu_name), source=source)
            unit.ID = unit_id
            units_dict[unit_id] = unit
        else:
            unit = units_dict[unit_id]
        unit.add_family(pangenome.get_gene_family(unit_row["geneFam"].decode()),
                        unit_row["metadata_source"].decode(), int(unit_row["metadata_id"]))

    def read_system(sys_row: tables.Table.row, sys_dict: Dict[str, System], unit_dict: Dict[str, SystemUnit]) -> System:
        """
        Global function to read a line in system table

        Args:
            sys_row: the row with system information to read
            sys_dict: Dictionary of read systems, with identifier as key and system as value
        """
        sys_id = sys_row["ID"].decode()
        if sys_id not in sys_dict:
            model = models.get_model(sys_row["name"].decode())
            sys = System(system_id=sys_id, model=model, source=source)
            sys_dict[sys_id] = sys
        else:
            sys = sys_dict[sys_id]
        sys.add_unit(unit_dict[sys_row["unit"]])
        return sys

    source = source_group._v_name
    fu2model = {fu.name: model.name for model in models for fu in model.func_units}
    system_table = source_group.systems
    unit_table = source_group.units
    systems = {}
    with tqdm(total=system_table.nrows + unit_table.nrows, unit="line",
              desc="Read pangenome systems", disable=disable_bar) as progress:
        units = {}
        for row in read_chunks(unit_table):
            read_system_unit(row, units)
            progress.update()
        for row in read_chunks(system_table):
            read_system(row, systems, units)
    logging.getLogger("PANORAMA").debug(f"Number of systems read: {len(systems)}")
    if read_canonical:
        canonic = {}
        canon_table = source_group.canonic
        canon_unit_table = source_group.canonic_units
        sys2canonical_table = source_group.system_to_canonical

        with tqdm(total=canon_table.nrows + canon_unit_table.nrows + sys2canonical_table.nrows,
                  unit="line", desc="Read pangenome canonical systems", disable=disable_bar) as progress:
            canon_units = {}
            for row in read_chunks(canon_unit_table):
                read_system_unit(row, canon_units)
                progress.update()

            for row in read_chunks(canon_table):
                read_system(row, canonic, canon_units)
                progress.update()

            for row in read_chunks(sys2canonical_table):
                systems[row["system"].decode()].add_canonical(canonic[row["canonic"].decode()])
                progress.update()
            logging.getLogger("PANORAMA").debug(f"Number of canonical system: {len(canonic)}")

    logging.getLogger("PANORAMA").info(f"Add system from {source} to pangenome...")
    for system in tqdm(sorted(systems.values(), key=lambda x: (len(x.model.canonical), -len(x), -x.number_of_families)),
                       disable=False if logging.getLogger().level == logging.DEBUG else True,
                       desc="Add system to pangenome"):
        pangenome.add_system(system)
    logging.getLogger("PANORAMA").info(f"Add {pangenome.number_of_systems(source=source, with_canonical=False)} "
                                       f"systems from {source} to pangenome...")


def read_systems(pangenome: Pangenome, h5f: tables.File, models: List[Models], sources: List[str],
                 disable_bar: bool = False):
    """
    Read information about systems in the pangenome HDF5 file and add them to the pangenome object.

    Args:
        pangenome (Pangenome): Pangenome object.
        h5f (tables.File): Pangenome HDF5 file with gene families information.
        models (List[Models]): List of models for each source.
        sources (List[str]): List of different sources.
        disable_bar (bool): Whether to disable the progress bar.
    """
    systems_group = h5f.root.systems
    metadata_sources = set()
    for index, source in enumerate(sources):
        source_group = h5f.get_node(systems_group, source)
        metadata_sources |= source_group._v_attrs.metadata_sources
        logging.getLogger("PANORAMA").info(f"Read system from {source}...")
        read_systems_by_source(pangenome, source_group, models[index], disable_bar=disable_bar)
        logging.getLogger("PANORAMA").debug(f"{source} has been read and added")
    pangenome.status["systems"] = "Loaded"
    return metadata_sources


def read_gene_families_info(pangenome: Pangenome, h5f: tables.File, information: bool = False,
                            sequences: bool = False, disable_bar: bool = False):
    """
    Read information about gene families in the pangenome HDF5 file and add them to the pangenome object.

    Args:
        pangenome (Pangenome): Pangenome object.
        h5f (tables.File): Pangenome HDF5 file with gene families information.
        information (bool): Whether to read information.
        sequences (bool): Whether to read sequences.
        disable_bar (bool): Whether to disable the progress bar.
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
    Read gene families in the pangenome HDF5 file and add them to the pangenome object.

    Args:
        pangenome (Pangenome): Pangenome object.
        h5f (tables.File): Pangenome HDF5 file with gene families information.
        disable_bar (bool): Whether to disable the progress bar.
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
    Read hotspots in the pangenome HDF5 file and add them to the pangenome object.

    Args:
        pangenome (Pangenome): Pangenome object.
        h5f (tables.File): Pangenome HDF5 file with spots computed.
        disable_bar (bool): Whether to disable the progress bar.
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


def read_modules(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read modules in the pangenome HDF5 file and add them to the pangenome object.

    Args:
        pangenome (Pangenome): Pangenome object.
        h5f (tables.File): Pangenome HDF5 file with modules computed.
        disable_bar (bool): Whether to disable the progress bar.

    Raises:
        Exception: If gene families have not been loaded.
    """
    if pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
        raise Exception("It's not possible to read the modules if the gene families have not been loaded.")
    table = h5f.root.modules
    modules = {}  # id2mod
    for row in tqdm(read_chunks(table, chunk=20000), total=table.nrows, unit="module", disable=disable_bar):
        curr_module = modules.get(int(row["module"]))
        if curr_module is None:
            curr_module = Module(int(row["module"]))
            modules[row["module"]] = curr_module
        family = pangenome.get_gene_family(row["geneFam"].decode())
        curr_module.add(family)
    for module in modules.values():
        pangenome.add_module(module)
    pangenome.status["modules"] = "Loaded"


def read_pangenome(pangenome: Pangenome, annotation: bool = False, gene_families: bool = False, graph: bool = False,
                   rgp: bool = False, spots: bool = False, gene_sequences: bool = False, modules: bool = False,
                   metadata: bool = False, systems: bool = False, disable_bar: bool = False, **kwargs):
    """
    Read a previously written pangenome with all of its parts, depending on what is asked,
    and what is filled in the 'status' field of the HDF5 file.

    Args:
        pangenome (Pangenome): Pangenome object without some information.
        annotation (bool): Whether to read the annotation.
        gene_families (bool): Whether to read gene families.
        graph (bool): Whether to read the graph.
        rgp (bool): Whether to read RGP.
        spots (bool): Whether to read hotspots.
        gene_sequences (bool): Whether to read gene sequences.
        modules (bool): Whether to read modules.
        metadata (bool): Whether to read metadata.
        systems (bool): Whether to read systems.
        disable_bar (bool): Whether to disable the progress bar.
        **kwargs: Additional parameters to get attributes.

    Raises:
        FileNotFoundError: If the provided pangenome does not have an associated .h5 file.
        ValueError: If the required annotation, gene families, gene sequences, or RGP information is not present in the file.
        AttributeError: If the required graph, spots, or modules information is not present in the file.
        KeyError: If the required metadata information is not present in the file.
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

    if systems:
        metadata = True
        metadata_sources = read_systems(pangenome, h5f, kwargs["models"], kwargs["systems_sources"], disable_bar)
        if "meta_sources" in kwargs:
            kwargs["metatypes"].add("families")
            kwargs["meta_sources"] |= metadata_sources
        else:
            kwargs["metatypes"] = {"families"}
            kwargs["meta_sources"] = metadata_sources

    if metadata:
        for metatype in kwargs["metatypes"]:
            if h5f.root.status._v_attrs.metadata:
                metastatus = h5f.root.status._f_get_child("metastatus")
                metasources = h5f.root.status._f_get_child("metasources")

                metatype_sources = set(metasources._v_attrs[metatype]) & kwargs["sources"]
                if "meta_sources" in kwargs:
                    metatype_sources &= kwargs["meta_sources"]
                if metastatus._v_attrs[metatype] and len(metatype_sources) > 0:
                    logging.getLogger("PPanGGOLiN").info(
                        f"Reading the {metatype} metadata from sources {metatype_sources}...")
                    read_metadata(pangenome, h5f, metatype, metatype_sources, disable_bar=disable_bar)
            else:
                raise KeyError(
                    f"The pangenome in file '{pangenome.file}' does not have metadata associated to {metatype}, ")
    h5f.close()


def check_pangenome_info(pangenome, need_families_info: bool = False, need_families_sequences: bool = False,
                         need_systems: bool = False, models: List[Models] = None,
                         systems_sources: List[str] = None, disable_bar: bool = False, **kwargs):
    """
    Defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`.

    Args:
        pangenome (Pangenome): Pangenome object without some information.
        need_families_info (bool): Whether gene families info is needed.
        need_families_sequences (bool): Whether gene families sequences are needed.
        need_systems (bool): Whether systems are needed.
        models (List[Models]): List of models.
        systems_sources (List[str]): List of systems sources.
        disable_bar (bool): Whether to disable the progress bar.
        **kwargs: Additional parameters to get attributes.

    Raises:
        AssertionError: If gene families need to be loaded to load either information or sequences.
    """
    need_info = get_need_info(pangenome, **kwargs)

    if need_families_info or need_families_sequences:
        if kwargs.get('gene_families'):
            raise AssertionError("Gene families need to be loaded to load either information or sequences.")
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

    Args:
        name (str): The name of the pangenome.
        path (Path): The path to the pangenome file.
        taxid (int): The taxonomic ID associated with the pangenome.
        need_info (Dict[str, bool]): A dictionary containing information required to load in the Pangenome object.
        check_function (Callable[[Pangenome, ...], None], optional): Function to check the pangenome before loading information.
        disable_bar (bool): Whether to disable the progress bar.
        **kwargs: Additional parameters to get attributes.

    Returns:
        Pangenome: The pangenome object with the loaded information.

    Raises:
        Exception: If an error occurs during the pangenome check.
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

    Args:
        pangenome_list (Path): Path to the pangenomes list files.
        need_info (dict): A flag indicating what information is needed during pangenome loading.
        check_function (callable, optional): Function to check the pangenome before loading information.
        max_workers (int): The maximum number of worker processes to use in the process pool executor.
        lock (Lock, optional): A multiprocessing lock used for synchronization.
        disable_bar (bool): Whether to disable the tqdm progress bar.
        **kwargs: Additional parameters to get attributes.

    Returns:
        Pangenomes: List of loaded pangenomes with required information.
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
