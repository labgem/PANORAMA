#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from tqdm import tqdm
from typing import List
from pathlib import Path

# installed libraries
import tables
from ppanggolin.formats import read_chunks
from ppanggolin.formats import check_pangenome_info as check_pp
from ppanggolin.formats import get_status as super_get_status

# local libraries
from panorama.system import System
from panorama.models import Models
from panorama.pangenomes import Pangenome


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



def check_pangenome_info(pangenome: Pangenome, sources: List[str] = None, source: str = None,
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
             need_gene_sequences=need_gene_sequences, need_modules=need_modules, disable_bar=disable_bar)

    if need_metadata:
        assert not all(x is not None for x in [source, sources])
        assert not all(x is None for x in [source, sources])
        if source is not None:
            check_pp(pangenome=pangenome, need_metadata=need_metadata, metatype=kwargs["metatype"], source=source,
                     disable_bar=disable_bar)
        elif sources is not None:
            for q_source in sources:
                check_pp(pangenome=pangenome, need_metadata=need_metadata, metatype=kwargs["metatype"], source=q_source,
                         disable_bar=disable_bar)

    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")

    if need_systems:
        assert models is not None and sources is not None
        read_systems(pangenome, h5f, models, sources, disable_bar)

    h5f.close()
