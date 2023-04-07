#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from typing import List
from pathlib import Path

# installed libraries
import tables
from tqdm import tqdm
from ppanggolin.formats import read_chunks
from ppanggolin.formats import check_pangenome_info as check_pp
from ppanggolin.formats import get_status as super_get_status

# local libraries
from panorama.annotation import Annotation
from panorama.models import Models
from panorama.system import System
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


def read_gene_families_annotations_by_source(pangenome: Pangenome, source_table, disable_bar: bool = False):
    for row in tqdm(read_chunks(source_table), total=source_table.nrows, unit="annotation", disable=disable_bar):
        gf = pangenome.get_gene_family(row["geneFam"].decode())
        annotation = Annotation(source=source_table.name, name=row['name'].decode(),
                                accession=row["accession"].decode(),
                                secondary_names=row["secondary_names"].decode(), score=row["score"],
                                description=row["description"].decode(), e_val=row["e_val"], bias=row["bias"])
        gf.add_annotation(source=source_table.name, annotation=annotation)


def read_gene_families_annotations(pangenome: Pangenome, h5f: tables.File, sources: List[str] = None,
                                   disable_bar: bool = False):
    """Read information about gene families in pangenome hdf5 file to add in pangenome object

    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    annotations_group = h5f.root.geneFamiliesAnnot
    if sources is None:
        sources = pangenome.status["annotations_sources"]

    for source in sources:
        annotations_table = h5f.get_node(annotations_group, source)
        logging.getLogger().info(f"Read annotations from {source}...")
        read_gene_families_annotations_by_source(pangenome, annotations_table, disable_bar)
        logging.getLogger().debug(f"{source} has been read")
    pangenome.status["annotations"] = "Loaded"


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
    logging.getLogger().info(f"Add system from {system_table.name} to pangenome...")
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
        logging.getLogger().info(f"Read system from {source}...")
        read_systems_by_source(pangenome, systems_table, models, disable_bar)
        logging.getLogger().debug(f"{source} has been read and added")
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

    check_pp(pangenome=pangenome, disable_bar=disable_bar, **kwargs)

    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")

    if need_systems:
        assert models is not None and sources is not None
        read_systems(pangenome, h5f, models, sources, disable_bar)

    h5f.close()
