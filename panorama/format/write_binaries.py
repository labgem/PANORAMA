#!/usr/bin/env python3
# coding:utf-8

"""
This module provides functions to write, update and erase a pangenome data from HDF5 files.
"""

# default libraries
import logging
from typing import Tuple
# installed libraries
import tables
from tqdm import tqdm
from ppanggolin.formats.writeBinaries import write_status as super_write_status
from ppanggolin.formats.writeBinaries import erase_pangenome as super_erase_pangenome
from ppanggolin.formats.writeBinaries import write_pangenome as super_write_pangenome

# local libraries
from panorama.pangenomes import Pangenome


def system_desc(max_id_len: int = 1, max_name_len: int = 1, max_gf_name_len: int = 1,
                max_metadata_source: int = 1, with_canonic: bool = False) -> dict:
    """

    Args:
        max_id_len: Maximum size of system name
        max_name_len: Maximum size of gene family name
        max_gf_name_len: Maximum size of gene family name
        max_metadata_source: Maximum size of annotation source name
        with_canonic: If true, include number of canonical

    Returns:
        Formated table
    """
    desc = {
        "ID": tables.StringCol(itemsize=max_id_len),
        "name": tables.StringCol(itemsize=max_name_len),
        "geneFam": tables.StringCol(itemsize=max_gf_name_len),
        "metadata_id": tables.Int64Col(),
        "metadata_source": tables.StringCol(itemsize=max_metadata_source),
    }
    if with_canonic:
        desc["canonical"] = tables.Int64Col()
    return desc


def get_system_len(pangenome: Pangenome, source: str
                   ) -> Tuple[Tuple[int, int, int, int, int], Tuple[int, int, int, int, int], int]:
    """
    Get maximum size of gene families information
    Args:
        pangenome: Pangenome filled with systems
        source: Name of the system source

    Returns:
        Maximum size of each element
    """

    def compare_len(sys, max_id_len, max_name_len, max_gf_name_len, max_annot_source_len) -> Tuple[int, int, int, int]:
        """
        Compare the length of elements to known maximum

        Args:
            sys: system to get info from
            max_id_len: Maximum size of the system identifier
            max_name_len: Maximum size of system name
            max_gf_name_len: Maximum size of gene family name
            max_annot_source_len: Maximum size of annotation source name

        Returns:
            Maximum length of each element
        """

        if len(sys.ID) > max_id_len:
            max_id_len = len(sys.name)
        if len(sys.name) > max_name_len:
            max_name_len = len(sys.name)
        for annot_source in sys.annotation_sources():
            if len(annot_source) > max_annot_source_len:
                max_annot_source_len = len(annot_source)
        for gf in sys.families:
            if len(gf.name) > max_gf_name_len:
                max_gf_name_len = len(gf.name)
        return max_id_len, max_name_len, max_gf_name_len, max_annot_source_len

    max_len_sys = (1, 1, 1, 1)
    max_len_canonical = (1, 1, 1, 1)
    expected_rows_sys, expected_rows_can, expected_rows_sys2can = (0, 0, 0)

    for system in pangenome.get_system_by_source(source):
        max_len_sys = compare_len(system, *max_len_sys)
        expected_rows_sys += len(system)
        for canonical in system.canonical:
            max_len_canonical = compare_len(canonical, *max_len_canonical)
            expected_rows_can += len(canonical)

    return (*max_len_sys, expected_rows_sys), (*max_len_canonical, expected_rows_can), expected_rows_sys2can


def write_systems(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing all systems detected in pangenome

   Args:
        pangenome: Pangenome with systems detected
        h5f: HDF5 file to write systems
        source: source of the systems
        disable_bar: Flag to disable progress bar
    """
    if '/systems' not in h5f:
        systems_group = h5f.create_group("/", "systems", "Detected systems")
    else:
        systems_group = h5f.root.systems
    system_len, canonical_len, expected_rows_sys2can = get_system_len(pangenome, source)
    source_group = h5f.create_group(systems_group, source, f"Detected systems from source: {source}")
    source_group._v_attrs.metadata_sources = pangenome.systems_sources_to_metadata_source()[source]
    system_table = h5f.create_table(source_group, "system", description=system_desc(*system_len[:-1], True),
                                    expectedrows=system_len[-1])
    canonical_table = h5f.create_table(source_group, "canonic", description=system_desc(*canonical_len[:-1]),
                                       expectedrows=canonical_len[-1])
    sys2canonical_table = h5f.create_table(source_group, "system_to_canonical",
                                           description={"system": tables.StringCol(itemsize=system_len[0]),
                                                        "canonic": tables.StringCol(itemsize=canonical_len[0])},
                                           expectedrows=expected_rows_sys2can)
    system_row = system_table.row
    canonical_row = canonical_table.row
    sys2canonical_row = sys2canonical_table.row
    canonic_seen = set()
    with tqdm(total=pangenome.number_of_systems(source=source, with_canonical=False), unit="system",
              disable=disable_bar) as progress:
        for system in pangenome.systems:
            for gf in system.families:
                system_row["ID"] = system.ID
                system_row["name"] = system.name
                system_row["geneFam"] = gf.name
                system_row["metadata_source"], system_row["metadata_id"] = system.get_metainfo(gf)
                system_row["canonical"] = len(system.canonical)
                system_row.append()
            progress.update()
            for canonical in system.canonical:
                if canonical.ID not in canonic_seen:
                    canonic_seen.add(canonical.ID)
                    for gf in canonical.families:
                        canonical_row["geneFam"] = gf.name
                        canonical_row["ID"] = canonical.ID
                        canonical_row["name"] = canonical.name
                        canonical_row["metadata_source"], system_row["metadata_id"] = canonical.get_metainfo(gf)
                        canonical_row.append()
                sys2canonical_row["system"] = system.ID
                sys2canonical_row["canonic"] = canonical.ID
                sys2canonical_row.append()
        system_table.flush()
        canonical_table.flush()
        sys2canonical_table.flush()
    logging.getLogger("PANORAMA").debug(f"Write {pangenome.number_of_systems(source=source, with_canonical=False)} "
                                        f"systems, with {len(canonic_seen)} canonical systems.")


def write_status(pangenome: Pangenome, h5f: tables.File):
    """
    Write pangenome status in HDF5 file

    Args:
        pangenome: Pangenome object
        h5f: Pangenome file
    """
    super_write_status(pangenome, h5f)  # call write_status from ppanggolin
    status_group = h5f.root.status

    if "systems" in pangenome.status:
        if not hasattr(status_group._v_attrs, 'systems_sources'):
            status_group._v_attrs.systems_sources = set()
        if pangenome.status["systems"] in ["Computed", "Loaded"]:
            status_group._v_attrs.systems = True
            status_group._v_attrs.systems_sources |= pangenome.systems_sources
    else:
        status_group._v_attrs.systems = False
        status_group._v_attrs.systems_sources = set()


def erase_pangenome(pangenome: Pangenome, graph: bool = False, gene_families: bool = False, partition: bool = False,
                    rgp: bool = False, spots: bool = False, modules: bool = False, metadata: bool = False,
                    systems: bool = False, source: str = None):
    """
    Erases tables from a pangenome .h5 file

    Args:
        pangenome: Pangenome to erase information
        graph: remove graph information
        gene_families: remove gene families information
        partition: remove partition information
        rgp: remove rgp information
        spots: remove spots information
        modules: remove modules information
        metadata: remove metadata
        systems: remove systems
        source: source of system or metadata
    """

    h5f = tables.open_file(pangenome.file, "a")
    status_group = h5f.root.status

    super_erase_pangenome(pangenome, graph, gene_families, partition, rgp, spots, modules,
                          metadata=metadata, metatype="families" if metadata else None, source=source)

    if '/systems' in h5f and systems:
        assert source is not None
        systems_group = h5f.root.systems
        if source in systems_group:
            logging.getLogger("PANORAMA").info(f"Erasing the formerly computed systems from source {source}")
            h5f.remove_node("/systems", source, recursive=True)
            status_group._v_attrs.systems_sources.remove(source)
            pangenome.status["systems_sources"].remove(f"{source}")
        if len(status_group._v_attrs.systems_sources) == 0:
            h5f.remove_node("/", "systems")
            status_group._v_attrs.systems = False
            pangenome.status["systems"] = "No"
    h5f.close()


def write_pangenome(pangenome: Pangenome, file_path: str, source: str = None,
                    force: bool = False, disable_bar: bool = False):
    """
    Writes or updates a pangenome file

    Args:
        pangenome: Pangenome object to write in file
        file_path: HDF5 file to save pangenome if not given the original file is used
        source: source of systems or metadata
        force: Flag to force overwrite
        disable_bar: Flag to disable progress bar
    """
    super_write_pangenome(pangenome, file_path, force, disable_bar)

    h5f = tables.open_file(file_path, "a")

    if "systems" in pangenome.status and pangenome.status["systems"] == "Computed":
        assert source is not None
        logging.getLogger("PANORAMA").info("Writing detected systems...")
        write_systems(pangenome=pangenome, h5f=h5f, source=source, disable_bar=disable_bar)
        pangenome.status["systems"] = "Loaded"

    write_status(pangenome, h5f)

    h5f.close()
