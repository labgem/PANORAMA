#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from tqdm import tqdm
from typing import List
from pathlib import Path
# installed libraries
import tables
from ppanggolin.formats.writeBinaries import write_status as super_write_status
from ppanggolin.formats.writeBinaries import erase_pangenome as super_erase_pangenome

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.system import System
from panorama.pangenomes import Pangenome


def annot_desc(max_name_len: int = 1, max_accession_len: int = 1,
               max_secondary_names_len: int = 1, max_description_len: int = 1,
               max_gf_name_len: int = 1) -> dict:
    """
    Create a formated table for gene families annotation description
    :param max_name_len:
    :param max_accession_len:
    :param max_secondary_names_len:
    :param max_description_len:
    :param max_gf_name_len:

    :return: Formated table
    """
    return {
        "accession": tables.StringCol(itemsize=max_accession_len),
        "name": tables.StringCol(itemsize=max_name_len),
        "secondary_names": tables.StringCol(itemsize=max_secondary_names_len),
        "description": tables.StringCol(itemsize=max_description_len),
        "score": tables.Float64Col(),
        "e_val": tables.Float64Col(),
        "bias": tables.Float64Col(),
        "geneFam": tables.StringCol(itemsize=max_gf_name_len)
    }


def get_annot_len(select_gf: List[GeneFamily], source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_gf: selected gene families from source
    :param source: Name of the annotation source
    :return: Maximum size of each element
    """
    max_name_len, max_accession_len, max_secondary_names_len, max_description_len = (1, 1, 1, 1)
    max_gf_name_len, expected_rows = (1, 0)

    for gf in select_gf:
        if len(gf.name) > max_gf_name_len:
            max_gf_name_len = len(gf.name)
        for annotation in gf.get_source(name=source):
            if len(annotation.name) > max_name_len:
                max_name_len = len(annotation.name)
            if len(annotation.accession) > max_accession_len:
                max_accession_len = len(annotation.accession)
            if len(annotation.secondary_names) > max_secondary_names_len:
                max_secondary_names_len = len(annotation.secondary_names)
            if len(annotation.description) > max_description_len:
                max_description_len = len(annotation.description)
            expected_rows += 1

    return max_name_len, max_accession_len, max_secondary_names_len, max_description_len, max_gf_name_len, expected_rows


def write_gene_fam_annot(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param source: name of the annotation source
    :param disable_bar: Disable progress bar
    """
    if '/geneFamiliesAnnot' not in h5f:
        annot_group = h5f.create_group("/", "geneFamiliesAnnot", "Gene families functional annotations")
    else:
        annot_group = h5f.root.geneFamiliesAnnot
    select_gf = list(pangenome.get_gf_by_sources(source=source))
    annot_len = get_annot_len(select_gf, source)
    source_table = h5f.create_table(annot_group, source, annot_desc(*annot_len[:-1]),
                                    expectedrows=annot_len[-1])
    annot_row = source_table.row
    for gf in tqdm(select_gf, unit='Gene family', desc=f'Source = {source}', disable=disable_bar):
        for annotation in gf.get_source(name=source):
            annot_row["name"] = annotation.name
            annot_row["accession"] = annotation.accession
            annot_row["secondary_names"] = annotation.secondary_names
            annot_row["description"] = annotation.description
            annot_row["score"] = annotation.score
            annot_row["e_val"] = annotation.e_val
            annot_row["bias"] = annotation.bias
            annot_row["geneFam"] = gf.name
            annot_row.append()
    source_table.flush()


def system_desc(max_id_len: int = 1, max_name_len: int = 1, max_canonical_len: int = 1,
                max_gf_name_len: int = 1) -> dict:
    """
    Create a formated table for detected systems in pangenome
    :param max_name_len:
    :param max_accession_len:
    :param max_secondary_names_len:
    :param max_description_len:
    :param max_gf_name_len:

    :return: Formated table
    """
    return {
        "ID": tables.StringCol(itemsize=max_id_len),
        "name": tables.StringCol(itemsize=max_name_len),
        "canonical": tables.StringCol(itemsize=max_canonical_len),
        "geneFam": tables.StringCol(itemsize=max_gf_name_len)
    }


def get_system_len(pangenome: Pangenome, source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param select_gf: selected gene families from source
    :param source: Name of the annotation source
    :return: Maximum size of each element
    """

    def compare_len(system: System, max_id_len, max_name_len, max_gf_name_len):
        if len(system.ID) > max_id_len:
            max_id_len = len(system.name)
        if len(system.name) > max_name_len:
            max_name_len = len(system.name)
        for gf in system.gene_families:
            if len(gf.name) > max_gf_name_len:
                max_gf_name_len = len(gf.name)
        return max_id_len, max_name_len, max_gf_name_len

    max_id_len, max_name_len, max_canonical_len, max_gf_name_len, expected_rows = (1, 1, 1, 1, 0)

    for system in pangenome.get_system_by_source(source):
        max_id_len, max_name_len, max_gf_name_len = compare_len(system, max_id_len, max_name_len, max_gf_name_len)
        canonical_name_len = 0
        for canonical in system.canonical:
            canonical_name_len += len(canonical.name)
            max_id_len, max_name_len, max_gf_name_len = compare_len(canonical, max_id_len,
                                                                    max_name_len, max_gf_name_len)
            expected_rows += len(canonical.gene_families)
        if canonical_name_len > max_canonical_len:
            max_canonical_len = canonical_name_len + len(system.canonical) - 1
        expected_rows += len(system.gene_families)

    return max_id_len, max_name_len, max_canonical_len, max_gf_name_len, expected_rows


def write_systems(pangenome: Pangenome, h5f: tables.File, source: str, disable_bar: bool = False):
    """
    Writing a table containing all systems detected in pangenome

    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param disable_bar: Disable progress bar
    """
    if '/systems' not in h5f:
        systems_group = h5f.create_group("/", "systems", "Detected systems")
    else:
        systems_group = h5f.root.systems
    system_len = get_system_len(pangenome, source)
    source_table = h5f.create_table(systems_group, source, description=system_desc(*system_len[:-1]),
                                    expectedrows=system_len[-1])
    source_row = source_table.row
    for system in tqdm(pangenome.systems, total=len(list(pangenome.get_system_by_source(source))),
                       unit="system", disable=disable_bar):
        for gf in system.gene_families:
            source_row["ID"] = system.ID
            source_row["name"] = system.name
            source_row["geneFam"] = gf.name
            source_row["canonical"] = ",".join([canonical.name for canonical in system.canonical])
            source_row.append()
        for canonical in system.canonical:
            for gf in canonical.gene_families:
                source_row["geneFam"] = gf.name
                source_row["ID"] = canonical.ID
                source_row["name"] = canonical.name
                source_row.append()
    source_table.flush()


def write_status(pangenome: Pangenome, h5f: tables.File):
    """
    Write pangenome status in HDF5 file
    :param pangenome: Pangenome object
    :param h5f: Pangenome file
    """
    super_write_status(pangenome, h5f)  # call write_status from ppanggolin
    status_group = h5f.root.status

    status_group._v_attrs.annotations = True if pangenome.status["annotations"] in ["Computed", "Loaded",
                                                                                    "inFile"] else False
    if hasattr(status_group._v_attrs, "annotations_sources") and pangenome.status["annotations"] in ["Computed",
                                                                                                     "Loaded"]:
        if len(status_group._v_attrs.annotations_sources) > 0:
            status_group._v_attrs.annotations_sources.add(*pangenome.annotations_sources)
        else:
            status_group._v_attrs.annotations_sources = pangenome.systems_sources

    status_group._v_attrs.systems = True if pangenome.status["systems"] in ["Computed", "Loaded",
                                                                            "inFile"] else False
    if hasattr(status_group._v_attrs, "systems_sources") and pangenome.status["systems"] in ["Computed", "Loaded"]:
        if len(status_group._v_attrs.systems_sources) > 0:
            status_group._v_attrs.systems_sources.add(*pangenome.systems_sources)
        else:
            status_group._v_attrs.systems_sources = pangenome.systems_sources


def erase_pangenome(pangenome: Pangenome, graph: bool = False, gene_families: bool = False, partition: bool = False,
                    rgp: bool = False, spots: bool = False, modules: bool = False, annotations: bool = False,
                    systems: bool = False, source: str = None):
    """
        Erases tables from a pangenome .h5 file
        :param pangenome: Pangenome
        :param graph: remove graph information
        :param gene_families: remove gene families information
        :param partition: remove partition information
        :param rgp: remove rgp information
        :param spots: remove spots information
        :param modules: remove modules information
        :param annotations: remove annotation
        :param source: annotation source
        """

    h5f = tables.open_file(pangenome.file, "a")
    status_group = h5f.root.status

    super_erase_pangenome(pangenome, graph, gene_families, partition, rgp, spots, modules)

    if '/geneFamiliesAnnot' in h5f and annotations:
        assert source is not None
        annotations_group = h5f.root.geneFamiliesAnnot
        if source in annotations_group:
            logging.getLogger().info(f"Erasing the formerly computed annotations from source {source}")
            h5f.remove_node("/geneFamiliesAnnot", source)
            status_group._v_attrs.annotations_sources.remove(source)
            pangenome.status["annotations_sources"].remove(f"{source}")
        if len(status_group._v_attrs.annotations_sources) == 0:
            h5f.remove_node("/", "geneFamiliesAnnot")
            status_group._v_attrs.annotations = False
            pangenome.status["annotations"] = "No"
    if '/systems' in h5f and systems:
        assert source is not None
        systems_group = h5f.root.systems
        if source in systems_group:
            logging.getLogger().info(f"Erasing the formerly computed systems from source {source}")
            h5f.remove_node("/systems", source)
            status_group._v_attrs.systems_sources.remove(source)
            pangenome.status["systems_sources"].remove(f"{source}")
        if len(status_group._v_attrs.systems_sources) == 0:
            h5f.remove_node("/", "systems")
            status_group._v_attrs.systems = False
            pangenome.status["systems"] = "No"
    h5f.close()


def write_pangenome(pangenome: Pangenome, file_path: str, disable_bar: bool = False, source: str = None):
    """
    Writes or updates a pangenome file

    :param pangenome: pangenome object
    :param file_path: HDF5 file to save pangenome if not given the original file is used
    :param disable_bar: Allow to disable progress bar
    :param source: annotation source
    """

    h5f = tables.open_file(file_path, "a")

    if pangenome.status["annotations"] == "Computed":
        assert source is not None
        logging.getLogger().info("Writing gene families annotations...")
        write_gene_fam_annot(pangenome=pangenome, h5f=h5f, source=source, disable_bar=disable_bar)
        pangenome.status["annotations"] = "Loaded"

    if pangenome.status["systems"] == "Computed":
        assert source is not None
        logging.getLogger().info("Writing detected systems...")
        write_systems(pangenome=pangenome, h5f=h5f, source=source, disable_bar=disable_bar)
        pangenome.status["systems"] = "Loaded"

    write_status(pangenome, h5f)
    h5f.close()
