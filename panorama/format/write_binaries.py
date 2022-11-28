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
    :param disable_bar: Disable progress bar
    """
    if '/geneFamiliesAnnot' not in h5f:
        annot_group = h5f.create_group("/", "geneFamiliesAnnot", "Gene families functional annotation")
    else:
        annot_group = h5f.root.geneFamiliesAnnot
    if f'/geneFamiliesAnnot/{source}' in h5f:
        logging.getLogger().info(f"Erasing gene family annotations from {source}...")
        h5f.remove_node('/geneFamiliesAnnot', f'{source}')  # erasing the table, and rewriting a new one.
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

def write_status(pangenome: Pangenome, h5f: tables.File):
    """
    Write pangenome status in HDF5 file
    :param pangenome: Pangenome object
    :param h5f: Pangenome file
    """
    super_write_status(pangenome, h5f)  # call write_status from ppanggolin
    status_group = h5f.root.status

    status_group._v_attrs.annotation = True if pangenome.status["annotation"] in ["Computed", "Loaded",
                                                                                  "inFile"] else False
    status_group._v_attrs.annotation_source = pangenome.annotation_source

def erase_pangenome(pangenome: Pangenome, graph: bool = False, gene_families: bool = False, partition: bool = False,
                    rgp: bool = False, spots: bool = False, modules: bool = False, annotation: bool = False, **kargs):
    """
        Erases tables from a pangenome .h5 file
        :param pangenome: Pangenome
        :param graph: remove graph information
        :param gene_families: remove gene families information
        :param partition: remove partition information
        :param rgp: remove rgp information
        :param spots: remove spots information
        :param modules: remove modules information
        """

    h5f = tables.open_file(pangenome.file, "a")
    status_group = h5f.root.status
    info_group = h5f.root.info

    super_erase_pangenome(pangenome, graph, gene_families, partition, rgp, spots, modules)

    if '/geneFamiliesAnnot' in h5f and annotation:
        assert 'source' in kargs.keys()
        annotation_table = h5f.root.geneFamiliesAnnot
        if kargs['source'] in annotation_table:
            logging.getLogger().info(f"Erasing the formerly computed annotation from source {kargs['source']}")
            h5f.remove_node("/geneFamiliesAnnot", f"{kargs['source']}")
            status_group._v_attrs.annotation_source.remove(kargs['source'])
            pangenome.status["annotation_source"].remove(f"{kargs['source']}")
        if len(status_group._v_attrs.annotation_source) == 0:
            h5f.remove_node("/", "geneFamiliesAnnot")
            status_group._v_attrs.annotation = False
            pangenome.status["annotation"] = "No"
    h5f.close()


def write_pangenome(pangenome: Pangenome, file_path: Path, disable_bar: bool = False, **kargs):
    """
    Writes or updates a pangenome file

    :param pangenome: pangenome object
    :param file_path: HDF5 file to save pangenome if not given the original file is used
    :param force: force to write on pangenome if information already exist
    :param disable_bar: Allow to disable progress bar
    """

    h5f = tables.open_file(file_path, "a")

    if pangenome.status["annotation"] == "Computed":
        assert 'source' in kargs.keys()
        logging.getLogger().info("Writing gene families annotations...")
        write_gene_fam_annot(pangenome=pangenome, h5f=h5f, source=kargs['source'], disable_bar=disable_bar)
        pangenome.status["annotation"] = "Loaded"

    write_status(pangenome, h5f)
    h5f.close()
