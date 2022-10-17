#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging

from panorama.geneFamily import GeneFamily
from tqdm import tqdm

# installed libraries
import tables


# local libraries
from panorama.pangenomes import Pangenome


def gene_annot_desc(max_annotation_len: int = 1, max_fam_len: int = 1) -> dict:
    """
    Create a formated table for gene families description
    :param max_name_len: Maximum size of gene family name
    :param max_sequence_length: Maximum size of gene family representing gene sequences
    :param max_part_len: Maximum size of gene family partition
    :return: Formated table
    """
    return {
        "annotation": tables.StringCol(itemsize=max_annotation_len),
        "geneFam": tables.StringCol(itemsize=max_fam_len)
    }


def get_annot_len(pangenome: Pangenome, source: str) -> (int, int):
    """
    Get maximum size of gene families information
    :param pangenome: Pangenome with gene families computed
    :param source: Name of the annotation source
    :return: Maximum size of each element
    """
    max_annotation_len = 1
    max_fam_len = 1
    expected_rows = 0

    for genefam in pangenome.get_gf_by_annnotation(source=source):
        if len(genefam.name) > max_fam_len:
            max_fam_len = len(genefam.name)
        for annot in genefam.annotation[source]:
            if len(annot) > max_annotation_len:
                max_annotation_len = len(annot)
            expected_rows += 1

    return max_annotation_len, max_fam_len, expected_rows


def write_gene_fam_annot(pangenome: Pangenome, h5f: tables.File, force: bool = False, disable_bar: bool = False):
    """
    Writing a table containing the protein sequences of each family
    :param pangenome: Pangenome with gene families computed
    :param h5f: HDF5 file to write gene families
    :param force: force to write information if precedent information exist
    :param disable_bar: Disable progress bar
    """
    if '/geneFamiliesAnnot' in h5f and force is True:
        logging.getLogger().info("Erasing the formerly computed gene family annotations...")
        h5f.remove_node('/', 'geneFamiliesAnnot', recursive=True)  # erasing the table, and rewriting a new one.
    annot_group = h5f.create_group("/", "geneFamiliesAnnot", "Gene families functional annotation")
    for source in pangenome.annotation_source:
        max_annotation_len, max_fam_len, expected_rows = get_annot_len(pangenome, source)
        source_table = h5f.create_table(annot_group, source, gene_annot_desc(max_annotation_len, max_fam_len),
                                        expectedrows=expected_rows)
        annot_row = source_table.row
        for genefam in pangenome.get_gf_by_annnotation(source=source):
            for annot in genefam.annotation[source]:
                annot_row["annotation"] = annot
                annot_row["geneFam"] = genefam.name
                annot_row.append()
        source_table.flush()
