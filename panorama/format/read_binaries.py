#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from tqdm import tqdm

# installed libraries
import tables
from ppanggolin.formats import read_chunks
from ppanggolin.formats import check_pangenome_info as check_pp

# local libraries
from panorama.pangenomes import Pangenome


def read_gene_families_info(pangenome: Pangenome, disable_bar: bool = False):
    """
    Read information about gene families in pangenome hdf5 file to add in pangenome object
    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")

    annotation_group = h5f.root.geneFamiliesAnnot

    if pangenome.status["genesClustered"] != "Loaded":
        raise Exception("Gene families aren't loaded in pangenome.")

    for source_table in annotation_group:
        logging.getLogger().info(f"Reading annotation {source_table.name}")
        for row in tqdm(read_chunks(source_table, chunk=20000), total=source_table.nrows,
                        unit="gene family", disable=disable_bar):
            fam = pangenome.get_gene_family(row["geneFam"].decode())
            fam.add_annotation(source=source_table.name, annotation=row["annotation"].decode())


def check_pangenome_info(pangenome, need_annotations: bool = False, need_families: bool = False,
                         need_graph: bool = False, need_partitions: bool = False, need_rgp: bool = False,
                         need_spots: bool = False, need_gene_sequences: bool = False, need_modules: bool = False,
                         need_anntation_fam: bool = False, disable_bar: bool = False):
    """
    Defines what needs to be read depending on what is needed, and automatically checks if the required elements
    have been computed with regard to the `pangenome.status`
    :param pangenome: Pangenome object without some information
    :param need_annotations: get annotation
    :param need_families: get gene families
    :param need_graph: get graph
    :param need_partitions: get partition
    :param need_rgp: get RGP
    :param need_spots: get hotspot
    :param need_gene_sequences: get gene sequences
    :param need_modules: get modules
    :param disable_bar: Allow to disable the progress bar
    """
    if need_anntation_fam:
        if pangenome.status["genesClustered"] == "inFile":
            need_families = True
        elif pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no gene families. See the 'cluster' subcommand.")

    check_pp(pangenome, need_annotations, need_families, need_graph, need_partitions, need_rgp,
             need_spots, need_gene_sequences, need_modules, disable_bar=disable_bar)

    if need_anntation_fam:
        read_gene_families_info(pangenome)
