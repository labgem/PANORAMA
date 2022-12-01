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
from panorama.annotation import Annotation
from panorama.pangenomes import Pangenome


def get_status(pangenome, pangenome_file):
    """
        Checks which elements are already present in the file.
    """
    h5f = tables.open_file(pangenome_file, "r")
    logging.getLogger().info("Getting the current pangenome status")
    status_group = h5f.root.status
    if status_group._v_attrs.genomesAnnotated:
        pangenome.status["genomesAnnotated"] = "inFile"
    if status_group._v_attrs.genesClustered:
        pangenome.status["genesClustered"] = "inFile"
    if status_group._v_attrs.geneSequences:
        pangenome.status["geneSequences"] = "inFile"
    if status_group._v_attrs.geneFamilySequences:
        pangenome.status["geneFamilySequences"] = "inFile"
    if status_group._v_attrs.NeighborsGraph:
        pangenome.status["neighborsGraph"] = "inFile"

    if 'Partitionned' in status_group._v_attrs._f_list():
        # Partitionned keep working with older version
        if status_group._v_attrs.Partitionned:
            status_group._v_attrs.Partitioned = True
        del status_group._v_attrs.Partitionned

    if status_group._v_attrs.Partitioned:
        pangenome.status["partitioned"] = "inFile"

    if hasattr(status_group._v_attrs, "predictedRGP") and status_group._v_attrs.predictedRGP:
        pangenome.status["predictedRGP"] = "inFile"

    if hasattr(status_group._v_attrs, "spots") and status_group._v_attrs.spots:
        pangenome.status["spots"] = "inFile"

    if hasattr(status_group._v_attrs, "modules") and status_group._v_attrs.modules:
        pangenome.status["modules"] = "inFile"

    if "/info" in h5f:
        info_group = h5f.root.info
        pangenome.parameters = info_group._v_attrs.parameters

    if hasattr(status_group._v_attrs, "annotation") and status_group._v_attrs.annotation:
        pangenome.status["annotation"] = "inFile"
        pangenome.status["annotation_source"] = status_group._v_attrs.annotation_source
    h5f.close()


def read_gene_families_annotations(pangenome: Pangenome, h5f: tables.File, disable_bar: bool = False):
    """
    Read information about gene families in pangenome hdf5 file to add in pangenome object
    :param pangenome: Pangenome object without gene families information
    :param h5f: Pangenome HDF5 file with gene families information
    :param disable_bar: Disable the progress bar
    """
    annotation_table = h5f.root.geneFamiliesAnnot
    for source in annotation_table:
        for row in tqdm(read_chunks(source), total=source.nrows, unit="annotation",
                        desc=f"from source {source.name}", disable=disable_bar):
            gf = pangenome.get_gene_family(row["geneFam"].decode())
            annotation = Annotation(source=source.name, name=row['name'].decode(), accession=row["accession"].decode(),
                                    secondary_names=row["secondary_names"].decode(), score=row["score"],
                                    description=row["description"].decode(), e_val=row["e_val"], bias=row["bias"])
            gf.add_annotation(source=source.name, annotation=annotation)
    pangenome.status["annotation"] = "Loaded"


def check_pangenome_info(pangenome: Pangenome, need_annotations: bool = False, need_families: bool = False,
                         need_graph: bool = False, need_partitions: bool = False, need_rgp: bool = False,
                         need_spots: bool = False, need_gene_sequences: bool = False, need_modules: bool = False,
                         need_annotation_fam: bool = False, disable_bar: bool = False):
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
    if need_annotation_fam:
        if pangenome.status["genesClustered"] == "inFile":
            need_families = True
        elif pangenome.status["genesClustered"] not in ["Computed", "Loaded"]:
            raise Exception("Your pangenome has no gene families. See the 'cluster' subcommand of ppanggolin.")

    check_pp(pangenome, need_annotations, need_families, need_graph, need_partitions, need_rgp,
             need_spots, need_gene_sequences, need_modules, disable_bar=disable_bar)

    if hasattr(pangenome, "file"):
        filename = pangenome.file
    else:
        raise FileNotFoundError("The provided pangenome does not have an associated .h5 file")
    h5f = tables.open_file(filename, "r")
    if need_annotation_fam:
        read_gene_families_annotations(pangenome, h5f, disable_bar=disable_bar)
    h5f.close()
