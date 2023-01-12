#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from copy import copy
import itertools
import logging
from pathlib import Path
from tqdm import tqdm
from typing import List

# installed libraries
import numpy as np
import pandas as pd

# local libraries
from annotate.hmm_search import profile_gfs
from format.read_binaries import check_pangenome_info
from utils import check_tsv_sanity, mkdir
from system import System
from geneFamily import GeneFamily
from pangenomes import Pangenome


def write_annotation_to_families(pangenome: Pangenome, output: Path, disable_bar: bool = False):
    """ Write a tsv file with all annotation and sources present in pangenomes

    :param pangenome: Pangenome with all annotations
    :param output: Path to output directory
    :param disable_bar: allow to disable the progress bar
    """
    nb_source = len(pangenome.annotation_source)
    source_list = list(pangenome.annotation_source)
    column_name = np.array(
        f"Pangenome,Gene Family,{','.join([f'Annotation {source},Accession {source},Secondary names {source}' for source in source_list])}".split(
            ','))
    array_list = []
    for gf in tqdm(pangenome.gene_families, unit='gene families', disable=disable_bar):
        annot_array = np.empty((gf.max_annotation_by_source()[1], 2 + nb_source * 3), dtype=object)
        if annot_array.shape[0] > 0:
            annot_array[:, 0] = pangenome.name
            annot_array[:, 1] = gf.name
            index_source = 2
            for source in source_list:
                index_annot = 0
                if source in gf.sources:
                    for annotation in gf.get_source(source):
                        annot_array[index_annot, index_source] = annotation.name
                        annot_array[index_annot, index_source + 1] = annotation.accession
                        annot_array[
                            index_annot, index_source + 2] = annotation.secondary_names if annotation.secondary_names is not pd.NA else '-'
                        index_annot += 1
                index_source += 3
            array_list.append(annot_array)
    out_df = pd.DataFrame(np.concatenate(array_list), columns=column_name)
    out_df = out_df.sort_values(by=['Pangenome', 'Gene Family'] + list(column_name[range(3, len(column_name), 2)]))
    out_df.to_csv(f"{output}/families_annotations.tsv", sep="\t", header=True)


def write_hmm(gf: GeneFamily, output: Path):
    """ Write an HMM profile for a gene family

    :param gf: Gene family with HMM
    :param output: Path to output directory
    """
    with open(f"{output.absolute().as_posix()}/{gf.name}.hmm", 'wb') as hmm_file:
        gf.HMM.write(hmm_file)


def write_hmm_profile(pangenome: Pangenome, output: Path, threads: int = 1,
                      disable_bar: bool = False):
    """Write an HMM profile for all gene families in pangenome

    :param pangenome: Pangenome with gene families
    :param output: Path to output directory
    :param threads: number of threads available
    :param disable_bar: allow to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=pangenome.number_of_gene_families(), unit='file',
                  desc='write gene families hmm/profile', disable=disable_bar) as progress:
            futures = []
            for gf in pangenome.gene_families:
                if gf.HMM is not None:
                    future = executor.submit(write_hmm, gf, output)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

                for future in futures:
                    future.result()


def write_sys_to_organism(system, organism, projection, mandatory_family, accessory_family, fu) -> bool:
    """ Project system to organism

    :param system: Detected system
    :param organism: organism in pangenome
    :param projection: projection list output
    :param mandatory_family: list of mandatory families in system
    :param accessory_family: list of accessory family in system
    :param fu: functional unit

    :return: True if a projection occurs in organism
    """
    def search_diagonal(matrix: np.ndarray, min_necessary: int = None) -> List[List[int]]:
        """Search any possible diagonals in matrix to confirm system projection

        :param matrix: matrix with present family in system and pangenome
        :param min_necessary: minimum length of diagonal to project

        :return: list of diagonals found
        """
        r_diags = []
        min_necessary = matrix.shape[0] if min_necessary is None else min_necessary
        nonull_matrix = matrix[:, ~np.all(matrix == 0, axis=0)]
        diags = [nonull_matrix.diagonal(i).tolist() for i in range(-min_necessary, min_necessary, 1)
                 if len(nonull_matrix.diagonal(i).tolist()) >= min_necessary]
        for comb_cols in itertools.permutations(range(0, nonull_matrix.shape[1]), 2):
            copy_matrix = copy(nonull_matrix)
            copy_matrix[:, [comb_cols[1], comb_cols[0]]] = nonull_matrix[:, comb_cols]
            diags.extend(copy_matrix.diagonal(i).tolist() for i in range(-min_necessary, min_necessary, 1)
                         if len(copy_matrix.diagonal(i).tolist()) >= min_necessary)
        for diag in diags:
            for comb in [diag[i: i + min_necessary] for i in range(0, len(diag) - min_necessary + 1)]:
                if all(x > 0 for x in comb):
                    r_diags.append(comb)
        r_diags.sort()
        return list(r_diag for r_diag, _ in itertools.groupby(r_diags))
    if organism.name == "GCF_000773865.1_ASM77386v1_genomic" and system.name == "DMS_other":
        print("carapuce")
    org_projection = [system.ID, system.name, organism.name]
    select_families = system.gene_families.intersection(organism.families)
    mandatory_index = {mandatory: index for index, mandatory in enumerate(mandatory_family, start=1)}
    accessory_index = {accessory: index for index, accessory in enumerate(accessory_family,
                                                                          start=1 + len(mandatory_family))}
    nb_annotation = len(mandatory_family) + len(accessory_family)
    bool_projection = False
    for comb_gf in itertools.permutations(select_families, fu.parameters['min_total']):
        mandatory_matrix = np.zeros((len(comb_gf), len(mandatory_family)))
        all_matrix = np.concatenate((mandatory_matrix, np.zeros((len(comb_gf), len(accessory_family)))), axis=1)
        coordinate_dict = {}
        for index_gf, gf in enumerate(comb_gf):
            for index_annot, annot in enumerate(mandatory_family + accessory_family, start=1):
                coordinate_dict[index_gf * nb_annotation + index_annot] = (gf, annot)
            if gf.get_source(system.source) is not None:
                for annotation in [annot.name for annot in gf.get_source(system.source)]:
                    if annotation in mandatory_family:
                        mandatory_matrix[index_gf][mandatory_index[annotation] - 1] = index_gf * nb_annotation + mandatory_index[annotation]
                        all_matrix[index_gf][mandatory_index[annotation] - 1] = index_gf * nb_annotation + mandatory_index[annotation]
                    elif annotation in accessory_family:
                        all_matrix[index_gf][accessory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                accessory_index[annotation]
        if len(search_diagonal(mandatory_matrix, fu.parameters['min_mandatory'])) > 0:
            search_diag = search_diagonal(all_matrix)
            if len(search_diag) > 0:
                for diag in search_diag:
                    for coord in diag:
                        gf, annotation = coordinate_dict[coord]
                        genes = gf.get_genes_per_org(organism)
                        for gene in genes:
                            projection.append(org_projection + [gf.name, gf.named_partition, annotation,
                                                                gene.name, gene.start, gene.stop, gene.strand])
                bool_projection = True
    return bool_projection


def write_system_projection(system: System) -> pd.DataFrame:
    """Project system on all pangenome organisms

    :param system: system to project

    :return: Dataframe result of projection
    """
    projection = []
    mandatory_family = [fam.name for fam in system.model.families if fam.type == 'mandatory']
    accessory_family = [fam.name for fam in system.model.families if fam.type == 'accessory']
    system_orgs = set.union(*list([gf.organisms for gf in system.gene_families]))
    fu = list(system.model.func_units)[0]
    for organism in system_orgs:
        find_projection = write_sys_to_organism(system, organism, projection, mandatory_family, accessory_family, fu)
        if not find_projection:
            for canonical_system in system.canonical:
                canonical_mandatory_family = [fam.name for fam in canonical_system.model.families
                                              if fam.type == 'mandatory']
                canonical_accessory_family = [fam.name for fam in canonical_system.model.families
                                              if fam.type == 'accessory']
                canonical_fu = list(canonical_system.model.func_units)[0]
                write_sys_to_organism(canonical_system, organism, projection, canonical_mandatory_family,
                                      canonical_accessory_family, canonical_fu)
    projection_df = pd.DataFrame(projection).drop_duplicates()
    return projection_df


def write_systems_projection(pangenome: Pangenome, output: Path, threads: int = 1,
                             force: bool = False, disable_bar: bool = False):
    """ Write all systems in pangenomes and project on organisms

    :param pangenome: Pangenome to project
    :param output: Path to output directory
    :param threads:  Number of threads available
    :param force: Force to write into the output directory
    :param disable_bar: Allow to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.getLogger().info('Write system projection')
        with tqdm(total=pangenome.number_of_systems(), unit='system', disable=disable_bar) as progress:
            futures = []
            for system in pangenome.systems:
                future = executor.submit(write_system_projection, system)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

                results = []
                for future in futures:
                    result = future.result()
                    results.append(result)
    systems_projection = pd.DataFrame()
    for result in results:
        systems_projection = pd.concat([systems_projection, result], ignore_index=True)
    systems_projection.columns = ["system number", "system name", "organism", "gene family", "partition",
                                  "annotation name", "gene", "start", "stop", "strand"]
    systems_projection.sort_values(by=["system number", "system name", "organism", "start", "stop"],
                                   ascending=[True, True, True, True, True], inplace=True)
    mkdir(f"{output}/projection_systems1", force=force)
    for organism_name in systems_projection["organism"].unique():
        org_df = systems_projection.loc[systems_projection["organism"] == organism_name]
        org_df = org_df.iloc[:, [0, 1, 3, 4, 5, 6, 7, 8]]
        org_df.sort_values(by=["system number", "system name", "start", "stop"],
                           ascending=[True, True, True, True], inplace=True)
        org_df.to_csv(f"{output}/projection_systems1/{organism_name}.tsv", sep="\t", index=False)
    systems_projection.to_csv(f"{output}/systems1.tsv", sep="\t", index=False)


def write_flat_files(pangenome, output: Path, annotation: bool = False, hmm: bool = False,
                     disable_bar: bool = False, **kwargs):
    """Launcher to write flat file from pangenomes

    :param pangenome: Pangenome with information to write
    :param output: Path to output directory
    :param annotation: Launch annotation write function
    :param hmm: Launch hmm write function
    :param disable_bar: disable progress bar
    """
    need_annotations = False
    need_families = False
    need_graph = False
    need_partitions = False
    need_spots = False
    need_regions = False
    need_modules = False
    need_gene_sequences = False
    need_annotations_fam = False

    if annotation:
        need_families = True
        need_annotations_fam = True

    if hmm:
        need_families = True
        need_annotations = True

    check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                         need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_regions,
                         need_spots=need_spots, need_gene_sequences=need_gene_sequences, need_modules=need_modules,
                         need_annotation_fam=need_annotations_fam, disable_bar=disable_bar)

    if annotation:
        write_annotation_to_families(pangenome, output)

    if hmm:
        msa_path = kwargs.get('msa_path', None)
        msa_format = kwargs.get('msa_format', 'afa')
        threads = kwargs.get('threads', 1)
        profile_gfs(pangenome, msa_path, msa_format, threads, disable_bar)
        write_hmm_profile(pangenome, output, threads, disable_bar)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for pangenome_name, pangenome_info in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
        pangenome.add_file(pangenome_info["path"])
        write_flat_files(pangenome, output=args.output, annotation=args.annotation, hmm=args.hmm,
                         msa_path=args.msa, threads=args.threads, disable_bar=args.disable_prog_bar)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("write", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_write(parser)
    return parser


def parser_write(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--annotation", required=False, action="store_true",
                          help="Write all the annotations from families")
    optional.add_argument("--hmm", required=False, action="store_true",
                          help="Write an hmm for each gene families in pangenomes")
    optional.add_argument("--msa", required=False, type=Path, default=None,
                          help="To create a HMM profile for families, you can give a msa of each gene in families."
                               "This msa could be get from ppanggolin (See ppanggolin msa). "
                               "If no msa provide Panorama will launch one.")
    optional.add_argument("--msa-format", required=False, type=str, default="afa",
                          choices=["stockholm", "pfam", "a2m", "psiblast", "selex", "afa",
                                   "clustal", "clustallike", "phylip", "phylips"],
                          help="Format of the input MSA.")
    optional.add_argument("--threads", required=False, type=int, default=1)


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_write(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
