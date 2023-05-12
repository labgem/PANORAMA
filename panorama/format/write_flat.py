#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import pathlib
from concurrent.futures import ThreadPoolExecutor
from copy import copy
import collections
import itertools
import logging
from pathlib import Path
from tqdm import tqdm
from typing import Collection, Dict, List, Set

# installed libraries
import numpy as np
import pandas as pd
import ppanggolin.metadata
from tqdm import tqdm
from ppanggolin.region import Region, Spot

# local libraries
from panorama.annotate.hmm_search import profile_gfs
from panorama.format.read_binaries import check_pangenome_info
from panorama.format.write_proksee import write_proksee
from panorama.utils import check_tsv_sanity, mkdir
from panorama.system import System
from panorama.region import Module
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome

# For Quentin
from collections import OrderedDict
from panorama.format.figure import per_pan_heatmap, heatmap, upsetplot, hbar_ID_type


def check_flat_parameters(args):
    if args.hmm and (args.msa is None or args.msa_format is None):
        raise Exception("To write HMM you need to give msa files and format")
    if args.systems or args.systems_asso is not None or "systems" in args.proksee:
        if args.models is None or args.sources is None:
            raise Exception("To read system and write flat related, "
                            "it's necessary to give models and source.")
        if len(args.models) != len(args.sources):
            raise Exception("To read systems from different sources you need to give "
                            "one models directory by corresponding source.")


def write_annotations_to_families(pangenome: Pangenome, output: Path, sources: Set[str], disable_bar: bool = False):
    """ Write a tsv file with all annotations and sources present in pangenomes

    :param pangenome: Pangenome with all annotations
    :param output: Path to output directory
    :param disable_bar: allow to disable the progress bar
    """
    nb_source = len(sources)
    source_list = list(sources)
    column_name = np.array(f"Pangenome,Gene Family,{','.join([f'Annotation {source},Accession {source},Secondary names {source}' for source in source_list])}".split(','))
    array_list = []
    for gf in tqdm(pangenome.gene_families, unit='gene families', disable=disable_bar):
        annot_array = np.empty((gf.max_metadata_by_source()[1], 2 + nb_source * 3), dtype=object)
        if annot_array.shape[0] > 0:
            annot_array[:, 0] = pangenome.name
            annot_array[:, 1] = gf.name
            index_source = 2
            for source in source_list:
                index_annot = 0
                if source in gf.sources:
                    for annotation in gf.get_source(source):
                        annotation: ppanggolin.metadata.Metadata
                        annot_array[index_annot, index_source] = annotation.get("protein_name")
                        annot_array[index_annot, index_source + 1] = annotation.get("Accession")
                        if annotation.get("secondary_names", skip_error=True) not in [pd.NA, None]:
                            annot_array[index_annot, index_source + 2] = annotation.get("secondary_names")
                        else:
                            annot_array[index_annot, index_source + 2] = '-'
                        index_annot += 1
                index_source += 3
            array_list.append(annot_array)
    out_df = pd.DataFrame(np.concatenate(array_list), columns=column_name)
    out_df = out_df.sort_values(by=['Pangenome', 'Gene Family'] + list(column_name[range(3, len(column_name), 2)]))
    out_df.to_csv(f"{output}/{pangenome.name}/families_annotations.tsv", sep="\t", header=True)


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

    org_projection = [system.ID, system.name, organism.name]
    select_families = system.gene_families.intersection(organism.families)
    mandatory_index = {mandatory: index for index, mandatory in enumerate(mandatory_family, start=1)}
    accessory_index = {accessory: index for index, accessory in enumerate(accessory_family,
                                                                          start=1 + len(mandatory_family))}
    nb_annotation = len(mandatory_family) + len(accessory_family)
    bool_projection = False
    for comb_gf in itertools.permutations(select_families, fu.min_total):
        mandatory_matrix = np.zeros((len(comb_gf), len(mandatory_family)))
        all_matrix = np.concatenate((mandatory_matrix, np.zeros((len(comb_gf), len(accessory_family)))), axis=1)
        coordinate_dict = {}
        for index_gf, gf in enumerate(comb_gf):
            for index_annot, annot in enumerate(mandatory_family + accessory_family, start=1):
                coordinate_dict[index_gf * nb_annotation + index_annot] = (gf, annot)
            if gf.get_source(system.source) is not None:
                for annotation in [annot.get("protein_name") for annot in gf.get_source(system.source)]:
                    if annotation in mandatory_family:
                        mandatory_matrix[index_gf][mandatory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                      mandatory_index[annotation]
                        all_matrix[index_gf][mandatory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                mandatory_index[annotation]
                    elif annotation in accessory_family:
                        all_matrix[index_gf][accessory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                accessory_index[annotation]
        if len(search_diagonal(mandatory_matrix, fu.min_mandatory)) > 0:
            search_diag = search_diagonal(all_matrix)
            if len(search_diag) > 0:
                for diag in search_diag:
                    for coord in diag:
                        gf, annotation = coordinate_dict[coord]
                        genes = gf.get_genes_per_org(organism)
                        for gene in genes:
                            projection.append(org_projection + [gf.name, gf.named_partition, annotation,
                                                                gene.name, gene.start, gene.stop, gene.strand, gene.is_fragment])
                bool_projection = True
    return bool_projection


def write_system_projection(system: System) -> pd.DataFrame:
    """Project system on all pangenome organisms

    :param system: system to project

    :return: Dataframe result of projection
    """

    projection, mandatory_family, accessory_family = ([], [], [])
    for fam in system.model.families:
        if fam.presence == 'mandatory':
            mandatory_family.append(fam.name)
            mandatory_family += fam.exchangeables
        if fam.presence == 'accessory':
            accessory_family.append(fam.name)
            accessory_family += fam.exchangeables

    system_orgs = set.union(*list([gf.organisms for gf in system.gene_families]))
    fu = list(system.model.func_units)[0]
    for organism in system_orgs:
        find_projection = write_sys_to_organism(system, organism, projection, mandatory_family, accessory_family, fu)
        if not find_projection:
            for canonical_system in system.canonical:

                canonical_mandatory_family, canonical_accessory_family = ([], [])
                for fam in canonical_system.model.families:
                    if fam.presence == 'mandatory':
                        canonical_mandatory_family.append(fam.name)
                        canonical_mandatory_family += fam.exchangeables
                    if fam.presence == 'accessory':
                        canonical_accessory_family.append(fam.name)
                        canonical_accessory_family += fam.exchangeables

                canonical_fu = list(canonical_system.model.func_units)[0]
                write_sys_to_organism(canonical_system, organism, projection, canonical_mandatory_family,
                                      canonical_accessory_family, canonical_fu)
    projection_df = pd.DataFrame(projection).drop_duplicates()
    return projection_df


def write_systems_projection(name: str, pangenome: Pangenome, output: Path, source: str, threads: int = 1,
                             force: bool = False, disable_bar: bool = False):
    """ Write all systems in pangenomes and project on organisms

    :param pangenome: Pangenome to project
    :param output: Path to output directory
    :param threads:  Number of threads available
    :param force: Force to write into the output directory
    :param disable_bar: Allow to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.getLogger().info(f'Write system projection for source : {source}')
        with tqdm(total=pangenome.number_of_systems(), unit='system', disable=disable_bar) as progress:
            futures = []
            for system in pangenome.get_system_by_source(source):
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
                                  "annotation name", "gene", "start", "stop", "strand", "is_fragment"]
    systems_projection.sort_values(by=["system number", "system name", "organism", "start", "stop"],
                                   ascending=[True, True, True, True, True], inplace=True)
    mkdir(f"{output}/{name}/projection_{source}_{name}", force=force)
    for organism_name in systems_projection["organism"].unique():
        org_df = systems_projection.loc[systems_projection["organism"] == organism_name]
        org_df = org_df.iloc[:, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10]]
        org_df.sort_values(by=["system number", "system name", "start", "stop"],
                           ascending=[True, True, True, True], inplace=True)
        org_df.to_csv(f"{output}/{name}/projection_{source}_{name}/{organism_name}.tsv", sep="\t", index=False)
    systems_projection.to_csv(f"{output}/{name}/systems_{source}_{name}.tsv", sep="\t", index=False)
    systems_projection.insert(0, "pangenome name", name)
    return systems_projection

def pan_presence_percentage(name: str, systems_projection):
    systems_projection_filter = systems_projection[['system name', 'organism']]
    system_names = list(systems_projection_filter['system name'].unique())
    system_names.sort(key=str.casefold)
    pan_count = len(system_names)*[None]
    for i, system in enumerate(system_names):
        data_system = systems_projection_filter[(systems_projection_filter['system name'] == system)]
        data_system_drop = data_system.drop_duplicates(subset=['system name', 'organism'])
        pan_count[i] = data_system_drop.shape[0]
    df = pd.DataFrame(list(zip(system_names, pan_count)), columns=["system name", "nb of presence"])
    df.insert(0, "pangenome name", name)
    return df

def pangenome_number_system(name: str, systems_projection):
    systems_projection_filter = systems_projection[['system number', 'system name', 'organism']]
    system_names = list(systems_projection['system name'].unique())
    system_names.sort(key=str.casefold)
    count_system_id = [None]*len(system_names)
    count_system_id_org = [None] * len(system_names)
    for i, system in enumerate(system_names):
        systems_projection_filter_2 = systems_projection_filter[systems_projection_filter['system name'] == system]
        systems_projection_id = systems_projection_filter_2.drop_duplicates(subset=['system number', 'system name'],
                                                                            keep="first", inplace=False)
        systems_projection_id_org = systems_projection_filter_2.drop_duplicates(subset=['system number', 'system name', 'organism'],
                                                                                keep="first", inplace=False)
        count_system_id[i] = systems_projection_id.shape[0]
        count_system_id_org[i] = systems_projection_id_org.shape[0]
    df_id = pd.DataFrame(list(zip(system_names, count_system_id)), columns=["system name", "Number of IDs"])
    df_id.insert(0, "pangenome name", name)
    df_id_org = pd.DataFrame(list(zip(system_names, count_system_id_org)), columns=["system name", "Number of IDs per org"])
    df_id_org.insert(0, "pangenome name", name)
    return df_id, df_id_org


def systems_to_features(systems_proj, name, pangenome: Pangenome, output: Path, source: str, rgp: bool = False,
                        modules: bool = False, spots: bool = False, threads: int = 1, disable_bar: bool = False):
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=pangenome.number_of_systems(), unit='system', disable=disable_bar) as progress:
            futures = []
            for system in pangenome.systems:
                future = executor.submit(system_to_features, system, pangenome.regions, pangenome.modules, rgp, modules)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
                #for canonical in system.canonical:
                    #future_rgp = executor.submit(system_to_rgp, canonical, pangenome.regions)
                    #future_rgp.add_done_callback(lambda p: progress.update())
                    #futures_rgp.append(future_rgp)

            result_list = []
            dict_system_rgp = {}
            dict_system_completion_rgp = {}
            dict_system_modules = {}
            for future in futures:
                res = future.result()
                result_list.append(res)
                if len(res[1]) > 0:
                    dict_system_rgp[res[0]] = res[1]
                    dict_system_completion_rgp[res[0]] = res[2]
                if len(res[3]) > 0:
                    dict_system_modules[res[0]] = res[3]

        dict_system_organism_with_rgp = {}
        dict_system_organism_with_module = {}
        for result in result_list:
            if result[1] != [] or result[3] != []:
                rgp_org = []
                mod_org = []
                for result_rgp in result[1]:
                    rgp_org.append(result_rgp.organism.name)
                dict_system_organism_with_rgp[result[0]] = rgp_org

                for result_module in result[3]:
                    mod_org.extend(x for x in list(map(lambda x: x.name, result_module.organisms)) if x not in mod_org)
                dict_system_organism_with_module[result[0]] = mod_org

        dict_region_spot = {}
        for spot in pangenome.spots:
            for region in spot.regions:
                dict_region_spot[region.name] = spot.ID

    df_features = write_systems_to_features(systems_proj, name, pangenome, output, source, dict_system_rgp, dict_system_completion_rgp,
                              dict_system_modules, dict_region_spot, dict_system_organism_with_rgp, dict_system_organism_with_module,
                                            rgp, modules, spots)
    return df_features

def write_systems_to_features(systems_proj, name, pangenome: Pangenome, output: Path, source: str,  dict_system_rgp: Dict[System, Region],
                              dict_system_completion_rgp: Dict[System, Completion], dict_system_modules: Dict[System, Module],
                              dict_region_spots: Dict[Region, Spot], dict_system_organism_with_rgp: Dict[System, Organism],
                              dict_system_organism_with_module: Dict[System, Organism], rgp: bool = False, modules: bool = False, spots: bool = False):
    list_to_df = []
    df_collection = collections.namedtuple("system_rgp_module_spot",
                                           ['system_ID', 'system_name', 'mod_organism', 'module', 'rgp_organism', 'spot', 'rgp', 'completion_rgp'])
    systems_ID = list(OrderedDict.fromkeys(list(dict_system_rgp) + list(dict_system_modules)))
    systems_proj_filter = systems_proj[['system number', 'organism']]

    for system_ID in systems_ID:
        system = pangenome.get_system(system_ID)
        organism_list_of_rgp = dict_system_organism_with_rgp.get(system_ID, None)
        organism_list_of_mod = dict_system_organism_with_module.get(system_ID, None)
        systems_proj_filter_2 = systems_proj_filter[systems_proj_filter['system number'] == system_ID]
        organism_all = list(systems_proj_filter_2['organism'].unique())
        rgp_list = dict_system_rgp.get(system_ID, None)
        completion_rgp = dict_system_completion_rgp.get(system_ID, None)
        module_list = dict_system_modules.get(system_ID, None)

        spot_list = []
        if rgp_list is not None:
            for rgp_element in rgp_list:
                spot_list.append(dict_region_spots[rgp_element.name]) if rgp_element.name in dict_region_spots else spot_list.append("no_spot")
            rgp_list = list(map(lambda x: x.name, rgp_list))
        else:
            spot_list = None

        module_presence = None
        mod_org = None
        if module_list is not None:
            module_presence = []
            mod_org = []
            for module in module_list:
                for organism in organism_all:
                    if organism in list(map(lambda x: x.name, module.organisms)):
                        mod_org.append(organism)
                        module_presence.append(module.ID)

        list_to_df.append(df_collection(system_ID, system.name, mod_org, module_presence, organism_list_of_rgp, spot_list, rgp_list, completion_rgp))

    if rgp is True and modules is True and spots is True:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_modules_spots_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'mod_organism', 'module', 'rgp_organism', 'spot', 'rgp', 'completion_rgp']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if rgp is True and modules is True and spots is False:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_modules_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'rgp_organism', 'rgp', 'completion_rgp', 'mod_organism', 'module']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if rgp is True and modules is False and spots is True:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_spots_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name','rgp_organism', 'spot', 'rgp', 'completion_rgp']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if rgp is False and modules is True and spots is True:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_modules_spots_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'mod_organism', 'module', 'rgp_organism', 'spot']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if rgp is True and modules is False and spots is False:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'rgp_organism', 'rgp', 'completion_rgp']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if rgp is False and modules is True and spots is False:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_modules_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'mod_organism', 'module']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    return df


def system_to_features(system: System, regions: Set[Region], modules: Set[Module], rgp, module):

    fu = list(system.model.func_units)[0]
    mandatory_family = []
    accessory_family = []
    for fam in system.model.families:
        if fam.presence == 'mandatory':
            mandatory_family.append(fam.name)
            mandatory_family += fam.exchangeables
        if fam.presence == 'accessory':
            accessory_family.append(fam.name)
            accessory_family += fam.exchangeables

    rgp_present = []
    completion_list = []
    if rgp :
        for rgp in regions:
            completion = len(list(system.gene_families & rgp.families)) / len(system.gene_families)
            if completion != 0 :
                find_projection = write_sys_in_feature(system, rgp, mandatory_family, accessory_family, fu)
                if find_projection is True:
                    rgp_present.append(rgp)
                    completion_list.append(completion)

    modules_present = []
    if module :
        for module in modules:
            completion = len(list(system.gene_families & module.families)) / len(system.gene_families)
            if completion != 0 :
                find_projection = write_sys_in_feature(system, module, mandatory_family, accessory_family, fu)
                if find_projection is True:
                    module.add_system(system)
                    modules_present.append(module)

    return system.ID, rgp_present, completion_list, modules_present


def write_sys_in_feature(system, feature, mandatory_family, accessory_family, fu) -> bool:
    """ Project system to rgp

        :param system: Detected system
        :param rgp: rgp in system
        :param projection: projection list output
        :param mandatory_family: list of mandatory families in system
        :param accessory_family: list of accessory family in system
        :param fu: functional unit

        :return: True if a projection occurs in rgp
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

    select_families = system.gene_families.intersection(feature.families)
    mandatory_index = {mandatory: index for index, mandatory in enumerate(mandatory_family, start=1)}
    accessory_index = {accessory: index for index, accessory in enumerate(accessory_family,
                                                                          start=1 + len(mandatory_family))}
    nb_annotation = len(mandatory_family) + len(accessory_family)
    bool_projection = False
    for comb_gf in itertools.permutations(select_families, fu.min_total):
        mandatory_matrix = np.zeros((len(comb_gf), len(mandatory_family)))
        all_matrix = np.concatenate((mandatory_matrix, np.zeros((len(comb_gf), len(accessory_family)))), axis=1)
        coordinate_dict = {}
        for index_gf, gf in enumerate(comb_gf):
            for index_annot, annot in enumerate(mandatory_family + accessory_family, start=1):
                coordinate_dict[index_gf * nb_annotation + index_annot] = (gf, annot)
            if gf.get_source(system.source) is not None:
                for annotation in [annot.get("protein_name") for annot in gf.get_source(system.source)]:
                    if annotation in mandatory_family:
                        mandatory_matrix[index_gf][mandatory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                      mandatory_index[annotation]
                        all_matrix[index_gf][mandatory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                mandatory_index[annotation]
                    elif annotation in accessory_family:
                        all_matrix[index_gf][accessory_index[annotation] - 1] = index_gf * nb_annotation + \
                                                                                accessory_index[annotation]
        if len(search_diagonal(mandatory_matrix, fu.min_mandatory)) > 0:
            search_diag = search_diagonal(all_matrix)
            if len(search_diag) > 0:
                bool_projection = True
    return bool_projection


def write_flat_files(pan_to_path, output: Path, annotation: bool = False, systems: bool = False, hmm: bool = False,
                     systems_asso: List[str] = None, proksee: List[str] = None, proksee_template: Path = None,
                     organisms_list: List[str] = None, threads: int = 1, force: bool = False, disable_bar: bool = False, **kwargs):
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
    need_rgp = False
    need_modules = False
    need_gene_sequences = False
    need_metadata = False
    need_systems = False
    rgp = False
    modules = False
    spots = False

    if annotation:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            logging.getLogger().info('Write annotation')
            with tqdm(total=len(pan_to_path.keys()), unit='pangenome', disable=disable_bar) as progress:
                need_families = True
                need_metadata = True
                for pangenome_name, pangenome_info in pan_to_path.items():
                    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
                    pangenome.add_file(pangenome_info["path"])
                    kwargs['sources'] = kwargs['sources'] if kwargs['sources'] is not None else pangenome.status['metasources']["families"]
                    check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                         need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                         need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                                         need_modules=need_modules,
                                         need_metadata=need_metadata, need_systems=need_systems,
                                         models=kwargs["models"], sources=kwargs["sources"], metatype="families",
                                         disable_bar=disable_bar)
                    future = executor.submit(write_annotations_to_families, pangenome, output, sources=kwargs["sources"], disable_bar=disable_bar)
                    future.add_done_callback(lambda p: progress.update())
            logging.getLogger().info(f"Annotation has been written in {output}/{pangenome_name}/families_annotations.tsv")

    if hmm:
        need_families = True
        need_annotations = True
        for pangenome_name, pangenome_info in pan_to_path.items():
            pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
            pangenome.add_file(pangenome_info["path"])
            check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                 need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                 need_spots=need_spots, need_gene_sequences=need_gene_sequences,
                                 need_modules=need_modules,
                                 need_metadata=need_metadata, need_systems=need_systems,
                                 models=kwargs["models"], sources=kwargs["sources"], metatype="families",
                                 disable_bar=disable_bar)
            msa_path = kwargs.get('msa_path', None)
            msa_format = kwargs.get('msa_format', 'afa')
            threads = kwargs.get('threads', 1)
            profile_gfs(pangenome, msa_path, msa_format, threads, disable_bar)
            write_hmm_profile(pangenome, output, threads, disable_bar)


    if systems or systems_asso:
        need_annotations = True
        need_families = True
        need_systems = True
        need_metadata = True

        global_df = pd.DataFrame()
        global_percentage = pd.DataFrame()
        global_id = pd.DataFrame()
        global_id_org = pd.DataFrame()

        for pangenome_name, pangenome_info in pan_to_path.items():
            pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
            pangenome.add_file(pangenome_info["path"])
            if systems :
                check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                    need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                    need_spots=need_spots, need_gene_sequences=need_gene_sequences, need_modules=need_modules,
                                    need_metadata=need_metadata, need_systems=need_systems, models=kwargs["models"],
                                    sources=kwargs["sources"], metatype="families", disable_bar=disable_bar)
            if systems_asso:
                need_rgp = True
                need_modules = True
                need_spots = True
                check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                    need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                    need_spots=need_spots, need_gene_sequences=need_gene_sequences, need_modules=need_modules,
                                    need_metadata=need_metadata, need_systems=need_systems, models=kwargs["models"],
                                    sources=kwargs["sources"], metatype="families", disable_bar=disable_bar)
            logging.getLogger().info(f"Begin write systems projection for {pangenome_name}")
            for source in kwargs["sources"]:
                systems_proj = write_systems_projection(name=pangenome_name, pangenome=pangenome, output=output, source=source,
                                                        threads=threads, force=force, disable_bar=disable_bar)
                global_df = pd.concat([global_df, systems_proj])
                per_pan_heatmap(pangenome_name, systems_proj, output)
                logging.getLogger().info(f"Heatmap figure created for {pangenome_name}")

                presence_ratio = pan_presence_percentage(pangenome_name, systems_proj)
                global_percentage = pd.concat([global_percentage, presence_ratio])

                dataframe_ID, dataframe_ID_org = pangenome_number_system(pangenome_name, systems_proj)
                global_id = pd.concat([global_id, dataframe_ID])
                global_id_org = pd.concat([global_id_org, dataframe_ID_org])
                logging.getLogger().info(f"Projection written for {pangenome_name}")

                hbar_ID_type(pangenome_name, dataframe_ID, dataframe_ID_org, output)

                if systems_asso:
                    if systems_asso in ["rgp", "rgp-modules", "rgp-spots", "all"]:
                        rgp = True
                    if systems_asso in ["modules", "rgp-modules", "modules-spots", "all"]:
                        modules = True
                    if systems_asso in ["rgp-spots", "all"]:
                        spots = True
                    if systems_asso in ["modules-spots", "all"]:
                        spots = True
                    logging.getLogger().info("Begin write systems with features projection")
                    df_features = systems_to_features(systems_proj, name=pangenome_name, pangenome=pangenome, output=output, source=source,
                                        rgp=rgp, modules=modules, spots=spots, threads=threads, disable_bar=disable_bar)
                    logging.getLogger().info("Projection with features written")
                
                    upsetplot(pangenome_name, systems_proj, df_features, output)

        heatmap(global_df, output)
        logging.getLogger().info(f"Global heatmap figure created")
        global_percentage.to_csv(f"{output}/systems_presence_ratio.tsv", sep="\t", index=False)
        global_df.to_csv(f"{output}/global_systems.tsv", sep="\t", index=False)
        global_id.to_csv(f"{output}/global_systems_number_id.tsv", sep="\t", index=False)
        global_id_org.to_csv(f"{output}/global_systems_number_id_org.tsv", sep="\t", index=False)


    if proksee is not None:
        need_annotations = True
        need_gene_sequences = True
        need_families = True
        need_partitions = True
        need_metadata = True
        if "rgp" in proksee or "all" in proksee:
            need_rgp = True
        if "spots" in proksee or "all" in proksee:
            need_rgp = True
            need_spots = True
        if "modules" in proksee or "all" in proksee:
            need_modules = True
        if "annotations" in proksee or "all" in proksee:
            need_annotations = True
        if "systems" in proksee or "all" in proksee:
            need_systems = True
        for pangenome_name, pangenome_info in pan_to_path.items():
            pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
            pangenome.add_file(pangenome_info["path"])
            check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                                 need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_rgp,
                                 need_spots=need_spots, need_gene_sequences=need_gene_sequences, need_modules=need_modules,
                                 need_metadata=need_metadata, need_systems=need_systems,
                                 models=kwargs["models"], sources=kwargs["sources"], metatype="families",
                                 disable_bar=disable_bar)
            write_proksee(pangenome=pangenome, output=output, features=proksee, template=proksee_template,
                          organisms_list=organisms_list, threads=threads, disable_bar=disable_bar)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    check_flat_parameters(args)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    mkdir(args.output, force=args.force)
    write_flat_files(pan_to_path, output=args.output, annotation=args.annotations, systems=args.systems,
                     systems_asso=args.systems_asso, models=args.models, sources=args.sources,
                     proksee=args.proksee, proksee_template=args.proksee_template, organisms_list=args.organisms,
                     hmm=args.hmm, msa_path=args.msa, msa_format=args.msa_format,
                     threads=args.threads, force=args.force, disable_bar=args.disable_prog_bar)



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
    optional.add_argument("--annotations", required=False, action="store_true",
                          help="Write all the annotations from families")
    optional.add_argument("--systems", required=False, action="store_true",
                          help="Write all the systems in pangenomes and project on genomes")
    optional.add_argument("--systems_asso", required=False, type=str, default=None,
                          choices=["all", "rgp-modules", "rgp-spots", "modules-spots", "modules", "rgp"],
                          help="Write association between systems and others pangenomes elements")
    optional.add_argument('--models', required=False, type=Path, default=None, nargs="+",
                          help="Path to model directory")
    optional.add_argument("--sources", required=False, type=str, nargs="+", default=None,
                          help='Name of the annotation source where panorama as to select in pangenomes')
    optional.add_argument("--proksee", required=False, type=str, default=None, nargs='+',
                          choices=["all", "base", "modules", "rgp", "spots", "annotations", "systems"])
    optional.add_argument("--proksee_template", required=False, type=Path, default=None, nargs='?')
    optional.add_argument("--organisms", required=False, type=str, default=None, nargs='+')
    optional.add_argument("--hmm", required=False, action="store_true",
                          help="Write an hmm for each gene families in pangenomes")
    optional.add_argument("--msa", required=False, type=Path, default=None,
                          help="To create a HMM profile for families, you can give a msa of each gene in families."
                               "This msa could be get from ppanggolin (See ppanggolin msa). "
                               "If no msa provide Panorama will launch one.")
    optional.add_argument("--msa_format", required=False, type=str, default="afa",
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
