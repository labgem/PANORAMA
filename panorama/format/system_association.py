#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
from copy import copy
import collections
import itertools
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple

# installed libraries
import numpy as np
import pandas as pd
from pandas import DataFrame
from tqdm import tqdm
from ppanggolin.region import Region, Spot
from ppanggolin.genome import Organism
# local libraries

from panorama.utils import mkdir
from panorama.system import System
from panorama.region import Module
from panorama.pangenomes import Pangenome

# For Quentin
from collections import OrderedDict, Counter



def write_sys_to_organism(system: System, organism, projection, mandatory_family, accessory_family, fu) -> bool:
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
                             force: bool = False, disable_bar: bool = False) -> DataFrame:
    """ Write all systems in pangenomes and project on organisms

    :param name: Name of the pangenome
    :param pangenome: Pangenome to project
    :param output: Path to output directory
    :param source:  Annotation source
    :param threads:  Number of threads available
    :param force: Force to write into the output directory
    :param disable_bar: Allow to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.info(f'Write system projection for source : {source}')
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
    mkdir(output/f"{name}/projection_{source}_{name}", force=force)
    for organism_name in systems_projection["organism"].unique():
        org_df = systems_projection.loc[systems_projection["organism"] == organism_name]
        org_df = org_df.iloc[:, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10]]
        org_df.sort_values(by=["system number", "system name", "start", "stop"],
                           ascending=[True, True, True, True], inplace=True)
        org_df.to_csv(f"{output}/{name}/projection_{source}_{name}/{organism_name}.tsv", sep="\t", index=False)
    systems_projection.to_csv(f"{output}/{name}/systems_{source}_{name}.tsv", sep="\t", index=False)
    systems_projection.insert(0, "pangenome name", name)
    return systems_projection

def pan_distribution_system(name: str, systems_projection: DataFrame) -> DataFrame:
    """ Distribution of all types of systems in pangenome organisms

    :param name: Name of the pangenome
    :param systems_projection: Systems in the pangenome
    """

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

def pan_number_system(name: str, systems_projection: DataFrame) -> Tuple[DataFrame, DataFrame]:
    """ Return two DataFrames : one with the number of ID per type of systems and one with total number of systems per type in the pangenome

    :param name: Name of the pangenome
    :param systems_projection: Systems in the pangenome
    """

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
    df_total = pd.DataFrame(list(zip(system_names, count_system_id_org)), columns=["system name", "Total number of systems"])
    df_total.insert(0, "pangenome name", name)
    return df_id, df_total


def systems_to_features(name: str, pangenome: Pangenome, systems_projection: DataFrame, output: Path, source: str, bool_rgp: bool = False,
                        bool_modules: bool = False, bool_spots: bool = False, threads: int = 1, disable_bar: bool = False) -> DataFrame:
    """ Associate systems to features (RGPs, spots, modules) for the pangenome

    :param name: Name of the pangenome
    :param pangenome: pangenome object
    :param systems_projection: Systems in the pangenome
    :param output: Path to output directory
    :param source: Annotation source
    :param bool_rgp: Boolean to associate or not RGPs
    :param bool_modules: Boolean to associate or not modules
    :param bool_spots: Boolean to associate or not spots
    :param threads: Number of threads available
    :param disable_bar: Allow to disable progress bar

    :return: Dataframe with the association between systems and features
    """

    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=pangenome.number_of_systems(), unit='system', disable=disable_bar) as progress:
            futures = []
            for system in pangenome.systems:
                future = executor.submit(system_to_features, system, pangenome.regions, pangenome.modules, bool_rgp, bool_modules)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
                #for canonical in system.canonical:
                    #future_rgp = executor.submit(system_to_rgp, canonical, pangenome.regions)
                    #future_rgp.add_done_callback(lambda p: progress.update())
                    #futures_rgp.append(future_rgp)

            result_list = []
            dict_system_rgp = {}
            dict_system_modules = {}
            for future in futures:
                res = future.result()
                result_list.append(res)
                if len(res[1]) > 0:
                    dict_system_rgp[res[0]] = res[1]
                if len(res[2]) > 0:
                    dict_system_modules[res[0]] = res[2]

        dict_system_organism_with_rgp = {}
        dict_system_organism_with_module = {}
        for result in result_list:
            if result[1] != [] or result[2] != []:
                rgp_org = []
                mod_org = []
                for result_rgp in result[1]:
                    rgp_org.append(result_rgp.organism.name)
                dict_system_organism_with_rgp[result[0]] = rgp_org

                for result_module in result[2]:
                    mod_org.extend(x for x in list(map(lambda x: x.name, result_module.organisms)) if x not in mod_org)
                dict_system_organism_with_module[result[0]] = mod_org

        dict_region_spot = {}
        for spot in pangenome.spots:
            for region in spot.regions:
                dict_region_spot[region.name] = spot.ID

    df_features = write_systems_to_features(name=name, pangenome=pangenome, systems_projection=systems_projection,
                                            output=output, source=source, dict_system_rgp=dict_system_rgp,
                                            dict_system_modules=dict_system_modules, dict_region_spots=dict_region_spot,
                                            dict_system_organism_with_rgp=dict_system_organism_with_rgp,
                                            bool_rgp=bool_rgp, bool_modules=bool_modules, bool_spots=bool_spots)
    return df_features

def write_systems_to_features(name: str, pangenome: Pangenome, systems_projection: DataFrame, output: Path, source: str,
                              dict_system_rgp: Dict[System, Region], dict_system_modules: Dict[System, Module],
                              dict_region_spots: Dict[Region, Spot], dict_system_organism_with_rgp: Dict[System, Organism],
                              bool_rgp: bool = False, bool_modules: bool = False, bool_spots: bool = False) -> DataFrame:
    """ Write function to associate systems to features (RGPs, spots, modules) for the pangenome

    :param name: Name of the pangenome
    :param pangenome: pangenome object
    :param systems_projection: Systems in the pangenome
    :param output: Path to output directory
    :param source: Annotation source
    :param dict_system_rgp: dictionary with system as key and a list of RGP as value
    :param dict_system_modules: dictionary with system as key and a list of module as value
    :param dict_region_spots: dictionary with region as key and spot as value
    :param dict_system_organism_with_rgp: dictionary with system as key and a list of organism names from RGP
    :param bool_rgp: Boolean to associate or not RGPs
    :param bool_modules: Boolean to associate or not modules
    :param bool_spots: Boolean to associate or not spots

    :return: Dataframe with the association between systems and features
    """

    list_to_df = []
    df_collection = collections.namedtuple("system_rgp_module_spot",
                                           ['system_ID', 'system_name', 'mod_organism', 'module', 'rgp_organism', 'spot', 'rgp'])
    systems_ID = list(OrderedDict.fromkeys(list(dict_system_rgp) + list(dict_system_modules)))
    systems_projection_filter = systems_projection[['system number', 'organism']]

    for system_ID in systems_ID:
        system = pangenome.get_system(system_ID)
        organism_list_of_rgp = dict_system_organism_with_rgp.get(system_ID, None)
        systems_projection_filter_2 = systems_projection_filter[systems_projection_filter['system number'] == system_ID]
        organism_all = list(systems_projection_filter_2['organism'].unique())
        rgp_list = dict_system_rgp.get(system_ID, None)
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

        list_to_df.append(df_collection(system_ID, system.name, mod_org, module_presence, organism_list_of_rgp, spot_list, rgp_list))

    if bool_rgp is True and bool_modules is True and bool_spots is True:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_modules_spots_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'mod_organism', 'module', 'rgp_organism', 'spot', 'rgp']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if bool_rgp is True and bool_modules is True and bool_spots is False:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_modules_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'rgp_organism', 'rgp', 'mod_organism', 'module']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if bool_rgp is True and bool_modules is False and bool_spots is True:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_spots_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name','rgp_organism', 'spot', 'rgp']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if bool_rgp is False and bool_modules is True and bool_spots is True:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_modules_spots_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'mod_organism', 'module', 'rgp_organism', 'spot']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if bool_rgp is True and bool_modules is False and bool_spots is False:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_rgp_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'rgp_organism', 'rgp']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    if bool_rgp is False and bool_modules is True and bool_spots is False:
        df = pd.DataFrame(list_to_df)
        outpath = Path(output / name / "systems_to_modules_{0}_{1}".format(source, name)).with_suffix(".tsv")
        header = ['system_ID', 'system_name', 'mod_organism', 'module']
        df.to_csv(path_or_buf=outpath, sep="\t", index=False, columns = header)

    return df


def system_to_features(system: System, regions: Set[Region], modules: Set[Module], bool_rgp: bool, bool_modules: bool) -> Tuple[System, List[Region], List[Module]]:
    """ Associate systems to features (RGPs, spots, modules) for the pangenome

    :param system: system from the pangenome
    :param regions: set of region
    :param modules: set of module object
    :param bool_rgp: Boolean to associate or not RGPs
    :param bool_modules: Boolean to associate or not modules

    :return: System and region associated
    """

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
    if bool_rgp :
        for rgp in regions:
            completion = len(list(system.gene_families & rgp.families)) / len(system.gene_families)
            if completion != 0 :
                find_projection = write_sys_in_feature(system=system, feature=rgp, mandatory_family=mandatory_family,
                                                       accessory_family=accessory_family, fu=fu)
                if find_projection is True:
                    rgp_present.append(rgp)

    modules_present = []
    if bool_modules :
        for module in modules:
            completion = len(list(system.gene_families & module.families)) / len(system.gene_families)
            if completion != 0 :
                find_projection = write_sys_in_feature(system=system, feature=module, mandatory_family=mandatory_family,
                                                       accessory_family=accessory_family, fu=fu)
                if find_projection is True:
                    module.add_system(system)
                    modules_present.append(module)

    return system.ID, rgp_present, modules_present


def write_sys_in_feature(system: System, feature: Region or Module, mandatory_family: List, accessory_family: List, fu: List) -> bool:
    """ Project system in RGPs or modules

        :param system: Detected system
        :param feature: feature (RGP or module)
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


def spot2sys(name: str, pangenome: Pangenome, system_to_feature: DataFrame, output: Path) -> Tuple[DataFrame, Dict]:
    """Write systems associated to spots of pangenome and extract spot content

    :param name: pangenome name
    :param pangenome: pangenome object
    :param system_to_feature: Dataframe with modules, spots and RGPs associated to each system ID of pangenome
    :param output: Path to output directory

    :return: Dataframe with systems in each spots
    """

    system_to_feature = system_to_feature.drop(columns=['mod_organism', 'module', 'rgp'])
    system_to_feature = system_to_feature.where(pd.notnull(system_to_feature), '[]')
    spot_set = set()
    for index, value in system_to_feature['spot'].iteritems():
        spot_set = set(value) | spot_set
    spot_set.remove('no_spot')
    spot_list = sorted(list(spot_set))

    dict_spot_ID_system = {}
    dict_spot_system = {}
    dict_spot_org = {}
    for spot in spot_list:
        system_to_feature_2 = system_to_feature[system_to_feature['spot'].apply(lambda x: spot in x)]
        for index1, row in system_to_feature_2.iterrows():
            system_ID = row['system_ID']
            system = row['system_name']
            row_spots = row['spot']
            organism = row['rgp_organism']
            for index2, row_spot in enumerate(row_spots):
                if row_spot == spot:
                    dict_spot_ID_system.setdefault(spot, []).append(system_ID)
                    dict_spot_system.setdefault(spot, []).append(system)
                    dict_spot_org.setdefault(spot, []).append(organism[index2])
    df_spot = pd.DataFrame({'Spot': list(dict_spot_system.keys()), 'system_name': list(dict_spot_system.values()), 'organism': list(dict_spot_org.values())})

    number_sys_list = []
    number_org_list = []
    system_count_list = []
    org_count_list = []
    for index, row in df_spot.iterrows():
        number_sys = row['system_name']
        number_sys_list.append(len(number_sys))
        number_org = list(set(row['organism']))
        number_org_list.append(len(number_org))
        system_count = Counter(row['system_name'])
        org_count = Counter(row['organism'])
        system_count_list.append(system_count)
        org_count_list.append(org_count)

    df_spot['system_name'] = system_count_list
    df_spot['organism'] = org_count_list
    df_spot['nb_system'] = number_sys_list
    df_spot['nb_organism'] = number_org_list

    df_spot["system_name"] = df_spot["system_name"].apply(repr)
    df_spot["organism"] = df_spot["organism"].apply(repr)
    df_spot["system_name"] = df_spot["system_name"].str.replace("Counter", "")
    df_spot["system_name"] = df_spot["system_name"].str.replace(r"\(|\)", "", regex=True)
    df_spot["system_name"] = df_spot["system_name"].str.replace(r"{|}", "", regex=True)
    df_spot["system_name"] = df_spot["system_name"].str.replace('"', "")
    df_spot["system_name"] = df_spot["system_name"].str.replace("'", "")
    df_spot["organism"] = df_spot["organism"].str.replace("Counter", "")
    df_spot["organism"] = df_spot["organism"].str.replace(r"\(|\)", "", regex=True)
    df_spot["organism"] = df_spot["organism"].str.replace(r"{|}", "", regex=True)
    df_spot["organism"] = df_spot["organism"].str.replace("'", "")
    df_spot["organism"] = df_spot["organism"].str.replace('"', '')

    # Extraction spot content

    df_spot = extract_spot_content(pangenome=pangenome, df_spot=df_spot)

    df_spot.to_csv(output / name / f"spot_to_system_{name}.tsv", sep='\t', header=["Spot", "system_name", "organism",
                                                                                   "nb_system", "nb_organism",
                                                                                   "Average_number_genes", "spot_content" ])

    df_spot.insert(0, 'pangenome_name', name)
    return df_spot, dict_spot_org


def extract_spot_content(pangenome: Pangenome, df_spot: DataFrame) -> DataFrame:
    """Add information about spot content (average number of genes and product)

    :param pangenome: pangenome object
    :param df_spot: Systems in each spot for the pangenome

    :return: Dataframe with spot content
    """
    spot_list = df_spot["Spot"].tolist()
    list_size_spot = []
    list_spot_content = []
    for s in pangenome.spots:
        list_size_region = []
        list_rgp_content = []
        if s.ID in spot_list:
            for region in s.regions:
                list_size_region.append(len(list(region.genes))) # append number of genes of each rgp
                genes_region = region.genes
                for gene in genes_region:
                    list_rgp_content.append(gene.product)
            list_size_spot.append(sum(list_size_region)/len(list_size_region)) # append average number of genes
            list_spot_content.append(Counter(list_rgp_content))

    df_spot['Average_number_genes'] = list_size_spot
    df_spot['spot_content'] = list_spot_content

    df_spot["Average_number_genes"] = df_spot["Average_number_genes"].apply(repr)
    df_spot["spot_content"] = df_spot["spot_content"].apply(repr)
    df_spot["Average_number_genes"] = df_spot["Average_number_genes"].str.replace("Counter", "")
    df_spot["Average_number_genes"] = df_spot["Average_number_genes"].str.replace(r"\(|\)", "", regex=True)
    df_spot["Average_number_genes"] = df_spot["Average_number_genes"].str.replace(r"{|}", "", regex=True)
    df_spot["Average_number_genes"] = df_spot["Average_number_genes"].str.replace('"', "")
    df_spot["Average_number_genes"] = df_spot["Average_number_genes"].str.replace("'", "")
    df_spot["spot_content"] = df_spot["spot_content"].str.replace("Counter", "")
    df_spot["spot_content"] = df_spot["spot_content"].str.replace(r"\(|\)", "", regex=True)
    df_spot["spot_content"] = df_spot["spot_content"].str.replace(r"{|}", "", regex=True)
    df_spot["spot_content"] = df_spot["spot_content"].str.replace("'", "")
    df_spot["spot_content"] = df_spot["spot_content"].str.replace('"', '')

    return df_spot

def mod2sys(name: str, system_to_feature: DataFrame, output: Path) -> DataFrame:
    """Write systems associated to modules of pangenome

    :param name: pangenome name
    :param system_to_feature: Dataframe with modules, spots and RGPs associated to each system ID of pangenome
    :param output: Path to output directory

    :return: Dataframe with systems in each modules
    """

    system_to_feature = system_to_feature.drop(columns=['system_ID', 'rgp_organism', 'spot', 'rgp'])
    system_to_feature = system_to_feature.applymap(lambda x: [] if x is None else x)
    module_set = set()
    for index, value in system_to_feature['module'].iteritems():
        module_set = set(value) | module_set
    module_list = sorted(list(module_set))

    dict_module_system = {}
    dict_module_org = {}
    for module in module_list:
        system_to_feature_2 = system_to_feature[system_to_feature['module'].apply(lambda x: module in x)]
        for index1, row in system_to_feature_2.iterrows():
            system = row['system_name']
            row_modules = row['module']
            organism = row['mod_organism']
            for index2, row_module in enumerate(row_modules):
                if row_module == module:
                    dict_module_system.setdefault(module, []).append(system)
                    dict_module_org.setdefault(module, []).append(organism[index2])
    df_module = pd.DataFrame({'Module': list(dict_module_system.keys()), 'system_name': list(dict_module_system.values()),
                              'organism': list(dict_module_org.values())})

    number_sys_list = []
    number_org_list = []
    system_count_list = []
    org_count_list = []
    for index, row in df_module.iterrows():
        number_sys = row['system_name']
        number_sys_list.append(len(number_sys))
        number_org = list(set(row['organism']))
        number_org_list.append(len(number_org))
        system_count = Counter(row['system_name'])
        org_count = Counter(row['organism'])
        system_count_list.append(system_count)
        org_count_list.append(org_count)

    df_module['system_name'] = system_count_list
    df_module['organism'] = org_count_list
    df_module['nb_system'] = number_sys_list
    df_module['nb_organism'] = number_org_list

    df_module["system_name"] = df_module["system_name"].apply(repr)
    df_module["organism"] = df_module["organism"].apply(repr)
    df_module["system_name"] = df_module["system_name"].str.replace('Counter', '')
    df_module["system_name"] = df_module["system_name"].str.replace(r'\(|\)', '', regex=True)
    df_module["system_name"] = df_module["system_name"].str.replace(r'{|}', '', regex=True)
    df_module["system_name"] = df_module["system_name"].str.replace('"', '')
    df_module["system_name"] = df_module["system_name"].str.replace("'", "")
    df_module["organism"] = df_module["organism"].str.replace('Counter', '')
    df_module["organism"] = df_module["organism"].str.replace(r'\(|\)', '', regex=True)
    df_module["organism"] = df_module["organism"].str.replace(r'{|}', '', regex=True)
    df_module["organism"] = df_module["organism"].str.replace('"', '')
    df_module["organism"] = df_module["organism"].str.replace("'", "")

    df_module.to_csv(output / name / f"module_to_system_{name}.tsv", sep='\t', header=["Module", "system_name", "organism", "nb_system", "nb_organism"])
    df_module.insert(0, 'pangenome_name', name)

    return df_module
