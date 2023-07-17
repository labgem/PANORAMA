#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import os.path
import ast
from pathlib import Path
from typing import List

# installed libraries
import pandas as pd
from pandas import DataFrame
from collections import Counter
# local libraries
from panorama.pangenomes import Pangenome


def write_borders_spot(name: str, pangenome : Pangenome) -> DataFrame:
    """ Extract borders of pangenome spots

         :param name: Name of the pangenome
         :param pangenome: pangenome object

         :return: Dataframe with spots borders
         """

    df_borders = pd.DataFrame(columns=['pangenome_name', 'Spot', 'number_org', 'borders'])
    for spot in pangenome.spots:
        borders = spot.borders(set_size=pangenome.parameters["spots"]["set_size"], multigenics=pangenome.get_multigenics(dup_margin=0.05))
        for c, border in borders:
            famstring = ",".join([fam.name for fam in border[0]]) + "," + ",".join([fam.name for fam in border[1]])
            new_row = {'pangenome_name': name, 'Spot': spot.ID, 'number_org': c, 'borders': famstring}
            df_borders.loc[len(df_borders)] = new_row
    idx = df_borders.groupby('Spot')['number_org'].idxmax() # keep spot and border with highest number of genomes having the border
    df_borders = df_borders.loc[idx]
    df_borders = df_borders.reset_index(drop=True)

    return df_borders


def identical_spot(df_borders_global: DataFrame, df_spot_global: DataFrame, df_align: DataFrame, output: Path, threshold: int = 4):
    """ Detect identical spots between pangenomes

         :param df_borders_global: Spot borders for all pangenomes
         :param df_spot_global: Systems in each spot for all pangenomes
         :param df_align: Alignment results of families from pangenomes
         :param output: Path to output directory
         :param threshold: Number of families conserved between pangenomes
         """

    names = df_borders_global['pangenome_name'].unique().tolist()

# get a dict where keys are families and values are similar families in others pangenomes
    dict_similar_fam = {}
    for key, value in zip(df_align['query'], df_align['target']):
        if key in dict_similar_fam:
            dict_similar_fam[key].append(value)
        else:
            dict_similar_fam[key] = [value]

    conserved_spot_global = pd.DataFrame(columns=names)
    for index, name in enumerate(names):
        df_borders = df_borders_global[df_borders_global['pangenome_name'] == name]
        definition_spot = df_borders.set_index('Spot')['borders'].to_dict() # create dict with spot as key and border as value
        for key, value in definition_spot.items():
            borders = value.split(',')
            border_count_q = Counter(borders) # count number of duplicates and triplicates families in borders (Query)
            dupli_borders_q = [item for item, count in border_count_q.items() if count == 2]
            tripli_borders_q = [item for item, count in border_count_q.items() if count == 3]
            score_df = df_borders_global.drop(['number_org', 'borders'], axis=1) # create dataframe to assign a score if similar families are found
            score_df['score'] = 0
            for fam in borders:
                similar = dict_similar_fam[fam]
                similar_df = df_borders_global[df_borders_global['borders'].str.contains('|'.join(similar), case=False)]
                for i, row in similar_df.iterrows():
                    string_name = row['pangenome_name']
                    string_spot = row['Spot']
                    string_borders = row['borders']
                    string_borders = string_borders.split(',')
                    border_count_t = Counter(string_borders) # count number of duplicates and triplicates families in borders (Target)
                    dupli_borders_t = [item for item, count in border_count_t.items() if count == 2]
                    tripli_borders_t = [item for item, count in border_count_t.items() if count == 3]

                    point = 1 #add 1 to score by default

                    # if a duplicate in query and no duplicate in target
                    if fam in dupli_borders_q and all(elem not in dupli_borders_t for elem in similar) :
                        point = 0.5		# add 1 point rather than 2 (0.5 by duplicate value)

                    # if a triplicate in query and no duplicate or triplicate in target
                    if fam in tripli_borders_q and all(elem not in dupli_borders_t for elem in similar) and all(elem not in tripli_borders_t for elem in similar):
                        point = 0.34		# add 1.02 point rather than 3 (0.34 by duplicate value)

                    # if a triplicate in query and a duplicate in target
                    if fam in tripli_borders_q and any(elem in dupli_borders_t for elem in similar):
                        point = 0.67		# add 2.01 point rather than three (0.67 by duplicate value)

                    score_df.loc[(score_df['pangenome_name'] == string_name) & (score_df['Spot'] == string_spot), 'score'] += point

            good_score_df = score_df[score_df['score'] > threshold-1] # Keep spots passing threshold requirement (Default : 4)
            good_score_df['Spot'] = good_score_df['Spot'].astype(str)
            good_score_df['score'] = good_score_df['score'].astype(str)
            fused_score_df = good_score_df.groupby('pangenome_name')['Spot'].apply(lambda x: list(map(str, x))).reset_index() # join spot_id column values of each pangenome as one list

            indicators = good_score_df['pangenome_name'].unique()
            conserved_spot = pd.DataFrame(columns=names)
            for indicator in indicators:
                fused_score_df_f = fused_score_df[fused_score_df['pangenome_name'] == indicator]
                conserved_spot[indicator] = fused_score_df_f['Spot'].to_list()
            conserved_spot_global = pd.concat([conserved_spot_global, conserved_spot], ignore_index=True)

    df_spot_filter = filter_spot(names=names, conserved_spot_global=conserved_spot_global, df_spot_global=df_spot_global)
    df_link_system = regroup_spot(names=names, df=df_spot_filter)
    df_link_system.to_csv(output/f"conserved_spot_link_systems.tsv", sep='\t', index=False)

    df_all = regroup_spot(names=names, df=conserved_spot_global)
    df_all.to_csv(output/"all_conserved_spot.tsv", sep='\t', index=False)

    asso_conserved_spot_system(conserved_spot_regroup=df_link_system, df_spot_global=df_spot_global,
                               df_borders_global=df_borders_global, output=output)


def filter_spot(names: List[str], conserved_spot_global: DataFrame, df_spot_global: DataFrame) -> DataFrame:
    """ Remove spots not associated with a system

         :param names: Names of all pangenomes
         :param conserved_spot_global: Similar spots in pangenomes not filtered
         :param df_spot_global: Systems in each spots for all pangenomes

         :return: Dataframe with spots linked to a system
         """

    dict_pan_spot = {}
    for index, name in enumerate(names):
        df_spot = df_spot_global[df_spot_global['pangenome_name'] == name]
        spot = df_spot["Spot"].tolist()
        dict_pan_spot[name] = spot #dict with spots associated with system in each pangenome

    df_spot_involved = pd.DataFrame(columns=names)
    for i, row in conserved_spot_global.iterrows():
        boolean_system = False
        for name in names:
            if boolean_system is False:
                if isinstance(row[name], list):
                    spot_check = dict_pan_spot[name] #spots associated with system
                    spot_check2 = [str(elem) for elem in spot_check]
                    if any(elem in row[name] for elem in spot_check2): #spots to check if associated with system
                        boolean_system = True
                        df_spot_involved = df_spot_involved.append(row, ignore_index=True)

    return df_spot_involved

def regroup_spot(names: List[str], df: DataFrame) -> DataFrame:
    """ Regroup similar spots together

         :param names: Names of all pangenomes
         :param df: Dataframe with either all conserved spots or spots link to systems

         :return: Dataframe with similar spots in pangenomes regrouped
         """

    df = df.applymap(lambda x: x if isinstance(x, list) else [])  # empty list for NaN values

    # First loop to reunite conserved spots
    df2 = pd.DataFrame(columns=names)
    for i, row in df.iterrows():
        row_list = row.tolist()
        dict_1 = {key: [] for key in names}
        for name1, element in zip(names, row_list):
            filter_spot_inv = df[df[name1].apply(lambda x: any(item in x for item in element))]
            for name2 in names:
                liste = filter_spot_inv[name2].tolist()
                flat_list = [item for sublist in liste for item in sublist]
                flat_list = list(set(flat_list))
                flat_list = sorted(flat_list, key=int)
                dict_1[name2] += flat_list
            dict_1 = {key: list(set(values)) for key, values in dict_1.items()}
        df2 = df2.append(dict_1, ignore_index=True)

    # Second loop to reunite conserved spots
    df3 = pd.DataFrame(columns=names)
    for i, row in df2.iterrows():
        row_list = row.tolist()
        dict_2 = {key: [] for key in names}
        for name1, element in zip(names, row_list):
            filter_spot_inv = df2[df2[name1].apply(lambda x: any(item in x for item in element))]
            for name2 in names:
                liste = filter_spot_inv[name2].tolist()
                flat_list = [item for sublist in liste for item in sublist]
                flat_list = list(set(flat_list))
                flat_list = sorted(flat_list, key=int)
                dict_2[name2] += flat_list
            dict_2 = {key: list(set(values)) for key, values in dict_2.items()}
        df3 = df3.append(dict_2, ignore_index=True)

    df_str = df3.astype(str)
    df_str = df_str.drop_duplicates()

    return df_str


def asso_conserved_spot_system(conserved_spot_regroup: DataFrame, df_spot_global: DataFrame,
                               df_borders_global: DataFrame, output: Path):
    """ Write conserved spot and systems associated for each conserved spots detected

         :param conserved_spot_regroup: Similar spots in pangenomes regrouped
         :param df_spot_global: Systems in each spots for all pangenomes
         :param df_borders_global: Spot borders for all pangenomes
         :param output: Path to output directory
         """

    conserved_spot_regroup = conserved_spot_regroup.reset_index(drop=True)

    if not os.path.exists(f"{output}/projection/"):
        os.mkdir(f"{output}/projection/")

    df_spot_global["system_name"] = df_spot_global["system_name"].apply(repr)
    df_spot_global["organism"] = df_spot_global["organism"].apply(repr)
    df_spot_global["system_name"] = df_spot_global["system_name"].str.replace('"Counter', '')
    df_spot_global["system_name"] = df_spot_global["system_name"].str.replace(r'\(|\)', '')
    df_spot_global["system_name"] = df_spot_global["system_name"].str.replace(r'{|}', '')
    df_spot_global["system_name"] = df_spot_global["system_name"].str.replace('"', '')
    df_spot_global["organism"] = df_spot_global["organism"].str.replace('"Counter', '')
    df_spot_global["organism"] = df_spot_global["organism"].str.replace(r'\(|\)', '')
    df_spot_global["organism"] = df_spot_global["organism"].str.replace(r'{|}', '')
    df_spot_global["organism"] = df_spot_global["organism"].str.replace('"', '')

    list_ID = []
    list_pan = []
    list_spot = []
    list_sys = []
    list_org = []

    sum_spot = pd.DataFrame(columns=["ID", "Number_species", "Number_spots", "Number_systems", "Number_organisms"])
    for i, row1 in conserved_spot_regroup.iterrows():
        save_spot = pd.DataFrame()
        for col, row2 in row1.items():
            row2 = ast.literal_eval(row2)
            if row2:
                filter_global_spot = df_spot_global[df_spot_global['pangenome_name'] == col]
                df_borders_global_filter = df_borders_global[df_borders_global['pangenome_name'] == col]
                for row3 in row2:
                    df_borders_global_filter2 = df_borders_global_filter[df_borders_global_filter['Spot'] == int(row3)]
                    filter_global_spot["Spot"] = filter_global_spot["Spot"].astype(str)
                    filter_global_spot_2 = filter_global_spot[filter_global_spot['Spot'] == row3]
                    check_spot_present = ~filter_global_spot_2['Spot'].str.contains(row3).any() #check if spot is not present in spot list
                    if check_spot_present:
                        new_row = {'pangenome_name': col, 'Spot': row3, 'system_name': '-', 'organism': '-',
                                   'borders': df_borders_global_filter2['borders'].values[0],
                                   'nb_system': 0, 'nb_organism': 0, 'mean_nb_genes': '-', 'max_nb_genes': '-',
                                   'min_nb_genes': '-'}
                        save_spot = save_spot.append(new_row, ignore_index=True)
                    save_spot = pd.concat([save_spot, filter_global_spot_2], ignore_index=True)

        save_spot.to_csv(f"{output}/projection/conserved_spot_{i}.tsv", sep='\t', index=False)

        save_global_spot = save_spot.drop(columns=['system_name', 'organism'], axis=1)
        list_ID.append(f"conserved_spot_{i}")
        nb_pan = save_global_spot['pangenome_name'].nunique()
        list_pan.append(nb_pan)
        nb_spot = save_global_spot.shape[0]
        list_spot.append(nb_spot)
        system_sum = save_global_spot['nb_system'].sum()
        list_sys.append(system_sum)
        org_sum = save_global_spot['nb_organism'].sum()
        list_org.append(org_sum)

    sum_spot['ID'] = list_ID
    sum_spot['Number_species'] = list_pan
    sum_spot['Number_spots'] = list_spot
    sum_spot['Number_systems'] = list_sys
    sum_spot['Number_organisms'] = list_org

    sum_spot.to_csv(f"{output}/summary_conserved_spot.tsv", sep='\t', index=False)
