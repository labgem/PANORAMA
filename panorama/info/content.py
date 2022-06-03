#!/usr/bin/env python3
# coding:utf-8

# installed libraries
from pathlib import Path

import tables
from ppanggolin.pangenome import Pangenome
import pandas as pd
from bokeh.plotting import save
from bokeh.models import ColumnDataSource, DataTable, TableColumn


def get_content_info(pangenome: Pangenome) -> dict:
    """ Get content information from pangenome file

    :param pangenome:
    :return:
    """
    h5f = tables.open_file(pangenome.file, "r")
    if "/info" in h5f:
        info_group = h5f.root.info
        req_info_list = ['numberOfGenes', 'numberOfOrganisms', 'numberOfClusters', 'numberOfEdges', 'numberOfCloud',
                         'numberOfPersistent', 'persistentStats', 'numberOfShell', 'shellStats', 'cloudStats',
                         'numberOfPartitions', 'numberOfSubpartitions', 'StatOfFamiliesInModules']
        if all(info in info_group._v_attrs._f_list() for info in req_info_list):
            content_dic = {'Number of Genes': info_group._v_attrs['numberOfGenes'],
                           'Number of Genomes': info_group._v_attrs['numberOfOrganisms'],
                           'Number of Families': info_group._v_attrs['numberOfClusters'],
                           'Number of Edges': info_group._v_attrs['numberOfEdges'],
                           'Number of Persistent': info_group._v_attrs['numberOfPersistent'],
                           'Minimum of Persistent': info_group._v_attrs['persistentStats']['min'],
                           'Maximum of Persistent': info_group._v_attrs['persistentStats']['min'],
                           'Mean of Persistent': info_group._v_attrs['persistentStats']['mean'],
                           'Sd of Persistent': info_group._v_attrs['persistentStats']['sd'],
                           'Number of Shell': info_group._v_attrs['numberOfShell'],
                           'Minimum of Shell': info_group._v_attrs['shellStats']['min'],
                           'Maximum of Shell': info_group._v_attrs['shellStats']['min'],
                           'Mean of Shell': info_group._v_attrs['shellStats']['mean'],
                           'Sd of Shell': info_group._v_attrs['shellStats']['sd'],
                           'Number of Cloud': info_group._v_attrs['numberOfCloud'],
                           'Minimum of Cloud': info_group._v_attrs['cloudStats']['min'],
                           'Maximum of Cloud': info_group._v_attrs['cloudStats']['min'],
                           'Mean of Cloud': info_group._v_attrs['cloudStats']['mean'],
                           'Sd of Cloud': info_group._v_attrs['cloudStats']['sd']
                           }
        else:
            raise Exception(f"No information about modules find in {pangenome.file}"
                            f"Please use ppanggolin module -p {pangenome.file} to compute module and "
                            f"ppanggolin metrics --info_modules {pangenome.file} to compute information about module")
        if 'genome_fluidity' in info_group._v_attrs._f_list():
            for part, value in info_group._v_attrs['genomes_fluidity'].items():
                content_dic[f'Genome Fluidity {part}'] = value
        if 'family_fluidity' in info_group._v_attrs._f_list():
            for part, value in info_group._v_attrs['genomes_fluidity'].items():
                content_dic[f'Family Fluidity {part}'] = value
        if 'numberOfRGP' in info_group._v_attrs._f_list():
            content_dic["Number of RGPs"] = info_group._v_attrs['numberOfRGP']
        if 'numberOfSpots' in info_group._v_attrs._f_list():
            content_dic["Number of Spots"] = info_group._v_attrs['numberOfSpots']
        if 'numberOfModules' in info_group._v_attrs._f_list():
            content_dic["Number of Modules"] = info_group._v_attrs['numberOfModules']
            content_dic["Families in Modules"] = info_group._v_attrs['numberOfFamiliesInModules']
        return content_dic
    else:
        raise Exception(f"No information find in {pangenome.file}")


def export_content(content_dict: dict, output: Path):
    """ Export Pangenome content

    :param content_dict:
    :return:
    """
    df = pd.DataFrame.from_dict(content_dict, orient='index').reset_index().rename(columns={'index': "Pangenomes"})
    source = ColumnDataSource(df)
    columns = [TableColumn(field=col, title=col) for col in df.columns]
    dt = DataTable(source=source, columns=columns, index_position=None,
                   autosize_mode='fit_columns', sizing_mode='stretch_both')
    save(dt, filename=f"{output.absolute().as_posix()}/content_info.html", title="Pangenome content Information",
         resources="cdn")
