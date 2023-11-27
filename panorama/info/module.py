#!/usr/bin/env python3
# coding:utf-8

# installed libraries
from pathlib import Path

import tables
from ppanggolin.pangenome import Pangenome
import pandas as pd
from bokeh.plotting import save
from bokeh.models import ColumnDataSource, DataTable, TableColumn


def get_module_info(pangenome: Pangenome) -> dict:
    """ Get module information from pangenome file

    :param pangenome: Pangenome with module computed

    :return: Dictionary with module information
    """
    h5f = tables.open_file(pangenome.file, "r")
    if "/info" in h5f:
        info_group = h5f.root.info
        module_info_list = ['numberOfModules', 'numberOfFamiliesInModules', 'PersistentSpecInModules',
                            'ShellSpecInModules', 'CloudSpecInModules', 'StatOfFamiliesInModules']
        if all(info in info_group._v_attrs._f_list() for info in module_info_list):
            module_dic = {'Number of Module': info_group._v_attrs['numberOfModules'],
                          'Number of families in modules': info_group._v_attrs['numberOfFamiliesInModules'],
                          "Percent of persistent families":  info_group._v_attrs['PersistentSpecInModules']['percent'],
                          "Percent of shell families": info_group._v_attrs['ShellSpecInModules']['percent'],
                          "Percent of cloud families": info_group._v_attrs['CloudSpecInModules']['percent'],
                          "Minimum Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['min'],
                          "Maximum Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['max'],
                          "Sd Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['sd'],
                          "Mean Families per Modules": info_group._v_attrs['StatOfFamiliesInModules']['mean'],
                          }
            return module_dic
        else:
            raise Exception(f"No information about modules find in {pangenome.file} "
                            f"Please use ppanggolin module -p {pangenome.file} to compute module and "
                            f"ppanggolin metrics --info_modules {pangenome.file} to compute information about module")
    else:
        raise Exception(f"No information find in {pangenome.file}")


def export_modules(modules_dict: dict, output: Path):
    """ Export modules information

    :param modules_dict: Content information
    :param output: Path to output directory
    """
    df = pd.DataFrame.from_dict(modules_dict, orient='index').T
    df.to_csv(path_or_buf=output/"modules_info.tsv")
    df = df.reset_index().rename(columns={"index": "Information"})
    source = ColumnDataSource(df)
    columns = [TableColumn(field=col, title=col) for col in df.columns]
    dt = DataTable(source=source, columns=columns, index_position=None,
                   sizing_mode='stretch_both')
    save(dt, filename=output/"modules_info.html", title="Pangenome modules Information",
         resources="cdn")
