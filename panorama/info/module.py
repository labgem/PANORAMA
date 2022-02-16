#!/usr/bin/env python3
# coding:utf-8

# installed libraries
import tables
from ppanggolin.pangenome import Pangenome
import pandas as pd
from bokeh.plotting import show
from bokeh.models import ColumnDataSource, DataTable, TableColumn


def get_module_info(pangenome: Pangenome) -> dict:
    h5f = tables.open_file(pangenome.file, "r")
    if "/info" in h5f:
        info_group = h5f.root.info
        if 'numberOfModules' in info_group._v_attrs._f_list():
            module_dic = {'Number of Module': info_group._v_attrs['numberOfModules'],
                          'Number of families in modules': info_group._v_attrs['numberOfFamiliesInModules'],
                          'Minimum of families': info_group._v_attrs['StatOfFamiliesInModules']['min'],
                          "Maximum of families": info_group._v_attrs['StatOfFamiliesInModules']['max'],
                          "sd families": info_group._v_attrs['StatOfFamiliesInModules']['sd'],
                          "mean families": info_group._v_attrs['StatOfFamiliesInModules']['mean'],
                          'Percent of families in shell': info_group._v_attrs['ShellSpecInModules']['percent'],
                          'sd families in shell': info_group._v_attrs['ShellSpecInModules']['sd'],
                          'mean families in shell': info_group._v_attrs['ShellSpecInModules']['mean'],
                          'Percent of families in cloud': info_group._v_attrs['CloudSpecInModules']['percent'],
                          'sd families in cloud': info_group._v_attrs['CloudSpecInModules']['sd'],
                          'mean families in cloud': info_group._v_attrs['CloudSpecInModules']['mean']
                          }
            return module_dic
        else:
            raise Exception(f"No information about modules find in {pangenome.file}"
                            f"Please use ppanggolin module -p {pangenome.file} to compute module")
    else:
        raise Exception(f"No information find in {pangenome.file}")


def export_modules(modules_dict: dict):
    for pangenome, info in modules_dict.items():
        print("pika", info)
    df = pd.DataFrame.from_dict(modules_dict, orient='index').reset_index().rename(columns={'index': "Pangenomes"})
    source = ColumnDataSource(df)
    columns = [TableColumn(field=col, title=col) for col in df.columns]
    dt = DataTable(source=source, columns=columns, index_position=None)
    show(dt)
