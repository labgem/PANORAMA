#!/usr/bin/env python3
# coding:utf-8

# installed libraries
import tables
from ppanggolin.pangenome import Pangenome
import pandas as pd


def get_module_info(pangenome: Pangenome) -> dict:
    h5f = tables.open_file(pangenome.file, "r")
    if "/info" in h5f:
        info_group = h5f.root.info
        module_dic = {}
        if 'numberOfModules' in info_group._v_attrs._f_list():
            module_dic['Number of Module'] = info_group._v_attrs['numberOfModules']
            module_dic['Families'] = {'Number of families in mofules': info_group._v_attrs['numberOfFamiliesInModules'],
                                      'Minimum of families': info_group._v_attrs['StatOfFamiliesInModules']['min'],
                                      "Maximum of families": info_group._v_attrs['StatOfFamiliesInModules']['max'],
                                      "sd": info_group._v_attrs['StatOfFamiliesInModules']['sd'],
                                      "mean": info_group._v_attrs['StatOfFamiliesInModules']['mean']
                                      }
            module_dic['Families sheel specific'] = {
                'Percent of families': info_group._v_attrs['ShellSpecInModules']['percent'],
                'sd': info_group._v_attrs['ShellSpecInModules']['sd'],
                'mean': info_group._v_attrs['ShellSpecInModules']['mean']
            }
            module_dic['Families cloud specific'] = {
                'Percent of families': info_group._v_attrs['CloudSpecInModules']['percent'],
                'sd': info_group._v_attrs['CloudSpecInModules']['sd'],
                'mean': info_group._v_attrs['CloudSpecInModules']['mean']
            }
            return module_dic
        else:
            raise Exception(f"No information about modules find in {pangenome.file}"
                            f"Please use ppanggolin module -p {pangenome.file} to compute module")
    else:
        raise Exception(f"No information find in {pangenome.file}")