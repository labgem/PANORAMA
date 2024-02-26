#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from math import floor, ceil
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union
from copy import copy

# installed libraries
from tqdm import tqdm
import pandas as pd
import tables
from bokeh.models import Button, CheckboxGroup, ColumnDataSource, DataTable, Div, RadioButtonGroup, RangeSlider, TableColumn
from bokeh.layouts import column, row
from bokeh.io import curdoc
from bokeh.models.callbacks import CustomJS
from bokeh.embed import file_html
from ppanggolin.info.info import read_status, read_metadata_status
from ppanggolin.formats import read_info, get_pangenome_parameters

# local libraries
from panorama.utils import check_tsv_sanity


def get_info(pangenomes_path: Dict[str, Dict[str, Union[int, str]]], status: bool = False, content: bool = False,
             parameters: bool = False, metadata: bool = False, disable_bar: bool = False) -> dict:
    if not (status or content or parameters or metadata):
        status, content, parameters, metadata = (True, True, True, True)
    info_dict = defaultdict(dict)
    for pangenome_name, pan_info in tqdm(pangenomes_path.items(), unit="Pangenome", disable=disable_bar):
        h5f = tables.open_file(pan_info["path"], "r+")
        if status:
            info_dict['status'][pangenome_name] = read_status(h5f)["Status"]
        if content:
            info_dict['content'][pangenome_name] = read_info(h5f)["Content"]
        if parameters:
            step_to_parameters = get_pangenome_parameters(h5f)
            print(step_to_parameters)
            # info_dict['parameters'][pangenome_name] =
        if metadata:
            info_dict['metadata'][pangenome_name] = read_metadata_status(h5f)
        h5f.close()
    return info_dict


def export_status(status_dict: Dict[str, Union[bool, str]], output: Path):
    logging.getLogger("PANORAMA").debug("Exporting status")
    status_df = pd.DataFrame.from_dict(status_dict, orient='index')
    status_df = status_df.reset_index().rename(columns={"index": "Pangenome name"})
    source = ColumnDataSource(status_df)

    columns = [TableColumn(field=status, title=status) for status in status_df.columns]
    status_dt = DataTable(source=source, columns=columns, index_position=None,
                          width=1900, height = 30*status_df.shape[0])
    radio_buttons = []
    for column_name in source.column_names:
        if source.data[column_name].dtype == 'bool':
            radio_buttons.append(RadioButtonGroup(labels=['all', 'true', 'false'], active=0, name=column_name))

    filter_button = Button(label="Filter", button_type="primary")
    filter_button.js_on_event('button_click',
                              CustomJS(args=dict(source=source, radio_buttons=radio_buttons, original_source=copy(source.data)),
                                                 code=open(Path(__file__).parent/'filterData.js').read() +
                                                      "filterDataBool(source, radio_buttons);"))

    download_button = Button(label="Download", button_type="success")
    download_button.js_on_event('button_click',
                             CustomJS(args=dict(source=source, filename="pangenomes_status.tsv"),
                                      code=open(Path(__file__).parent/'download.js').read()))


    layout = column(status_dt,
                    row(Div(text='', width=160), row(*radio_buttons, spacing=29), Div(text='', width=20),
                        filter_button, download_button)
                    )

    curdoc().add_root(layout)

    html_content = file_html(layout, "cdn", title="Pangenome status Information")
    with open(output / 'status_info.html', "w") as file:
        file.write(html_content)

def unpack_content_dict(content_dict: Dict[str, Dict[str, Union[int, float, Dict[str, Union[int, float]]]]]):
    unpack_dict = defaultdict(dict)
    for name, content in content_dict.items():
        for key, value in content.items():
            if isinstance(value, dict):
                for k, v in value.items():
                    if isinstance(v, dict):
                        unpack_dict[name].update({f'{k}.{i}': j for i, j in v.items()})
                    elif k == "Number_of_modules" or re.match(r'.*_count', k):
                        unpack_dict[name][key] = v
                    else:
                        unpack_dict[name][k] = v
            else:
                unpack_dict[name][key] = value
    return unpack_dict

def export_content(content_dict: Dict[str, Dict[str, Union[int, float, Dict[str, Union[int, float]]]]], output: Path):
    column_names = {"min_genomes_frequency": "Genomes_frequency.min",
                    "max_genomes_frequency": "Genomes_frequency.max",
                    "mean_genomes_frequency": "Genomes_frequency.mean",
                    "sd_genomes_frequency": "Genomes_frequency.sd"}
    logging.getLogger("PANORAMA").debug("Exporting content")
    content_df = pd.DataFrame.from_dict(unpack_content_dict(content_dict), orient="index")
    content_df = content_df.reset_index().rename(columns={"index": "Pangenome name"})
    content_df = content_df.rename(columns=column_names)
    column_headers = content_df.columns.values.tolist()
    column_order = (column_headers[:3] + column_headers[6:8] + column_headers[9: 7 :-1] + column_headers[3:6] +
                    column_headers[10:13] + column_headers[17:20] + column_headers[13:17] + column_headers[20:])
    content_df = content_df.reindex(columns = column_order)
    source = ColumnDataSource(content_df)

    visible_index = [0, 1, 2, 7, 8, 9, 10, 11, 16, 17, 18]
    columns = [TableColumn(field=status, title=status, visible=True if index in visible_index else False)
               for index, status in enumerate(content_df.columns)]
    content_dt = DataTable(source=source, columns=columns, index_position=None,
                           width=1900, height = 32*content_df.shape[0])
    checkbox_group = CheckboxGroup(labels=source.column_names[1:], active=visible_index)
    checkbox_group.js_on_change('active',
                                CustomJS(args=dict(source=source, columns=columns, checkbox_group=checkbox_group),
                                          code=open(Path(__file__).parent / 'hide_show.js').read()))
    sliders = []
    for column_name in source.column_names[1:]:
        if source.data[column_name].dtype in (int, float):
            start = min(source.data[column_name]) if source.data[column_name].dtype == int else floor(min(source.data[column_name]))
            end = max(source.data[column_name]) if source.data[column_name].dtype == int else ceil(max(source.data[column_name]))
            step = 1 if source.data[column_name].dtype == int else .1
            if start == end:
                if start > 1:
                    start -= 1
                end += 1
            sliders.append(RangeSlider(start=start, end=end, value=(start, end), step = step, title=column_name))

    for idx, slider in enumerate(sliders):
        other_sliders = sliders[:idx] + sliders[idx+1:]
        slider.js_on_change('value_throttled',
                            CustomJS(args=dict(source=source, slider=slider, other_sliders=other_sliders,
                                               original_source=copy(source.data)),
                                                 code=open(Path(__file__).parent/'filterData.js').read() +
                                                      "filterDataSliders(source, slider, other_sliders);"))
    download_button = Button(label="Download", button_type="success")
    download_button.js_on_event('button_click',
                                CustomJS(args=dict(source=source, filename="pangenomes_content.tsv"),
                                         code=open(Path(__file__).parent / 'download.js').read()))

    layout = column(content_dt,
                    row(checkbox_group, column(*sliders[:6]), column(*sliders[6:12]),
                        column(*sliders[12:18]), column(*sliders[18:],
                                                        download_button),
                        spacing=20),
                    spacing=50)

    curdoc().add_root(layout)

    html_content = file_html(layout, "cdn", title="Pangenome content information")
    with open(output / 'content_info.html', "w") as file:
        file.write(html_content)



def export_info(info_dict: dict, output: Path):
    """ Export information to HTML file

    :param info_dict: Dictionary with information readable
    :param output: Path to output directory
    """
    for info, value in info_dict.items():
        if info == "status":
            export_status(value, output)
        elif info == "content":
            export_content(value, output)
        elif info == "parameters":
            raise NotImplementedError('Sorry I did not implement this for the time')
        elif info == "metadata":
            raise NotImplementedError('Sorry I did not implement this for the time')
        else:
            raise KeyError("The key to save information is not recognized")


def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    logging.getLogger("PANORAMA").debug("launch info command")
    pangenomes_to_path = check_tsv_sanity(args.pangenomes)
    info_dict = get_info(pangenomes_path=pangenomes_to_path, status=args.status, content=args.content,
                         parameters=args.parameters, metadata=args.metadata, disable_bar=args.disable_prog_bar)
    export_info(info_dict, args.output)
    logging.getLogger("PANORAMA").info("Done")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("info", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_info(parser)
    return parser


def parser_info(parser):
    """
    Parser for specific argument of info command

    :param parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help="A list of pangenome .h5 files")
    required.add_argument('-o', '--output', required=True, type=Path, nargs='?')
    optional = parser.add_argument_group(title="Information Display Options (default: all)")
    optional.add_argument("-a", "--parameters", required=False, action="store_true",
                          help="Create a file to compare parameters used or computed at "
                               "each step of pangenome generation for all pangenomes")
    optional.add_argument("-c", "--content", required=False, action="store_true",
                          help="Create a file to compare pangenome content between all pangenomes")
    optional.add_argument("-s", "--status", required=False, action="store_true",
                          help="Create a file to compare statuses of different elements between pangenomes")
    optional.add_argument("-m", "--metadata", required=False, action="store_true",
                          help="Create a file to summary the metadata saved in all the pangenome")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
