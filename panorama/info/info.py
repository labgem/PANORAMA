#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Union
from copy import copy

# installed libraries
from tqdm import tqdm
import pandas as pd
import tables
from bokeh.models import ColumnDataSource, DataTable, TableColumn, RadioButtonGroup, Button, Div
from bokeh.layouts import column, row
from bokeh.io import curdoc, output_file
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
    labels = []
    for column_name in source.column_names:
        if source.data[column_name].dtype == 'bool':
            radio_buttons.append(RadioButtonGroup(labels=['all', 'true', 'false'], active=0, name=column_name))
            labels.append(Div(text=f"<b>{column_name}</b>"))

    filter_button = Button(label="Filter", button_type="primary")
    filter_button.js_on_event('button_click',
                              CustomJS(args=dict(source=source, radio_buttons=radio_buttons, original_source=copy(source.data)),
                                                 code=open(Path(__file__).parent/'filterData.js').read()))

    download_button = Button(label="Download", button_type="success")
    download_button.js_on_event('button_click',
                             CustomJS(args=dict(source=source, filename="pangenomes_status.tsv"),
                                      code=open(Path(__file__).parent/'download.js').read()))


    layout = column(status_dt,
                    # row(*labels,  sizing_mode="scale_width", spacing=8),
                    row(Div(text='', width=160), row(*radio_buttons, spacing=29), Div(text='', width=20),
                        filter_button, download_button)
                    )

    curdoc().add_root(layout)

    html_content = file_html(layout, "cdn", title="Pangenome status Information")
    with open(output / 'status_info.html', "w") as file:
        file.write(html_content)


def export_info(info_dict: dict, output: Path):
    """ Export information to HTML file

    :param info_dict: Dictionary with information readable
    :param output: Path to output directory
    """
    for info, value in info_dict.items():
        if info == "status":
            export_status(value, output)

def launch(args: argparse.Namespace):
    """
    Command launcher

    :param args: All arguments provide by user
    """
    logging.debug("launch info command")
    pangenomes_to_path = check_tsv_sanity(args.pangenomes)
    info_dict = get_info(pangenomes_path=pangenomes_to_path, status=args.status, content=args.content,
                         parameters=args.parameters, metadata=args.metadata, disable_bar=args.disable_prog_bar)
    export_info(info_dict, args.output)
    logging.info("Done")


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
