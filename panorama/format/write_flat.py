#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
from pathlib import Path
# installed libraries
import numpy as np
import pandas as pd

# local libraries
from panorama.format.read_binaries import check_pangenome_info
from panorama.utils import check_tsv_sanity
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome


# def filter_df(df: pd.DataFrame):
#     current_gf = ''
#     index_access = range(4, df.shape[1] + 1, 2)
#     for row in df.itertuples():
#         if current_gf != row[2]:
#             current_gf = row[2]
#             access_dict = {i: set() for i in index_access}
#             for index in index_access:
#




def write_annotation_to_families(pangenome: Pangenome, output: Path):
    gf: GeneFamily
    # out_df = pd.DataFrame({"Pangenome": [], "Family": []})
    nb_source = len(pangenome.annotation_source)
    source_list = list(pangenome.annotation_source)
    column_name = np.array(f"Pangenome,Gene Family,{','.join([f'Annotation {source},Accession {source},Secondary names {source}' for source in source_list])}".split(','))
    array_list = []
    for gf in pangenome.gene_families:
        annot_array = np.empty((gf.max_annotation_by_source()[1], 2+nb_source*3), dtype=object)
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
                        annot_array[index_annot, index_source + 2] = annotation.secondary_names if annotation.secondary_names is not pd.NA else '-'
                        index_annot += 1
                index_source += 3
            array_list.append(annot_array)
    out_df = pd.DataFrame(np.concatenate(array_list), columns=column_name)
    out_df = out_df.sort_values(by=['Pangenome', 'Gene Family'] + list(column_name[range(3, len(column_name), 2)]))
    out_df.to_csv(f"{output}/families_annotations.tsv", sep="\t", header=True)


    # for source in pangenome.annotation_source:
    #     rows = []
    #     for gf in pangenome.gene_families:
    #         row_base = [pangenome.name, gf.name]
    #         if gf.get_source(source) is not None:
    #             for annotation in gf.get_source(source):
    #                 rows.append(row_base + [annotation.name, annotation.accession])
    #     out_df = out_df.merge(pd.DataFrame(rows, columns=["Pangenome", "Family", f"Annotation_{source}",
    #                                                       f"Accession_{source}"]), how='outer',
    #                           on=["Pangenome", "Family"])
    # filter_df(out_df)
    # out_df.sort_values(by=['Pangenome', 'Family']).to_csv(f"{output}/families_annotations.tsv", sep="\t", header=True)


def write_flat_files(pangenome, output: Path, annotation: bool = False, disable_bar: bool = False):
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

    check_pangenome_info(pangenome, need_annotations=need_annotations, need_families=need_families,
                         need_graph=need_graph, need_partitions=need_partitions, need_rgp=need_regions,
                         need_spots=need_spots, need_gene_sequences=need_gene_sequences, need_modules=need_modules,
                         need_annotation_fam=need_annotations_fam, disable_bar=disable_bar)

    if annotation:
        write_annotation_to_families(pangenome, output)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    # check_parameter(args)
    # systems_to_path = args.systems.absolute()
    # systems = read_systems(systems_to_path)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    for pangenome_name, pangenome_info in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
        pangenome.add_file(pangenome_info["path"])
        write_flat_files(pangenome, output=args.output, annotation=args.annotation, disable_bar=args.disable_prog_bar)


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
