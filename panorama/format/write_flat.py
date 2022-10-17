#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import csv
import json
import logging
from pathlib import Path
import tempfile
import re

import tables
from ppanggolin.genome import Organism
from tqdm import tqdm
# installed libraries
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.utils import check_tsv_sanity
from panorama.format.read_binaries import check_pangenome_info


def write_annotation_to_families(pangenome: Pangenome, output: Path):
    rows = []
    family: GeneFamily
    for family in pangenome.gene_families:
        row_base = [pangenome.name, family.name]
        for source, annotations in family.annotation.items():
            row_source = row_base + [source]
            for annotation in annotations:
                rows.append(row_source + [annotation])

    out_df = pd.DataFrame(rows, columns=["Pangenome", "Family", "Source", "Annotation"])
    out_df.to_csv(f"{output}/families_annotations.tsv", sep="\t", header=True)


def write_flat_files(pangenome, output: Path, annotation: bool = False, disable_bar: bool = False):
    needAnnotations = False
    needFamilies = False
    needGraph = False
    needPartitions = False
    needSpots = False
    needRegions = False
    needModules = False
    needGeneSequences=False
    needAnnotationsFam = False

    if annotation:
        needFamilies = True
        needAnnotationsFam = True

    check_pangenome_info(pangenome, need_annotations=needAnnotations, need_families=needFamilies, need_graph=needGraph,
                         need_partitions=needPartitions, need_rgp=needRegions, need_spots=needSpots,
                         need_modules=needModules, need_gene_sequences=needGeneSequences,
                         need_anntation_fam=needAnnotationsFam, disable_bar=disable_bar)

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
    for pangenome_name, pangenome_file in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name)
        pangenome.add_file(pangenome_file)
        write_flat_files(pangenome, output=args.output, annotation=args.annotation, disable_bar=args.disable_prog_bar)
    #     annot2fam = annot_pangenome(pangenome=pangenome, hmm=args.hmm, tsv=args.tsv, source=args.source,
    #                                 meta=args.meta, e_value=args.e_value, prediction_size=args.prediction_size,
    #                                 tmpdir=args.tmpdir, threads=args.threads, force=args.force,
    #                                 disable_bar=args.disable_prog_bar)
    #     search_system(systems, annot2fam, args.disable_prog_bar)
    #     logging.getLogger().info("Annotation Done")
    #     logging.getLogger().info(f"Write Annotation in pangenome {pangenome_name}")
    #     h5f = tables.open_file(pangenome_file, "a")
    #     write_gene_fam_annot(pangenome, h5f, force=args.force, disable_bar=args.disable_prog_bar)
    #     h5f.close()



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
