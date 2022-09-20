#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from tqdm import tqdm
from pathlib import Path

# installed librairies
from ppanggolin.formats.readBinaries import check_pangenome_info
# local librairies
from panorama.pangenomes import Pangenome
from panorama.utils import check_tsv_sanity
from panorama.dbGraph.translatePan import create_dict
from panorama.panDB import Neo4jDB
from panorama.dbGraph.exportPan import load_data


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pan_to_path = check_tsv_sanity(args.pangenomes)
    neo4j_db = Neo4jDB(uri=args.uri,
                       user=args.user,
                       pwd=args.pwd)
    if args.clean:
        neo4j_db.clean()
    for pangenome_name, pangenome_path in tqdm(pan_to_path.items(), unit='pangenome'):
        pangenome = Pangenome(name=pangenome_name)
        pangenome.add_file(pangenome_path)
        check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True)
        logging.getLogger().info(f"Begin {pangenome_name} translation for Neo4J")
        translate_pangenome = create_dict(pangenome)
        logging.getLogger().info(f"Begin {pangenome_name} load into Neo4J DB")
        load_data(db=neo4j_db, pangenome_data=translate_pangenome, threads=args.threads)
        neo4j_db.close()
        logging.getLogger().info("Load finished and connection to Neo4J DB closed")
    neo4j_db.close()


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("graph-db", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_pangraph_db(parser)
    return parser


def parser_pangraph_db(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    neo4j_arg = parser.add_argument_group(title="Required arguments for Neo4J DB")
    neo4j_arg.add_argument('--uri', required=True, type=str, nargs='?',
                           help='URI to connect to the database')
    neo4j_arg.add_argument('--user', required=True, type=str, nargs='?',
                           help='user name to connect to the database')
    neo4j_arg.add_argument('--pwd', required=True, type=str, nargs='?',
                           help='password to connect to the database')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--clean", required=False, action='store_true', default=False,
                          help="Clean the database before to add pangenomes. "
                               "WARNING, all the pangenomes will be removed.")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of av available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_pangraph_db(main_parser)
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
