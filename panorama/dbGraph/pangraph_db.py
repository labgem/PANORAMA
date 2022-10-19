#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging

import py2neo
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, CancelledError, Executor, Future
from multiprocessing import Manager, Lock
import functools

# installed librairies
from py2neo import Graph
from graphio import RelationshipSet
import pandas as pd

# local librairies
from panorama.pangenomes import Pangenome
from panorama.utils import check_tsv_sanity
from panorama.dbGraph.translatePan import create_dict
from panorama.panDB import PangenomeLoader, SimilarityLoader
from panorama.format.read_binaries import check_pangenome_info

logger = logging.getLogger(__name__)

db_loading_lock: Lock = None


def invert_edges_query(edge_label: str):
    return f"""
    MATCH (f)-[r:{edge_label}]->(s)
    CALL apoc.refactor.invert(r)
    yield input, output
    RETURN input, output"""

def invert_edges(edge_label: str):
    query = invert_edges_query(edge_label)
    try:
        tx = graph.begin()
        tx.run(query)
        tx.commit()
    except Exception:
        raise Exception("Pb to invert edges")


def load_similarities(is_similar_to, batch_size):
    is_similar_to.merge(graph=graph, batch_size=batch_size)


def load_similarities_mp(tsv: Path, cpu: int = 1, batch_size: int = 1000):
    global fam_index, fam_known, fam_list
    df = pd.read_csv(filepath_or_buffer=tsv, sep="\t", header=None,
                     names=["Family_1", "Family_2", "identity", "covery"])

    # manager = Manager()
    # lock = manager.Lock()
    # init_db_lock(lock)
    is_similar_list = []
    for start_label in ['Persistent', 'Shell', 'Cloud']:
        for end_label in ['Persistent', 'Shell', 'Cloud']:
            is_similar_to = RelationshipSet('IS_SIMILAR_TO', [start_label], [end_label], ['name'], ['name'])
            for row in df.iterrows():
                is_similar_to.add_relationship(start_node_properties={"name": row[1]['Family_1']},
                                               end_node_properties={"name": row[1]['Family_2']},
                                               properties={"identity": row[1]['identity'],
                                                           "covery": row[1]['covery']})
            is_similar_list.append(is_similar_to)
    with ProcessPoolExecutor(max_workers=cpu) as executor:
        list(tqdm(executor.map(load_similarities, is_similar_list, [batch_size]*len(is_similar_list)),
             unit="similarities_batch", total=len(is_similar_list)))


def init_db_lock(lock):
    global db_loading_lock
    if db_loading_lock is None:
        db_loading_lock = lock


def load_pangenome_done(pangenome_name, executor: Executor, bar, future: Future):
    try:
        future.result()
    except CancelledError:
        # canceled by god or anyone
        logging.getLogger().info(f"Load {pangenome_name} cancelled")
        return
    except Exception as error:
        # if config.CANCEL_WHOLE_IMPORT_IF_A_WORKER_FAILS:
        logging.getLogger().warning(f"Load {pangenome_name} failed. Cancel all tasks and stop workers...")
        executor.shutdown(cancel_futures=True)

        logging.getLogger().info(f"{pangenome_name} failed")
        logging.getLogger().exception(f"Load {pangenome_name} raised {error}")
        raise error
    bar.update()
    logging.getLogger().info(f"Load {pangenome_name} finished")


def load_pangenome(pangenome_name, pangenome_info):
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True,
                         need_rgp=True, need_spots=True, need_modules=True, need_anntation_fam=True,
                         disable_bar=False)
    data = create_dict(pangenome)
    loader = PangenomeLoader(pangenome_name, data, db_loading_lock)
    loader.load(graph)


def load_pangenome_mp(pangenomes: dict, cpu: int = 1):
    manager = Manager()
    lock = manager.Lock()
    futures = []
    bar = tqdm(total=len(pangenomes), unit='pangenome')
    with ProcessPoolExecutor(max_workers=cpu, initializer=init_db_lock, initargs=(lock,)) as executor:
        for pangenome_name, pangenome_info in pangenomes.items():
            logging.getLogger().info(f"Add {pangenome_name} to schedule")
            future = executor.submit(load_pangenome, pangenome_name, pangenome_info)

            future.add_done_callback(
                functools.partial(
                    load_pangenome_done,
                    pangenome_name,
                    executor,
                    bar
                )
            )
            futures.append(future)
    bar.close()


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    global graph

    pangenomes = check_tsv_sanity(args.pangenomes)
    graph = Graph(uri=args.uri, user=args.user, password=args.pwd)
    if args.clean:
        graph.delete_all()
        # tx = graph.begin()
        # tx.run("create constraint on f:Family assert f.hash_id is unique;")
        # graph.commit(tx)
    load_pangenome_mp(pangenomes, args.cpu)
    logging.getLogger().info("All pangenomes loaded...")
    logging.getLogger().info("Bengin load of similarities...")
    load_similarities_mp(args.similarities)
    logging.getLogger().info("Invert edges...")
    labels2invert = ["IS_IN_PANGENOME", "IS_IN_MODULE", "IS_IN_PERSISTENT", "IS_IN_SHELL",
                     "IS_IN_CLOUD", "IS_IN_CONTIG", "IS_IN_GENOME", "IS_IN_SPOT", "IS_IN_RGP"]
    for edge_label in labels2invert:
        logging.getLogger().debug(f"Invert: {edge_label}")
        invert_edges(edge_label)


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
    optional.add_argument("--cpu", required=False, nargs='?', type=int, default=1,
                          help="Number of available cpu")
    optional.add_argument("--similarities", required=False, nargs='?', type=Path, default=None,
                          help="PANORAMA can perform an alignment of all gene families. In this case,"
                               "3 identiy levels will be used to have complete information."
                               "If you want you can give your own alignment file.")


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
