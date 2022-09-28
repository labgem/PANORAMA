#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, CancelledError, Executor, Future
from multiprocessing import Manager
import functools

# installed librairies
from ppanggolin.formats.readBinaries import check_pangenome_info
from py2neo import Graph
# local librairies
from panorama.pangenomes import Pangenome
from panorama.utils import check_tsv_sanity
from panorama.dbGraph.translatePan import create_dict
from panorama.panDB import PangenomeLoader


logger = logging.getLogger(__name__)


def init_db_lock(lock):
    global db_loading_lock
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


def load_pangenome(pangenome_name, pangenome_path):
    pangenome = Pangenome(name=pangenome_name)
    pangenome.add_file(pangenome_path)
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True,
                         need_modules=True, disable_bar=True)
    data = create_dict(pangenome)
    loader = PangenomeLoader(pangenome_name, data, db_loading_lock)
    loader.load(graph)


def load_pangenome_mp(pangenomes: dict, cpu: int = 1):
    manager = Manager()
    lock = manager.Lock()
    futures = []
    bar = tqdm(total=len(pangenomes), unit='pangenome')
    with ProcessPoolExecutor(max_workers=cpu, initializer=init_db_lock, initargs=(lock,)) as executor:
        for pangenome_name, pangenome_path in pangenomes.items():
            logging.getLogger().info(f"Add {pangenome_name} to schedule")
            future = executor.submit(load_pangenome, pangenome_name, pangenome_path)

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


# pandas.read_csv(config.METADATA_FILE)

# Simple singleprocessed loading
# def load_data():
#     dataloader = Dataloader(config.METADATA_FILE)
#     dataloader.parse()


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    global graph

    pan_to_path = check_tsv_sanity(args.pangenomes)
    graph = Graph(uri=args.uri, user=args.user, password=args.pwd)
    if args.clean:
        graph.delete_all()
    load_pangenome_mp(pan_to_path, args.cpu)



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
