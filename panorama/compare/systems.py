#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple
from shutil import rmtree
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock

# installed libraries
from tqdm import tqdm
import networkx as nx

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.utils import mkdir, init_lock
from panorama.utility.utility import check_models
from panorama.systems.system import System, ClusterSystems
from panorama.systems.write_systems import check_pangenome_write_systems
from panorama.compare.utils import parser_comparison, common_launch, cluster_on_frr, compute_frr


def check_compare_systems_args(args):
    need_info = {"need_annotations": True, "need_families": True, "need_families_info": True,  # "need_graph": True,
                 "need_metadata": True, "metatypes": ["families"], "need_systems": True,
                 "systems_sources": args.sources, "read_canonical": args.canonical}

    if len(args.sources) != len(args.models):
        raise argparse.ArgumentError(argument=None, message="Number of sources and models are different.")

    return need_info


def add_info_systems(pangenomes: Pangenomes, graph: nx.graph):
    for name, pangenome in pangenomes.items():
        for system in pangenome.systems:
            sys_info = {"pangenome": name,
                        "system_name": system.name,
                        "system_id": system.ID,
                        "families_models_count": system.number_of_model_gene_families,
                        "families_count": system.number_of_families}
            sys_hash = hash((name, system.name, system.ID))

            if graph.has_node(sys_hash):
                node_attributes = graph.nodes[sys_hash]
                node_attributes.update(sys_info)


def create_pangenome_system_graph(pangenome):
    graph = nx.Graph()
    sys2pangenome = {}
    syshash2sys = {}
    for system in pangenome.systems:
        sys_hash = hash((pangenome.name, system.name, system.ID))
        graph.add_node(sys_hash, system_id=system.ID, system_name=system.name, pangenome=pangenome.name)
        sys2pangenome[sys_hash] = pangenome.name
        syshash2sys[sys_hash] = system
    return graph, sys2pangenome, syshash2sys


def create_systems_graph(pangenomes: Pangenomes, threads: int = 1, lock: Lock = None, disable_bar: bool = False):
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit="Pangenome", disable=disable_bar) as pbar:
            futures = []
            for pangenome in pangenomes:
                logging.getLogger("PANORAMA").debug(f"Add spots for pangenome {pangenome.name}")
                future = executor.submit(create_pangenome_system_graph, pangenome)
                future.add_done_callback(lambda p: pbar.update())
                futures.append(future)
            systems_graph = nx.Graph()
            systems2pangenome = {}
            systemhash2system = {}
            for future in futures:
                res = future.result()
                systems_graph.add_nodes_from(res[0])
                systems2pangenome.update(res[1])
                systemhash2system.update(res[2])
    return systems_graph, systems2pangenome, systemhash2system


def compute_frr_edges(graph: nx.Graph, sys2pangenome: Dict[int, str], systemhash2system: Dict[int, System],
                      frr_cutoff: Tuple[float, float] = (0.8, 0.8), frr_models_cutoff: Tuple[float, float] = (0.8, 0.8),
                      disable_bar: bool = False):
    """
    Compute spots graph edges with frr score

    Args:
        graph: Spots graph with node only
        spots2pangenome: Dictionary of spot to pangenome from which they belong
        frr_cutoff: frr cutoff for frr score (default: 0.8)
        disable_bar: Flag to disable progress bar (default: False)
    """
    sys_pair = [(sys1_hash, sys2_hash) for sys1_hash, sys2_hash in combinations(graph.nodes, 2)
                if sys2pangenome[sys1_hash] != sys2pangenome[sys2_hash]]
    with tqdm(total=len(sys_pair), unit='system pair', desc="Compute frr", disable=disable_bar) as pbar:
        for sys1_hash, sys2_hash in sys_pair:
            if not graph.has_edge(sys1_hash, sys2_hash):
                sys1, sys2 = systemhash2system[sys1_hash], systemhash2system[sys2_hash]
                min_frr_models, max_frr_models, shared_models_gf = compute_frr(set(sys1.models_families),
                                                                               set(sys2.models_families))
                if min_frr_models > frr_models_cutoff[0] and max_frr_models > frr_models_cutoff[1]:
                    min_frr, max_frr, shared_gf = compute_frr(set(sys1.families), set(sys2.families))
                    if min_frr > frr_cutoff[0] and max_frr > frr_cutoff[1]:
                        graph.add_edge(sys1_hash, sys2_hash, min_frr_models=min_frr_models, max_frr_models=max_frr_models,
                                       shared_model_families=shared_models_gf, min_frr=min_frr, max_frr=max_frr,
                                       shared_families=shared_gf)
            pbar.update()


def write_conserved_systems(pangenomes: Pangenomes, output: Path, cs_graph: nx.Graph, graph_formats: List[str]):
    add_info_systems(pangenomes, cs_graph)

    if "gexf" in graph_formats:
        # writing graph in gexf format
        graph_file_name = output / "conserved_systems.gexf"
        logging.info(f"Writing graph in gexf format in {graph_file_name}.")
        nx.readwrite.gexf.write_gexf(cs_graph, graph_file_name)

    if "graphml" in graph_formats:
        graph_file_name = output / "conserved_systems.graphml"
        logging.info(f"Writing graph in graphml format in {graph_file_name}.")
        nx.readwrite.graphml.write_graphml(cs_graph, graph_file_name)

    outfile = output / "conserved_systems.tsv"
    logging.info(f"Writing rgp clusters in tsv format in {outfile}")


def compare_systems(pangenomes: Pangenomes, frr_metrics: str = "min_frr_models",
                    frr_cutoff: Tuple[float, float] = (0.8, 0.8), frr_models_cutoff: Tuple[float, float] = (0.8, 0.8),
                    threads: int = 1, lock: Lock = None, disable_bar: bool = False):
    systems_graph, systems2pangenome, systemhash2system = create_systems_graph(pangenomes, threads, lock, disable_bar)

    compute_frr_edges(systems_graph, systems2pangenome, systemhash2system, frr_cutoff, frr_models_cutoff)

    partitions = cluster_on_frr(systems_graph, frr_metrics)
    for cs_id, cluster_sys in enumerate(partitions, start=1):
        if len(cluster_sys) > 1:
            cs_sys = set()
            for sys_hash in cluster_sys:
                sys = systemhash2system[sys_hash]
                pangenome_name = systems2pangenome[sys_hash]
                node_attributes = systems_graph.nodes[sys_hash]
                del node_attributes[f"{frr_metrics}_cluster"]
                node_attributes.update({"cluster_systems_id": cs_id, "system_id": sys.ID, "system_name": sys.name,
                                        "pangenome": pangenome_name})
                cs_sys.add(sys)
            pangenomes.add_cluster_systems(ClusterSystems(cs_id, *cs_sys))
        else:
            systems_graph.remove_nodes_from(cluster_sys)

    return systems_graph


def launch(args):
    """
    Launch functions to align gene families from pangenomes

    Args:
        args: argument given in CLI
    """
    need_info = check_compare_systems_args(args)

    models_list = []
    for models in args.models:
        models_list.append(check_models(models, disable_bar=args.disable_prog_bar))

    need_info["models"] = models_list

    pangenomes, tmpdir, _, _ = common_launch(args, check_pangenome_write_systems, need_info, sources=args.sources)

    output = mkdir(args.output, force=args.force)

    cs_graph = compare_systems(pangenomes, frr_metrics=args.frr_metrics, frr_cutoff=args.frr_cutoff,
                               frr_models_cutoff=args.frr_models_cutoff, disable_bar=args.disable_prog_bar)

    write_conserved_systems(pangenomes, output, cs_graph, args.graph_formats)

    if not args.keep_tmp:
        rmtree(tmpdir, ignore_errors=True)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for cluster command

    Returns:
        argparse.ArgumentParser: parser arguments for cluster command
    """
    parser = sub_parser.add_parser("compare_systems",
                                   description='Comparison of systems among pangenomes')

    parser_comparison_systems(parser)
    return parser


def parser_comparison_systems(parser):
    """
    Add argument to parser for system comparison command

    Args:
        parser: parser for cluster argument
    """
    required, compare_opt, optional = parser_comparison(parser)
    required.add_argument('-m', '--models', required=True, type=Path, nargs="+",
                          help="Path to model list file. You can specify multiple models from different source. "
                               "For that separate the model list files by a space and "
                               "make sure you give them in the same order as the sources.")
    required.add_argument("-s", "--sources", required=True, type=str, nargs="+",
                          help="Name of the systems sources. You can specify multiple sources. "
                               "For that separate names by a space and "
                               "make sure you give them in the same order as the sources.")
    compare_opt.add_argument('--frr_metrics', required=False, type=str, default="min_frr_models",
                             choices=["min_frr_models", "max_frr_models", "min_frr", "max_frr"],
                             help="Metrics used to computed conserved systems cluster.")
    compare_opt.add_argument('--frr_models_cutoff', required=False, type=tuple, default=(0.2, 0.2), nargs=2,
                             help="The frr_models (Families Repertoire Relatedness) is used to assess the similarity "
                                  "between two systems based on their gene families in the models.\n"
                                  "\tThe 'min_frr_models': Computes the number of models gene families shared "
                                  "between the two systems and divides it by the smaller number of models gene families"
                                  " among the two systems.\n"
                                  "\tThe 'max_frr_models': Computes the number of models gene families shared "
                                  "between the two systems and divides it by the larger number of models gene families "
                                  "among the two systems.\n"
                             )
    optional.add_argument("--canonical", required=False, action="store_true",
                          help="Write the canonical version of systems too.")
