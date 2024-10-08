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

# installed libraries
from tqdm import tqdm
import networkx as nx

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.utils import mkdir
from panorama.utility.utility import check_models
from panorama.systems.system import System
from panorama.systems.write_systems import check_pangenome_write_systems
from panorama.compare.utils import parser_comparison, common_launch, cluster_on_frr


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

            node_attributes = graph.nodes[sys_hash]
            node_attributes.update(sys_info)


def compute_frr(queries: Set[GeneFamily], targets: Set[GeneFamily]) -> Tuple[float, float, int]:
    akins = {query_gf.akin.ID for query_gf in queries for target_gf in targets if query_gf.akin == target_gf.akin}
    min_frr = len(akins) / min(len(queries), len(targets))
    max_frr = len(akins) / max(len(queries), len(targets))
    if min_frr > 1 or max_frr > 1:
        print("hello")

    return min_frr, max_frr, len(akins)


def compute_edge_metrics(query: System, target: System):
    min_frr_models, max_frr_models, shared_models_gf = compute_frr(set(query.models_families),
                                                                   set(target.models_families))
    min_frr, max_frr, shared_gf = compute_frr(set(query.families), set(target.families))

    return {"min_frr_models": min_frr_models,
            "max_frr_models": max_frr_models,
            "shared_model_families": shared_models_gf,
            "min_frr": min_frr,
            "max_frr": max_frr,
            "shared_families": shared_gf}


def compare_pangenomes_pair_system(query: Pangenome, target: Pangenome, system2id: Dict[str, Set[int]],
                                   graph: nx.Graph[int], min_frr_cutoff: float = 0.5, max_frr_cutoff: float = 0.8,
                                   min_frr_models_cutoff: float = 0.5, max_frr_models_cutoff: float = 0.8):
    for query_sys in query.systems:
        query_hash = hash((query.name, query_sys.name, query_sys.ID))
        if query_sys.name in system2id:
            for sys_id in system2id[query_sys.name]:
                target_sys = target.get_system(sys_id)
                target_hash = hash((target.name, target_sys.name, target_sys.ID))
                if not graph.has_edge(query_hash, target_hash):
                    edges_metrics = compute_edge_metrics(query_sys, target_sys)
                    if (edges_metrics["min_frr_models"] > min_frr_models_cutoff and edges_metrics[
                        "max_frr_models"] > max_frr_models_cutoff
                            and edges_metrics["min_frr"] > min_frr_cutoff and edges_metrics[
                                "max_frr"] > max_frr_cutoff):
                        graph.add_edge(query_hash, target_hash, **edges_metrics)


def compare_systems(pangenomes: Pangenomes, output: Path, tmpdir: Path, graph_formats: List[str],
                    frr_metrics: str = "min_frr_models", frr_cutoff: Tuple[float, float] = (0.8, 0.8),
                    frr_models_cutoff: Tuple[float, float] = (0.8, 0.8), cpus: int = 1, disable_bar: bool = False):
    def search_sys_id(target_pan):
        system2id[target_pan.name] = defaultdict(set)
        for system in target_pan.systems:
            system2id[target_pan.name][system.name].add(system.ID)

    system2id = {}
    conserved_systems_graph = nx.Graph()
    for pangenome in pangenomes:
        conserved_systems_graph.add_nodes_from({hash((pangenome.name, system.name, system.ID))
                                                for system in pangenome.systems})

    for query, target in tqdm(combinations(pangenomes, 2), total=(len(pangenomes)*(len(pangenomes)-1)/2),
                              unit="pair", disable=disable_bar):
        if target.name not in system2id:
            search_sys_id(target)
        compare_pangenomes_pair_system(query, target, system2id[target.name], conserved_systems_graph,
                                       min_frr_cutoff=frr_cutoff[0], max_frr_cutoff=frr_cutoff[1],
                                       min_frr_models_cutoff=frr_models_cutoff[0],
                                       max_frr_models_cutoff=frr_models_cutoff[1])

    cluster_on_frr(conserved_systems_graph, frr_metrics)

    add_info_systems(pangenomes, conserved_systems_graph)

    if "gexf" in graph_formats:
        # writing graph in gexf format
        graph_file_name = output / "conserved_systems.gexf"
        logging.info(f"Writing graph in gexf format in {graph_file_name}.")
        nx.readwrite.gexf.write_gexf(conserved_systems_graph, graph_file_name)

    if "graphml" in graph_formats:
        graph_file_name = output / "conserved_systems.graphml"
        logging.info(f"Writing graph in graphml format in {graph_file_name}.")
        nx.readwrite.graphml.write_graphml(conserved_systems_graph, graph_file_name)

    outfile = output / "conserved_systems.tsv"
    logging.info(f"Writing rgp clusters in tsv format in {outfile}")

    # write_conserved_systems(outfile, conserved_systems_graph, rgps_in_graph, grr_metric, rgp_to_spot)


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

    compare_systems(pangenomes, output, tmpdir=args.tmpdir, graph_formats=args.graph_formats,
                    frr_metrics=args.frr_metrics, frr_cutoff=args.frr_cutoff, frr_models_cutoff=args.frr_models_cutoff,
                    cpus=args.cpus, disable_bar=args.disable_prog_bar)

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

    parser_comparison_context(parser)
    return parser


def parser_comparison_context(parser):
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
    compare_opt.add_argument('--frr_cutoff', required=False, type=tuple, default=(0.2, 0.2), nargs=2,
                             help="The frr (Families Repertoire Relatedness) is used to assess the similarity between two "
                                  "systems based on their gene families.\n"
                                  "\tThe 'min_frr': Computes the number of gene families shared between the two elements "
                                  "and divides it by the smaller number of gene families among the two elements.\n"
                                  "\tThe 'max_frr': Computes the number of gene families shared between the two elements "
                                  "and divides it by the larger number of gene families among the two elements."
                             )
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
    optional.add_argument('--graph_formats', required=False, type=str, choices=['gexf', "graphml"], nargs="+",
                          default=['gexf', 'graphml'], help="Format of the output graph.")
