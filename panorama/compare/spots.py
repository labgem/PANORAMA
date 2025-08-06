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
import numpy as np
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.geneFamily import GeneFamily
from panorama.region import Spot, ConservedSpots
from panorama.utils import mkdir, init_lock
from panorama.compare.utils import parser_comparison, common_launch, cluster_on_frr, compute_frr
from panorama.utility.utility import check_models


def check_compare_spots_args(args):
    need_info = {"need_annotations": True, "need_families": True, "need_families_info": True,  # "need_graph": True,
                 "need_rgp": True, "need_spots": True}
    if args.systems:
        if args.models is not None:
            if args.sources is None:
                raise argparse.ArgumentError(argument=None, message="You required to add systems to comparison but "
                                                                    "did not provide sources")
            else:
                models_list = []
                for models in args.models:
                    models_list.append(check_models(models, disable_bar=args.disable_prog_bar))
                need_info.update({"need_metadata": True, "metatypes": ["families"], "need_systems": True,
                                  "systems_sources": args.sources, "models": models_list,
                                  "read_canonical": args.canonical})
        else:
            if args.sources is not None:
                raise argparse.ArgumentError(argument=None, message="You required to add systems to comparison but "
                                                                    "did not provide models")
            else:
                raise argparse.ArgumentError(argument=None, message="You required to add systems to comparison but "
                                                                    "did not provide models and sources")
    return need_info


def check_pangenome_cs(pangenome: Pangenome, sources: List[str] = None):
    if pangenome.status["predictedRGP"] not in ["inFile", "Computed", "Loaded"]:
        raise ValueError(f"RGPs did not have been predicted in your pangenome {pangenome.name}."
                         "Please look at 'ppanggolin rgp' command")
    if pangenome.status["spots"] not in ["inFile", "Computed", "Loaded"]:
        raise ValueError(f"Spots did not have been predicted in your pangenome {pangenome.name}."
                         "Please look at 'ppanggolin spot' command")
    if sources is not None:
        if pangenome.status["systems"] != "inFile":
            raise AttributeError("Systems have not been detected."
                                 "Use 'panorama detect' subcommand to detect systems in pangenomes.")
        else:
            for systems_source in sources:
                if systems_source not in pangenome.status["systems_sources"]:
                    logging.getLogger("PANORAMA").error(f"Systems in pangenome {pangenome.name} are: "
                                                        f"{pangenome.status['systems_sources']}")
                    raise KeyError(
                        f"There is no systems in pangenome {pangenome.name}, for the source: {systems_source}."
                        f"Look at 'panorama detect' subcommand to detect systems in {pangenome.name}.")


def create_pangenome_spots_graph(pangenome: Pangenome, dup_margin: float = 0.05
                                 ) -> Tuple[nx.Graph, Dict[int, Set[GeneFamily]], Dict[int, str], Dict[int, Spot]]:
    """
    Create a graph of spots belonging to one pangenome

    Args:
        pangenome: Pangenome associated with spots
        dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated (default = 0.05)

    Returns:
        The spot graph without edges and a dictionary of spots links to their bordering gene families
    """

    def get_borders_families() -> Set[GeneFamily]:
        """
        Get all bordering gene families from a spot
        Returns:
            Set of bordering gene families
        """
        borders_families = set()
        borders = spot.borders(pangenome.parameters["spot"]["set_size"], multigenic)
        for _, border in borders:
            borders_families |= set(border[0])
            borders_families |= set(border[1])
        return borders_families

    graph = nx.Graph()
    spots2borders = {}
    spots2pangenome = {}
    spothash2spot = {}
    multigenic = pangenome.get_multigenics(dup_margin=dup_margin)
    for spot in pangenome.spots:
        spot_hash = hash((spot.ID, pangenome.name))
        graph.add_node(spot_hash, spot_id=spot.ID, pangenome=pangenome.name)
        spots2borders[spot_hash] = get_borders_families()
        spots2pangenome[spot_hash] = pangenome.name
        spothash2spot[spot_hash] = spot
    return graph, spots2borders, spots2pangenome, spothash2spot


def create_spots_graph(pangenomes: Pangenomes, dup_margin: float = 0.05, threads: int = 1, lock: Lock = None,
                       disable_bar: bool = False
                       ) -> Tuple[nx.Graph, Dict[int, Set[GeneFamily]], Dict[int, str], Dict[int, Spot]]:
    """
    Create a graph with the spots from all pangenomes as nodes. There are no edges computed.

    Args:
        pangenomes: Pangenomes object containing all pangenomes
        dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated (default = 0.05)
        threads: Available threads (default = 1)
        lock: Lock object (default = None)
        disable_bar: Flag to disable progress bar (default = False)

    Returns:
        The spot graph without edges and a dictionary of spots links to their bordering gene families
    """
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit="Pangenome", disable=disable_bar) as pbar:
            futures = []
            for pangenome in pangenomes:
                logging.getLogger("PANORAMA").debug(f"Add spots for pangenome {pangenome.name}")
                future = executor.submit(create_pangenome_spots_graph, pangenome, dup_margin)
                future.add_done_callback(lambda p: pbar.update())
                futures.append(future)
            spots_graph = nx.Graph()
            spots2borders = {}
            spots2pangenome = {}
            spothash2spot = {}
            for future in futures:
                res = future.result()
                spots_graph.add_nodes_from(res[0])
                spots2borders.update(res[1])
                spots2pangenome.update(res[2])
                spothash2spot.update(res[3])
    return spots_graph, spots2borders, spots2pangenome, spothash2spot


def compute_frr_edges(graph: nx.Graph, spots2borders: Dict[int, Set[GeneFamily]], spots2pangenome: Dict[int, str],
                      min_frr_cutoff: float = 0.5, max_frr_cutoff: float = 0.8, disable_bar: bool = False):
    """
    Compute spots graph edges with frr score

    Args:
        graph: Spots graph with node only
        spots2borders: Dictionary of spot link to their bordering families
        spots2pangenome: Dictionary of spot to pangenome from which they belong
        frr_cutoff: frr cutoff for frr score (default: 0.8)
        disable_bar: Flag to disable progress bar (default: False)
    """
    spots_pair = [(spot1, spot2) for spot1, spot2 in combinations(graph.nodes, 2)
                  if spots2pangenome[spot1] != spots2pangenome[spot2]]
    with tqdm(total=len(spots_pair), unit='spots pair', desc="Compute frr", disable=disable_bar) as pbar:
        for spot1, spot2 in spots_pair:
            if not graph.has_edge(spot1, spot2):
                border1_fams = {fam for fam in spots2borders[spot1]}
                border2_fams = {fam for fam in spots2borders[spot2]}

                min_frr, max_frr, shared_gf = compute_frr(border1_fams, border2_fams)
                if min_frr > min_frr_cutoff and max_frr > max_frr_cutoff:
                    graph.add_edge(spot1, spot2, min_frr=min_frr, max_frr=max_frr, shared_families=shared_gf)
            pbar.update()


def add_systems_info(pangenomes, cs_graph):
    sys_names = set()
    node2nbsys = defaultdict(lambda: 0)
    for pangenome in pangenomes:
        for system in pangenome.systems:
            spots = set()
            sys_names.add(system.name)
            for gf in system.model_families:
                spots |= set(gf.spots)
            for spot in spots:
                spot_hash = hash((spot.ID, pangenome.name))
                if cs_graph.has_node(spot_hash):
                    nodes_attributes = cs_graph.nodes[spot_hash]
                    nodes_attributes[system.name] = True
                    node2nbsys[spot_hash] += 1

    for node in cs_graph.nodes:
        nodes_attributes = cs_graph.nodes[node]
        nodes_attributes.update({sys: False for sys in sys_names if sys not in nodes_attributes})
        nodes_attributes["#systems"] = node2nbsys[node]


def create_pangenome_system_graph(pangenome):
    graph = nx.Graph()
    sys2pangenome = {}
    syshash2sys = {}
    systemhash2conserved_spots = defaultdict(set)
    syshash2orgs = defaultdict(set)
    for system in pangenome.systems:
        for spot in system.spots:
            sys_hash = hash((pangenome.name, system.name, spot.ID))
            if sys_hash not in syshash2sys:
                syshash2orgs[sys_hash] |= set(system.organisms)
                graph.add_node(sys_hash, system_name=system.name, pangenome=pangenome.name, spot=spot.ID)
                sys2pangenome[sys_hash] = pangenome.name
                syshash2sys[sys_hash] = system
            if spot.conserved_id:
                systemhash2conserved_spots[sys_hash].add(spot.conserved_id)
    nx.set_node_attributes(graph, {sys_hash: (len(n_orgs) / pangenome.number_of_organisms) * 100
                                   for sys_hash, n_orgs in syshash2orgs.items()}, "percent_org")
    return graph, sys2pangenome, syshash2sys, systemhash2conserved_spots


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
            systemhash2conserved_spots = defaultdict(set)
            for future in futures:
                res = future.result()
                systems_graph.add_nodes_from(res[0].nodes(data=True))
                systems2pangenome.update(res[1])
                systemhash2system.update(res[2])
                systemhash2conserved_spots.update(res[3])
    return systems_graph, systems2pangenome, systemhash2system, systemhash2conserved_spots


def graph_systems_link_with_conserved_spots(pangenomes: Pangenomes, output: Path,
                                            graph_formats: List[str] = None, threads: int = 1,
                                            lock: Lock = None,
                                            disable_bar: bool = False):
    systems_graph, system2pangenome, systemhash2system, systemhash2conserved_spots = create_systems_graph(pangenomes,
                                                                                                          threads, lock,
                                                                                                          disable_bar)

    for sys_hash1, sys_hash2 in tqdm(combinations(systems_graph.nodes, 2),
                                     total=len(systems_graph.nodes)*(len(systems_graph.nodes)+1)/2):
        common_cs = systemhash2conserved_spots[sys_hash1].intersection(systemhash2conserved_spots[sys_hash2])
        if len(common_cs) > 0:
            systems_graph.add_edge(sys_hash1, sys_hash2, cluster_spots=",".join(map(str, sorted(common_cs))),
                                   weight=len(common_cs))
            node_attr1, node_attr2 = systems_graph.nodes[sys_hash1], systems_graph.nodes[sys_hash2]
            pangenome1, pangenome2 = system2pangenome[sys_hash1], system2pangenome[sys_hash2]
            sys1, sys2 = systemhash2system[sys_hash1], systemhash2system[sys_hash2]
            node_attr1.update({"system_name": sys1.name, "pangenome": pangenome1})
            node_attr2.update({"system_name": sys2.name, "pangenome": pangenome2})
            if "spots" in node_attr1:
                node_attr1["spots"] |= {spot.ID for cs in common_cs
                                        for spot in pangenomes.get_conserved_spots(cs).spots
                                        if spot.pangenome.name == pangenome1}
            else:
                node_attr1["spots"] = {spot.ID for cs in common_cs for spot in
                                       pangenomes.get_conserved_spots(cs).spots
                                       if spot.pangenome.name == pangenome1}
            if "spots" in node_attr2:
                node_attr2["spots"] |= {spot.ID for cs in common_cs
                                        for spot in pangenomes.get_conserved_spots(cs).spots
                                        if spot.pangenome.name == pangenome2}
            else:
                node_attr2["spots"] = {spot.ID for cs in common_cs for spot in
                                       pangenomes.get_conserved_spots(cs).spots
                                       if spot.pangenome.name == pangenome2}

    systems_graph.remove_nodes_from(list(nx.isolates(systems_graph)))

    for node in systems_graph.nodes:
        node_attr = systems_graph.nodes[node]
        if "spots" in node_attr:
            node_attr["spots"] = ",".join(map(str, sorted(node_attr["spots"])))

    louvain_graph = systems_graph.copy()
    partitions = nx.algorithms.community.louvain_communities(louvain_graph, weight="weight")

    cluster_id = 1
    for _, cluster_systems in enumerate(partitions, start=1):
        if len(cluster_systems) > 1:
            pan_set = set()
            for sys_hash in cluster_systems:
                node_attr = louvain_graph.nodes[sys_hash]
                pangenome = system2pangenome[sys_hash]
                sys = systemhash2system[sys_hash]
                node_attr.update({"system_ID": sys.ID, "system_name": sys.name,
                                  "pangenome": pangenome, "cluster_id": cluster_id})
                pan_set.add(pangenome)
            if len(pan_set) == 1:
                louvain_graph.remove_nodes_from(cluster_systems)
            else:
                cluster_id += 1
        else:
            louvain_graph.remove_nodes_from(cluster_systems)

    if graph_formats is not None:
        # add_info_spots(pangenomes, cs_graph)
        if "gexf" in graph_formats:
            # writing graph in gexf format
            graph_file_name = output / "systems_link_with_conserved_spots_louvain.gexf"
            logging.info(f"Writing graph in gexf format in {graph_file_name}.")
            nx.readwrite.gexf.write_gexf(louvain_graph, graph_file_name)

        if "graphml" in graph_formats:
            graph_file_name = output / "systems_link_with_conserved_spots_louvain.graphml"
            logging.info(f"Writing graph in graphml format in {graph_file_name}.")
            nx.readwrite.graphml.write_graphml(louvain_graph, graph_file_name)

    # Build the Minimum Spanning Tree (MST)
    mst = nx.minimum_spanning_tree(systems_graph, weight="weight")

    # Extract edge weights from the MST
    weights = np.array([d['weight'] for _, _, d in mst.edges(data=True)])

    # Analyze weights to determine an automatic threshold
    sorted_weights = np.sort(weights)  # Sort the weights
    diffs = np.diff(sorted_weights)  # Compute differences between consecutive weights
    max_jump_index = np.argmax(diffs)  # Identify the largest jump in weights
    threshold = sorted_weights[max_jump_index]  # Define the threshold as the weight before the largest jump

    # Remove edges heavier than the threshold
    to_remove = [(u, v) for u, v, d in mst.edges(data=True) if d["weight"] > threshold]
    mst.remove_edges_from(to_remove)

    for cluster_id, cluster_systems in enumerate(nx.connected_components(mst), start=1):
        if len(cluster_systems) > 1:
            for sys_hash in cluster_systems:
                node_attr = mst.nodes[sys_hash]
                pangenome = system2pangenome[sys_hash]
                sys = systemhash2system[sys_hash]
                node_attr.update({"system_name": sys.name, "pangenome": pangenome, "cluster_id": cluster_id})
        else:
            mst.remove_nodes_from(cluster_systems)

    if graph_formats is not None:
        # add_info_spots(pangenomes, cs_graph)
        if "gexf" in graph_formats:
            # writing graph in gexf format
            graph_file_name = output / "systems_link_with_conserved_spots_mst.gexf"
            logging.info(f"Writing graph in gexf format in {graph_file_name}.")
            nx.readwrite.gexf.write_gexf(mst, graph_file_name)

        if "graphml" in graph_formats:
            graph_file_name = output / "systems_link_with_conserved_spots_mst.graphml"
            logging.info(f"Writing graph in graphml format in {graph_file_name}.")
            nx.readwrite.graphml.write_graphml(mst, graph_file_name)


def write_conserved_spots(pangenomes, output: Path, graph_formats: List[str] = None, cs_graph: nx.Graph = None,
                          force: bool = False, disable_bar: bool = False):
    """
    Write conserved spots into files

    Args:
        pangenomes: Pangenomes associated to conserved spots
        output: Path to the output directory
        force: Flag to overwrite existing files (default: False)
        disable_bar: Flag to disable progress bar (default: False)
    """
    all_cs = []
    cs_dir = mkdir(output / "conserved_spots", force=force, erase=force)
    for conserved_spot in tqdm(pangenomes.conserved_spots, total=pangenomes.number_of_conserved_spots,
                               disable=disable_bar):
        by_cs = []
        for spot in conserved_spot.spots:
            all_cs.append([conserved_spot.ID, spot.ID, spot.pangenome.name, len(spot), spot.number_of_families])
            for rgp in spot.regions:
                by_cs.append([spot.ID, spot.pangenome.name, rgp.name, ",".join([fam.name for fam in rgp.families])])
            by_cs_df = pd.DataFrame(by_cs, columns=['Spot', 'Pangenome', "RGP", "Families"])
            by_cs_df = by_cs_df.sort_values(by=['Spot', 'Pangenome', 'RGP'])
            by_cs_df.to_csv(cs_dir / f"conserved_spots_{conserved_spot.ID}.tsv", sep="\t", header=True, index=False)
    conserved_df = pd.DataFrame(all_cs, columns=['Conserved ID', 'Spot ID', 'Pangenome', "#RGP", "#Families"])
    conserved_df = conserved_df.sort_values(by=['Conserved ID', 'Spot ID', 'Pangenome', '#RGP', "#Families"])
    conserved_df.to_csv(output / "all_conserved_spots.tsv", sep="\t", header=True, index=False)

    add_systems_info(pangenomes, cs_graph)

    if graph_formats is not None:
        assert cs_graph is not None, AssertionError("No graph to write spots cluster")
        # add_info_spots(pangenomes, cs_graph)
        if "gexf" in graph_formats:
            # writing graph in gexf format
            graph_file_name = output / "conserved_spots.gexf"
            logging.info(f"Writing graph in gexf format in {graph_file_name}.")
            nx.readwrite.gexf.write_gexf(cs_graph, graph_file_name)

        if "graphml" in graph_formats:
            graph_file_name = output / "conserved_spots.graphml"
            logging.info(f"Writing graph in graphml format in {graph_file_name}.")
            nx.readwrite.graphml.write_graphml(cs_graph, graph_file_name)


def compare_spots(pangenomes: Pangenomes, dup_margin: float = 0.05, frr_metrics: str = "min_frr",
                  frr_cutoff: Tuple[float, float] = (0.8, 0.8), threads: int = 1,
                  lock: Lock = None, disable_bar: bool = False):
    """
    Main function to identify conserved spots between pangenomes and add them into pangenomes.

    Args:
        pangenomes: Pangenomes object containing pangenome
        dup_margin: minimum ratio of organisms in which family must have multiple genes to be considered duplicated (default = 0.05)
        threads: Available threads (default = 1)
        lock: Lock object (default = None)
        disable_bar: Flag to disable progress bar (default = False)
    """
    spots_graph, spots2borders, spots2pangenome, spothash2spot = create_spots_graph(pangenomes, dup_margin, threads,
                                                                                    lock, disable_bar)

    compute_frr_edges(spots_graph, spots2borders, spots2pangenome, min_frr_cutoff=frr_cutoff[0],
                      max_frr_cutoff=frr_cutoff[1], disable_bar=disable_bar)

    partitions = cluster_on_frr(spots_graph, frr_metrics)
    for cs_id, cluster_spots in enumerate(partitions, start=1):
        if len(cluster_spots) > 1:
            cs_spots = set()
            for spot_hash in cluster_spots:
                spot = spothash2spot[spot_hash]
                pangenome_name = spots2pangenome[spot_hash]
                node_attributes = spots_graph.nodes[spot_hash]
                del node_attributes[f"{frr_metrics}_cluster"]
                node_attributes.update({"cluster_spots_id": cs_id, "spot_id": spot.ID, "pangenome": pangenome_name})
                cs_spots.add(spot)
            pangenomes.add_conserved_spots(ConservedSpots(cs_id, *cs_spots))
        else:
            spots_graph.remove_nodes_from(cluster_spots)

    return spots_graph


def launch(args):
    """
    Launch functions to align gene families from pangenomes

    Args:
        args: argument given in CLI
    """
    need_info = check_compare_spots_args(args)

    pangenomes, tmpdir, _, lock = common_launch(args, check_pangenome_cs, need_info)

    output = mkdir(args.output, force=args.force)

    spots_graph = compare_spots(pangenomes=pangenomes, dup_margin=args.dup_margin, frr_metrics=args.frr_metrics,
                                frr_cutoff=args.frr_cutoff, threads=args.cpus, lock=lock,
                                disable_bar=args.disable_prog_bar)

    write_conserved_spots(pangenomes, output, cs_graph=spots_graph, graph_formats=args.graph_formats,
                          force=args.force, disable_bar=args.disable_prog_bar)

    if args.systems:
        graph_systems_link_with_conserved_spots(pangenomes=pangenomes, output=output, graph_formats=args.graph_formats,
                                                threads=args.cpus, lock=lock, disable_bar=args.disable_prog_bar)
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
    parser = sub_parser.add_parser("compare_spots",
                                   description='Comparison of systems among pangenomes')

    parser_comparison_spots(parser)
    return parser


def parser_comparison_spots(parser):
    """
    Add argument to parser for system comparison command

    Args:
        parser: parser for cluster argument
    """
    _, compare_opt, optional = parser_comparison(parser)
    compare_opt.add_argument('--frr_metrics', required=False, type=str, default="min_frr",
                             choices=["min_frr", "max_frr"],
                             help="Metrics used to computed spots cluster.")
    optional.add_argument("--dup_margin", required=False, type=float, default=0.05,
                          help="minimum ratio of genomes in which the family must have multiple genes "
                               "for it to be considered 'duplicated'")
    systems = parser.add_argument_group(title="Add systems to conserved spots analyses")
    systems.add_argument('--systems', required=False, action='store_true', default=False,
                         help="Add systems to conserved spots analyses")
    systems.add_argument('-m', '--models', required=False, type=Path, nargs="+", default=None,
                         help="Path to model list file. You can specify multiple models from different source. "
                              "For that separate the model list files by a space and "
                              "make sure you give them in the same order as the sources.")
    systems.add_argument("-s", "--sources", required=False, type=str, nargs="+", default=None,
                         help="Name of the systems sources. You can specify multiple sources. "
                              "For that separate names by a space and "
                              "make sure you give them in the same order as the models.")
    systems.add_argument("--canonical", required=False, action="store_true",
                         help="Write the canonical version of systems too.")
