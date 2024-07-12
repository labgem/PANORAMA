#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from itertools import combinations
from pathlib import Path
from typing import Dict, Set, Tuple
from multiprocessing import Lock
from concurrent.futures import ThreadPoolExecutor

# installed libraries
from tqdm import tqdm
import pandas as pd
import networkx as nx
from ppanggolin.region import Spot
from ppanggolin.RGP.rgp_cluster import compute_grr

# local libraries
from panorama.utils import mkdir, init_lock
from panorama.geneFamily import GeneFamily
from panorama.region import ConservedSpots
from panorama.pangenomes import Pangenomes, Pangenome


def create_pangenome_spots_graph(pangenome: Pangenome,
                                 dup_margin: float = 0.05) -> Tuple[nx.Graph, Dict[Spot, Set[GeneFamily]]]:
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
    multigenic = pangenome.get_multigenics(dup_margin=dup_margin)
    for spot in pangenome.spots:
        graph.add_node(spot)
        spots2borders[spot] = get_borders_families()
    return graph, spots2borders


def create_spots_graph(pangenomes: Pangenomes, dup_margin: float = 0.05, threads: int = 1, lock: Lock = None,
                       disable_bar: bool = False) -> Tuple[nx.Graph, Dict[Spot, Set[GeneFamily]]]:
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
            for future in futures:
                res = future.result()
                spots_graph.add_nodes_from(res[0])
                spots2borders.update(res[1])
    return spots_graph, spots2borders


def compute_grr_edges(graph: nx.Graph, spots2borders: Dict[Spot, Set[GeneFamily]],
                      spots2pangenome: Dict[Spot, str], grr_cutoff: float = 0.8,
                      disable_bar: bool = False):
    """
    Compute spots graph edges with grr score

    Args:
        graph: Spots graph with node only
        spots2borders: Dictionary of spot link to their bordering families
        spots2pangenome: Dictionary of spot to pangenome from which they belong
        grr_cutoff: Grr cutoff for grr score (default: 0.8)
        disable_bar: Flag to disable progress bar (default: False)
    """
    spots_pair = [(spot1, spot2) for spot1, spot2 in combinations(graph.nodes, 2)
                  if spots2pangenome[spot1] != spots2pangenome[spot2]]
    with tqdm(total=len(spots_pair), unit='spots pair', desc="Compute GRR", disable=disable_bar) as pbar:
        for spot1, spot2 in spots_pair:
            border1_refs = {fam.akin.reference for fam in spots2borders[spot1]}
            border2_refs = {fam.akin.reference for fam in spots2borders[spot2]}

            grr = compute_grr(border1_refs, border2_refs, min)
            if grr > grr_cutoff:
                graph.add_edge(spot1, spot2, GRR=grr)
            pbar.update()


def write_conserved_spots(pangenomes, output: Path, force: bool = False, disable_bar: bool = False):
    """
    Write conserved spots into files

    Args:
        pangenomes: Pangenomes associated to conserved spots
        output: Path to the output directory
        force: Flag to overwrite existing files (default: False)
        disable_bar: Flag to disable progress bar (default: False)
    """
    all_cs = []
    cs_dir = mkdir(output/"conserved_spots", force=force, erase=force)
    for conserved_spot in tqdm(pangenomes.conserved_spots, total=pangenomes.number_of_conserved_spots,
                               disable=disable_bar):
        by_cs = []
        for spot in conserved_spot.spots:
            all_cs.append([conserved_spot.ID, spot.ID, spot.pangenome.name, len(spot), spot.number_of_families])
            for rgp in spot.regions:
                by_cs.append([spot.ID, spot.pangenome.name, rgp.name, ",".join([fam.name for fam in rgp.families])])
            by_cs_df = pd.DataFrame(by_cs, columns=['Spot', 'Pangenome', "RGP", "Families"])
            by_cs_df = by_cs_df.sort_values(by=['Spot', 'Pangenome', 'RGP'])
            by_cs_df.to_csv(cs_dir/f"conserved_spots_{conserved_spot.ID}.tsv", sep="\t", header=True, index=False)
    conserved_df = pd.DataFrame(all_cs, columns=['Conserved ID', 'Spot ID', 'Pangenome', "#RGP", "#Families"])
    conserved_df = conserved_df.sort_values(by=['Conserved ID', 'Spot ID', 'Pangenome', '#RGP', "#Families"])
    conserved_df.to_csv(output/"all_conserved_spots.tsv", sep="\t", header=True, index=False)


def identify_conserved_spot(pangenomes: Pangenomes, dup_margin: float = 0.05, threads: int = 1,
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
    spots_graph, spots2borders = create_spots_graph(pangenomes, dup_margin, threads, lock, disable_bar)
    spots2pangenome = {spot: pangenome.name for pangenome in pangenomes for spot in pangenome.spots}
    compute_grr_edges(spots_graph, spots2borders, spots2pangenome, disable_bar)
    cs_id = 0
    for cc_graph in [spots_graph.subgraph(cc).copy() for cc in nx.connected_components(spots_graph)]:
        if len(cc_graph.nodes()) > 1:
            cs_id += 1
            pangenomes.add_conserved_spots(ConservedSpots(cs_id, *cc_graph.nodes()))
