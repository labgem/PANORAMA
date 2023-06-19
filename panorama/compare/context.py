#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from typing import Dict, Tuple
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Lock
import logging
from typing import Dict, Union, List, Set, Iterator
from itertools import combinations
import networkx as nx

from itertools import product
from collections import defaultdict


# installed libraries
from tqdm import tqdm
import pandas as pd
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import search_gene_context_in_pangenome
from ppanggolin.cluster.cluster import read_tsv as read_clustering_table


# local libraries
from panorama.utils import mkdir, init_lock, load_multiple_pangenomes, load_pangenome
from panorama.pangenomes import Pangenome
from panorama.region import GeneContext

# def check_run_context_arguments(kwargs):
#     if "sequences" not in kwargs and "family" not in kwargs:
#         raise Exception("At least one of --sequences or --family option must be given")


# def search_context_mp(pangenome_name: str, pangenome_info: Dict[str, str], output: Path, tmpdir: Path,
#                       threads_per_task: int = 1, **kwargs):
#     pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
#     pangenome.add_file(pangenome_info["path"])
#     search_gene_context_in_pangenome(pangenome=pangenome, output=output, tmpdir=tmpdir,
#                                      cpu=threads_per_task, disable_bar=True, **kwargs)
#     return True

# def parse_cluster_family(cluster_family_file, family_of_interest=None):

#     family_to_cluster = {}
#     with open(cluster_family_file) as fh:
#         for line in fh:
#             cluster_id, familly = line.rstrip().split('\t')

#             family_to_cluster[familly] = cluster_id

#     return family_to_cluster

def parse_context_results(contexts_result_file_list: str) -> Dict[str, Path]:
    """
    Parse the context results file list.

    :param contexts_result_file_list: The path to the file containing the list of pangenome names and context file paths.
    :return: A dictionary mapping pangenome names to context file paths.
    """

    pangenome_to_context_file = {}

    with open(contexts_result_file_list) as fh:
        for line in fh:
            pan_name, context_file = line.rstrip().split('\t')
            pangenome_to_context_file[pan_name] = Path(context_file)

    return pangenome_to_context_file


def make_gene_context_from_context_table(pangenome: Pangenome, context_table: str) -> Set[GeneContext]:
    """
    Create gene contexts from a context table.

    :param pangenome: The Pangenome object.
    :param context_table: The path to the context table.
    :return: A set of GeneContext objects.
    """
    context_objs = set()
    df_context = pd.read_csv(context_table, sep='\t')

    df_context_grp = df_context.groupby(['GeneContext ID']).agg(
        {"Gene family name": set,  "Sequence ID": set})

    for gene_context_id in df_context_grp.index:

        gene_family_names = df_context_grp.loc[gene_context_id,
                                               "Gene family name"]
        
        gene_families = [pangenome.get_gene_family(
            f_name) for f_name in gene_family_names]

        context_obj = GeneContext(
            pangenome, gc_id=gene_context_id, families=gene_families)
        context_objs.add(context_obj)

    return context_objs


def make_gene_context_from_context_graph(pangenome: Pangenome, graph_file: str) -> Set[GeneContext]:
    """
    Create gene contexts from a context graph.

    :param pangenome: The Pangenome object.
    :param context_graph: The path to the context graph.
    :return: A set of GeneContext objects.
    """
    context_objs = set()
    
    if graph_file.suffix == ".graphml":
        contexts_graph  = nx.read_graphml(graph_file)

    elif graph_file.suffix == ".gexf":
        contexts_graph  = nx.read_gexf(graph_file)

    # connected commponents in the graph is a context
    # lets build a context object containing the set of gene families 
    # and graph of the context 
    for i, families_in_context in enumerate(nx.connected_components(contexts_graph)):
        gene_families = [pangenome.get_gene_family(f_name) for f_name in families_in_context]

        context_graph = nx.subgraph_view(contexts_graph, filter_node=lambda n: n in families_in_context).copy()

        
        # node are family id in the current graph. 
        # We may want family object instead to be similar to 
        # what ppanggolin context launch inside panorama would produce

        # gene_family_to_obj = {f_name: pangenome.get_gene_family(f_name) for f_name in families_in_context}
        # G = nx.Graph()
        # G.add_edges_from(((gene_family_to_obj[f1_name], gene_family_to_obj[f2_name]) for f1_name,f2_name in context_graph.edges()))


        context_obj = GeneContext(pangenome, gc_id=i, families=gene_families, graph=context_graph)
        context_objs.add(context_obj)
        
    return context_objs


def write_context_summary(gene_contexts: List[GeneContext], output_table: Path):
    """
    Write a summary of gene contexts to a table.

    :param gene_contexts: A list of GeneContext objects representing gene contexts to summarize.
    :param output_table: The path to the output table file where the summary will be written.
    """

    gene_context_summaries = [gc.summarize() for gc in gene_contexts]
    summary_df = pd.DataFrame(gene_context_summaries)

    summary_df.to_csv(output_table, sep='\t', index=False)


def load_pangenome_and_get_contexts_from_result(context_result_path: Path, pangenome_name: str, pangenome_path: Path, taxid: str) -> List[GeneContext]:
    """
    Retrieve gene contexts from a table and create GeneContext objects.

    :param context_result_path: The path to the context table file or graph.
    :param pangenome_name: The name of the pangenome.
    :param pangenome_path: The path to the pangenome file.
    :param taxid: The taxonomic ID associated with the pangenome.
    """

    pangenome = load_pangenome(pangenome_name, pangenome_path, taxid, {
                               "need_families": True, })

    if context_result_path.suffix == ".tsv":
        gene_contexts = make_gene_context_from_context_table(pangenome, context_result_path)
    elif context_result_path.suffix in [".graphml", ".gexf"]:
        gene_contexts = make_gene_context_from_context_graph(pangenome, context_result_path)
    else:
        # TODO File extension should be checked when parsing file list. 
        raise ValueError('The gene context result has not the correct extension. {}'
                         'Panorama expects "tsv" for context tables and "graphml" or "gexf" for graph contexts.')
    return gene_contexts


def get_gene_contexts_from_results_mp(pan_name_to_path: Dict[str, Dict[str, Union[str, int]]],
                                     pan_name_to_context_result: Dict[str, Path],
                                     max_workers: int, disable_bar: bool) -> List[Pangenome]:
    """
    Retrieve gene contexts from multiple result files using multiprocessing.

    :param pan_name_to_path: A dictionary mapping pangenome names to their path information.
    :param pan_name_to_context_result: A dictionary mapping pangenome names to their corresponding context tables or graphs.
    :param max_workers: The maximum number of workers to use for multiprocessing.
    :param disable_bar: A boolean value indicating whether to disable the progress bar.
    :return: A list of Pangenome objects containing the retrieved gene contexts.
    """

    with ProcessPoolExecutor(max_workers=max_workers) as executor:

        futures = []

        for pangenome_name, pangenome_path_info in pan_name_to_path.items():

            context_result = pan_name_to_context_result[pangenome_name]

            future = executor.submit(load_pangenome_and_get_contexts_from_result, context_result, pangenome_name, pangenome_path_info["path"],
                                     pangenome_path_info["taxid"])
            futures.append(future)

        gene_contexts = [gene_context for future in tqdm(
            futures, unit="pangenome", disable=disable_bar) for gene_context in future.result()]

    return gene_contexts


def compare_pair_of_contexts(context_pair: Tuple[GeneContext, GeneContext], min_jaccard:float) -> Tuple[GeneContext, GeneContext, float]:
    """
    Compares a pair of gene contexts and calculates the Jaccard similarity between their family clusters.

    :param context_pair: A tuple containing two GeneContext objects to be compared.
    :param min_jaccard: min jaccard cutoff to report a pair of contexts.
    :return: A tuple containing the two GeneContext objects and the Jaccard similarity between their family clusters.
    """

    contextA, contextB = context_pair
    contextA_clst_family =  {gf.family_cluster for gf in contextA.families}
    contextB_clst_family = {gf.family_cluster for gf in contextB.families}
    shared_family = len(contextA_clst_family & contextB_clst_family)

    clst_family_jaccard = shared_family / len(contextA_clst_family | contextB_clst_family)
    if clst_family_jaccard >= min_jaccard: 
        return contextA.ID, contextB.ID, {'jaccard':clst_family_jaccard, 
                                          "shared family":shared_family,
                                          "jaccard_edge":True}
    else: 
        return contextA.ID, contextB.ID, None


def create_metanodes(gfA_to_cf: dict, gfB_to_cf: dict) -> Tuple[List[Tuple[str, dict]], dict, dict]:
    """
    Create metanodes for a multigraph based on gene family mappings.

    :param gfA_to_cf: A dictionary mapping gene families in graph A to their cluster families.
    :param gfB_to_cf: A dictionary mapping gene families in graph B to their cluster families.
    :return: A tuple containing the metanodes, graph A node to metanodes mapping, and graph B node to metanodes mapping.

    Metanodes are created for gene families that have a common cluster family between the two graphs.

    When a cluster family is associated with more than one gene family in a graph, multiple metanodes are created,
    and each metanode is differentiated by adding "´" to its name.
    """

    clusters = set(gfB_to_cf.values()) | set(gfA_to_cf.values())

    meta_nodes = {}

    gA_node_2_meta_nodes = defaultdict(list)
    gB_node_2_meta_nodes = defaultdict(list)

    for cluster in clusters:
        clstr_gA_nodes = [n for n in gfA_to_cf if gfA_to_cf[n] == cluster]
        clstr_gB_nodes = [n for n in gfB_to_cf if gfB_to_cf[n] == cluster]

        for i, (n_gA, n_gB) in enumerate(product(clstr_gA_nodes, clstr_gB_nodes)):
            full_name_meta_node = f"{cluster}: ({n_gB}, {n_gA})"

            meta_node = f"{cluster}" if i == 0 else f"{cluster}{'´'*i}"

            meta_node_attr =  {"cluster":cluster, "node_gA":n_gA, "node_gB":n_gB, 'full_name':full_name_meta_node}
            meta_nodes[meta_node] = meta_node_attr
            gA_node_2_meta_nodes[n_gA].append(meta_node)
            gB_node_2_meta_nodes[n_gB].append(meta_node)
            
    return meta_nodes, gA_node_2_meta_nodes, gB_node_2_meta_nodes
            


def get_multigraph_edges(g: nx.Graph, g_node_2_meta_nodes: dict) -> List[Tuple[str, str]]:
    """
    Translate edges of a graph into edges linking metanodes of a multigraph.

    This function takes a graph and a mapping of graph nodes to their corresponding metanodes
    and translates the edges of the graph into metanodes of the multigraph.
    Only nodes that have metanodes are considered, and the edges are formed by combining the metanodes
    corresponding to the endpoints of each edge.

    :param g: The graph whose edges are to be translated.
    :param g_node_2_meta_nodes: A dictionary mapping graph nodes to their corresponding metanodes.
    :return: A list of edges in the multigraph.
    """
    g_multigraph_edges = []
    for n1, n2 in g.edges():
        if n1 in g_node_2_meta_nodes and n2 in g_node_2_meta_nodes:
            n1_meta_nodes = g_node_2_meta_nodes[n1]
            n2_meta_nodes = g_node_2_meta_nodes[n2]

            g_multigraph_edges += list(product(n1_meta_nodes, n2_meta_nodes))
    
    return g_multigraph_edges




def get_connected_components(nodes: List[str], edges: List[Tuple[str, str]]) -> Iterator[Set[str]]:
    """
    Get the connected components in a graph.

    :param nodes: List of nodes in the graph.
    :param edges: List of edges in the graph.
    :return: Iterator of sets representing the connected components.
    """
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    return nx.connected_components(G)


def compute_CCC(meta_nodes: List[str], g1_edges: List[Tuple[str, str]], 
                                           g2_edges: List[Tuple[str, str]]) -> List[Set[str]]:
    """
    Compute the conserved connected components (CCC) between two graphs. 
    The method implement here is taken from Boyer et al. article
    (https://doi.org/10.1093/bioinformatics/bti711).

    :param meta_nodes: List of meta-nodes representing the common family clusters.
    :param g1_edges: List of edges in graph 1.
    :param g2_edges: List of edges in graph 2.
    :return: List of sets representing the conserved connected components.
    """

    g1_cc = get_connected_components(meta_nodes, g1_edges)
    g2_cc = get_connected_components(meta_nodes, g2_edges)

    # Compute the intersection of all connected components
    cc_intersections = [cc1 & cc2 for cc1, cc2 in product(g1_cc, g2_cc)]

    # If only one intersection, return it
    if len(cc_intersections) == 1:
        return cc_intersections
    
    partitions = []
    for i, cc_inter in enumerate(cc_intersections):
        # Recursively compute connected components within the intersection
        g1_edges_inter = [(n, v) for n, v in g1_edges if n in cc_inter and v in cc_inter]
        g2_edges_inter = [(n, v) for n, v in g2_edges if n in cc_inter and v in cc_inter]
        partitions += compute_CCC(cc_inter, g1_edges_inter, g2_edges_inter)

    return partitions

def get_conserved_genomics_contexts(gcA_graph: nx.Graph, gcB_graph: nx.Graph,
                                    gene_fam_2_cluster_fam: Dict[str, str],
                                    min_cgc_size:int =2,
                                    return_multigraph:bool = False) -> List[Tuple[Set[str], Set[str]]]:
    """
    Get the conserved genomics contexts between two gene context graphs.

    :param gcA_graph: The gene context graph A.
    :param gcB_graph: The gene context graph B.
    :param gene_fam_2_cluster_fam: Dictionary mapping gene families to cluster families.
    :param min_cgc_size: Minimum size of a conserved genomic context to report it.
    :param return_multigraph: Flag indicating whether to return the multigraph representation.

    :return: Tuple containing a list of tuples representing the conserved genomics contexts and
             the multigraph representation (if return_multigraph is True).
             Each tuple in the list contains two sets: the gene nodes from graph A and the gene nodes from graph B.
    """


    multigraph = None

    gfA_to_cf = {gf:cf for gf, cf in gene_fam_2_cluster_fam.items() if gf in gcA_graph}
    gfB_to_cf = {gf:cf for gf, cf in gene_fam_2_cluster_fam.items() if gf in gcB_graph}

    if len(set(gfA_to_cf.values()) & set(gfB_to_cf.values())) < min_cgc_size:
        # in case the two graph share not enough cluster to reach minimum context size
        return [], None

    meta_nodes_2_attributes, gA_node_2_meta_nodes, gB_node_2_meta_nodes = create_metanodes(gfA_to_cf, gfB_to_cf)

    gA_multigraph_edges = get_multigraph_edges(gcA_graph, gA_node_2_meta_nodes)
    gB_multigraph_edges = get_multigraph_edges(gcB_graph, gB_node_2_meta_nodes)

    partitions = compute_CCC(meta_nodes_2_attributes.keys(), gA_multigraph_edges, gB_multigraph_edges)

    cgc_nodes = []

    for i, meta_nodes in enumerate(partitions):
        if len(meta_nodes) >= min_cgc_size:
            gB_nodes = {meta_nodes_2_attributes[meta_node]['node_gB']  for meta_node in meta_nodes}
            gA_nodes = {meta_nodes_2_attributes[meta_node]['node_gA']  for meta_node in meta_nodes}
            cgc_nodes.append((gA_nodes, gB_nodes))
            
            for meta_node in meta_nodes:
                meta_nodes_2_attributes[meta_node]["GCG"] = f"CGC_{i}"
    
    # Construct the multigraph if requested. 
    # This graph is used for visualization and verification.
    if return_multigraph:
        multigraph = nx.MultiGraph()
        multigraph.add_nodes_from(meta_nodes)

        multigraph.add_edges_from(gA_multigraph_edges, origin="graphA")
        multigraph.add_edges_from(gB_multigraph_edges, origin="graphB")


    return cgc_nodes, multigraph


def compare_pair_of_context_graphs(context_pair: Tuple[GeneContext, GeneContext], 
                                   gene_fam_2_cluster_fam: Dict[str,str], 
                                   return_multigraph:bool):

    contextA, contextB = context_pair
    # Get conserved genomic context from the two context graph 
    conserved_genomics_contexts, multigraph = get_conserved_genomics_contexts(contextA.graph, contextB.graph, gene_fam_2_cluster_fam, return_multigraph=return_multigraph)

    cgc_sizes = [max((len(nodesA_in_cgc), len(nodesB_in_cgc))) for nodesA_in_cgc, nodesB_in_cgc in conserved_genomics_contexts]

    # TODO: Score the CGCs

    # Compute simple metrics 
    contextA_clst_family =  {gf.family_cluster for gf in contextA.families}
    contextB_clst_family = {gf.family_cluster for gf in contextB.families}
    shared_family = len(contextA_clst_family & contextB_clst_family)
    clst_family_jaccard = shared_family / len(contextA_clst_family | contextB_clst_family)


    if clst_family_jaccard >= 0.5:
        return contextA.ID, contextB.ID, {'family_compo_jaccard':clst_family_jaccard,
                                          "shared family":shared_family,
                                          "cgc_count":len(conserved_genomics_contexts),
                                          "cgc_sizes":':'.join([str(size) for size in sorted(cgc_sizes)]),}, multigraph
    else:
        return contextA.ID, contextB.ID, None, multigraph

def launch_context_comparison(pack: tuple) -> tuple:
    """ 
    Allow to launch in multiprocessing the context comparison

    :param pack: Pack of argument for context comparison

    :return: edge metrics 
    """
    return compare_pair_of_contexts(*pack)

def launch_compare_pair_of_context_graphs(pack: tuple) -> tuple:
    """ 
    Allow to launch in multiprocessing the context comparison

    :param pack: Pack of argument for context comparison

    :return: edge metrics 
    """
    return compare_pair_of_context_graphs(*pack)



def compare_gene_contexts_on_cluster_families(gene_contexts: List[GeneContext], min_jaccard, max_workers: int, disable_bar: bool) -> List[GeneContext]:
    """
    Compares gene contexts by calculating the Jaccard similarity between their family clusters.

    :param gene_contexts: A list of GeneContext objects to be compared.
    :param min_jaccard: jaccard cutoff to report a pair of contexts.
    :param max_workers: The maximum number of worker processes for parallel execution.
    :param disable_bar: A boolean flag indicating whether to disable the progress bar.
    :return: A list of GeneContext objects.
    """
    context_graph = nx.Graph()
    for gc in gene_contexts:
        context_graph.add_node(gc.ID, pangenome=gc.pangenome)

    context_pairs =  combinations(gene_contexts, 2)
    comparison_arguments = ((p, min_jaccard) for p in context_pairs)
    pair_count = (len(gene_contexts) **2 - len(gene_contexts)) / 2

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for gcA, gcB, metrics in tqdm(executor.map(launch_context_comparison, comparison_arguments, chunksize=5), total=pair_count, 
                                      disable=disable_bar, unit="context pair"):
            if metrics:
                context_graph.add_edge( gcA, gcB, **metrics)

    logging.info(f'Context graph: {context_graph}')
    return context_graph

def compare_gene_contexts_graph_mp(gene_contexts: List[GeneContext], 
                                   gene_fam_2_cluster_fam: Dict[str, str], 
                                   max_workers: int, disable_bar: bool,
                                   return_multigraph: bool) -> List[GeneContext]:
    """
    Compares gene contexts by looking at their context graphs.

    :param gene_contexts: A list of GeneContext objects to be compared.
    :param gene_fam_2_cluster_fam: dict mapping gene family name to cluster family name
    :param max_workers: The maximum number of worker processes for parallel execution.
    :param disable_bar: A boolean flag indicating whether to disable the progress bar.
    :return: A list of GeneContext objects.
    """
    context_graph = nx.Graph()
    multigraphs = {}
    for gc in gene_contexts:
        context_graph.add_node(gc.ID, pangenome=gc.pangenome)

    context_pairs =  combinations(gene_contexts, 2)
    comparison_arguments = ((p, gene_fam_2_cluster_fam, return_multigraph) for p in context_pairs)
    pair_count = (len(gene_contexts) **2 - len(gene_contexts)) / 2

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for gcA, gcB, metrics, multigraph in tqdm(executor.map(launch_compare_pair_of_context_graphs, comparison_arguments, chunksize=5), total=pair_count, 
                                      disable=disable_bar, unit="context pair"):
            if metrics:
                context_graph.add_edge( gcA, gcB, **metrics)

                if multigraph:
                    multigraphs[(gcA, gcB)] = multigraph

    logging.info(f'Context graph: {context_graph}')
    return context_graph, multigraphs


def context_comparison(pangenome_to_path: Dict[str, Union[str, int]], contexts_results: Path, family_clusters: bool, min_jaccard,
                       lock: Lock, output: Path, tmpdir: Path, task: int = 1, threads_per_task: int = 1,
                       disable_bar: bool = False, force: bool = False, **kwargs):
    """
    Perform comparison of gene contexts and cluster families.

    :param pangenome_to_path: A dictionary mapping pangenome names to their corresponding paths.
    :param contexts_results: The path to the file containing the list of context results.
    :param family_clusters: A boolean indicating whether to use precomputed family clusters or run the cluster family function.
    :param lock: A Lock object for thread synchronization.
    :param output: The output directory path.
    :param tmpdir: The temporary directory path.
    :param task: The number of tasks for parallel processing (default: 1).
    :param threads_per_task: The number of threads per task (default: 1).
    :param disable_bar: A boolean indicating whether to disable progress bars (default: False).
    :param force: A boolean indicating whether to force overwriting existing output files (default: False).
    :param kwargs: Additional keyword arguments.
    """

    mkdir(output, force)

    if contexts_results:
        logging.info(f"Retrieving gene contexts from existing results: {contexts_results}")

        # TODO check consistency between pangenome_to_path and context results

        pan_name_to_context_file = parse_context_results(contexts_results)

        gene_contexts = get_gene_contexts_from_results_mp(
            pangenome_to_path, pan_name_to_context_file, task, disable_bar)
        

    else:
        # run ppanggolin context in parallel
        raise NotImplementedError

        check_run_context_arguments(kwargs)

        with ProcessPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
            with tqdm(total=len(pangenome_to_path), unit='pangenome', disable=disable_bar) as progress:
                futures = []
                for pangenome_name, pangenome_info in pangenome_to_path.items():
                    mkdir(output/pangenome_name, force)
                    future = executor.submit(search_context_mp, pangenome_name, pangenome_info, output, tmpdir,
                                             threads_per_task, **kwargs)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)
                for future in futures:
                    results = future.result()

    # write gene context summary
    summary_out_table = output / "gene_context_summary.tsv"
    logging.info(f'Writting gene context summary: {summary_out_table} ')
    write_context_summary(gene_contexts, summary_out_table)

    if family_clusters:
        logging.info(f"Retrieving family clusters from existing clustering: {family_clusters}")
        
        # Parse the given cluster family results
        gene_family_to_family_cluster, cluster2family = read_clustering_table(family_clusters)
        # remove fragmentation info
        gene_family_to_family_cluster = {gf:fc for gf, (fc, is_fragmented) in gene_family_to_family_cluster.items()} 
        
    else:
        # run cluster familly
        # family_clusters_file = panorama cluster function
        raise NotImplementedError

    # TODO Check that all gene families in context have a family cluster

    # add family cluster info in gene contexts 
    for gene_context in gene_contexts:
        for gene_family in gene_context.families:
            family_cluster = gene_family_to_family_cluster[gene_family.name]
            gene_family.add_family_cluster(family_cluster)

    # Compare gene contexts based on their family clusters  

    context_graph_clstr_families = compare_gene_contexts_on_cluster_families(gene_contexts, min_jaccard, task, disable_bar)


    # Compare gene contexts based on their synteny information  
    context_synteny_graph, context_pair_2_multigraphs = compare_gene_contexts_graph_mp(gene_contexts, gene_family_to_family_cluster, task, disable_bar, 
                                                                        return_multigraph=True)
    
    context_graph_merged = nx.compose(context_graph_clstr_families, context_synteny_graph)


    # add node attributes

    nx.set_node_attributes(context_graph_merged, {gc.ID:gc.summarize() for gc in gene_contexts})


    context_graph_file = output / f"context.graphml"
    logging.info(f'Writting gene context graph: {context_graph_file}')
    nx.readwrite.graphml.write_graphml(context_graph_merged, context_graph_file)

    # Writting multigraphs
    for (contextA, contextB), multigraph in context_pair_2_multigraphs.items():
        logging.debug(f"{(contextA, contextB)}, {multigraph}")
        multigraph_file = output / f"multigraph_{contextA}_vs_{contextB}.graphml"

        nx.readwrite.graphml.write_graphml(multigraph, multigraph_file)


def context_comparison_parser(parser):

    use_context_arg = parser.add_argument_group(
        "Use already computed contexts arguments")

    use_context_arg.add_argument('--context_results', type=Path, required=False,
                                 help="Tsv file with two columns: name of pangenome and path to the corresponding context results."
                                 "Results can be a table (tsv) or a graph (graphml or gexf)")
    use_context_arg.add_argument('--family_clusters', type=Path, required=False,
                                 help="A tab-separated file listing the cluster names, the family IDs,")
    
    use_context_arg.add_argument('--min_jaccard', type=restricted_float, required=False, default=0.5,
                                 help="Minimum value for jaccard index")
    
    # TODO: use context parser of ppanggolin to have correct args and prevent duplication

    # run_context_arg = parser.add_argument_group("PPanGGOLiN context arguments")

    # # required = parser._action_groups[2]
    # run_context_arg.add_argument('-S', '--sequences', type=Path, required=False,
    #                       help="Fasta file with the sequences of interest")

    # run_context_arg.add_argument('-F', '--family', type=Path, required=False,
    #                       help="List of family IDs of interest from the pan")

    # # optional = parser._action_groups[4]

    # run_context_arg.add_argument('--no_defrag', required=False, action="store_true",
    #                       help="DO NOT Realign gene families to link fragments with"
    #                            "their non-fragmented gene family.")

    # run_context_arg.add_argument('--identity', required=False, type=float, default=0.5,
    #                       help="min identity percentage threshold")

    # run_context_arg.add_argument('--coverage', required=False, type=float, default=0.8,
    #                       help="min coverage percentage threshold")

    # run_context_arg.add_argument("-t", "--transitive", required=False, type=int, default=4,
    #                       help="Size of the transitive closure used to build the graph. This indicates the number of "
    #                            "non related genes allowed in-between two related genes. Increasing it will improve "
    #                            "precision but lower sensitivity a little.")

    # run_context_arg.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
    #                       help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
    #                            "will improve precision but lower sensitivity a lot.")
