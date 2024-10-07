#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
from typing import Dict, Tuple, Union
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Lock
import logging
from typing import Dict, Union, List, Set, Iterator
from itertools import combinations
import networkx as nx
from shutil import rmtree
from itertools import product
from collections import defaultdict

# installed libraries
from tqdm import tqdm
import pandas as pd
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import search_gene_context_in_pangenome, check_pangenome_for_context_search
from ppanggolin.cluster.cluster import read_tsv as read_clustering_table

# local libraries
from panorama.utils import mkdir
from panorama.format.read_binaries import load_pangenome
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.region import GeneContext
from panorama.compare.utils import parser_comparison, common_launch


def check_context_comparison(args):
    """
    Checks the provided keyword arguments to ensure either 'sequences' or 'family' is present.

    Args:
        kwargs (dict): Keyword arguments to check.

    Raises:
        Exception: If neither 'sequences' nor 'family' is present in kwargs.
    """
    if args.context_results:
        if args.identity > 1 or args.coverage > 1:
            raise argparse.ArgumentError(message="Identity and coverage must be between 0 and 1", argument=None)
        if args.transitive < 0 or args.window < 0:
            raise argparse.ArgumentError(message="Transitivity and window must be positif", argument=None)
        if args.jaccard > 1:
            raise argparse.ArgumentError(message="Jaccard must be between 0 and 1", argument=None)


def launch_ppanggolin_context(pangenome_to_path: dict, ppanggolin_context_args: dict, output, tmpdir, max_workers: int,
                              threads_per_task: int, disable_bar):
    pangenome_name_to_graph_outfile = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:

        futures = {}

        for pangenome_name, pangenome_info in pangenome_to_path.items():
            future = executor.submit(search_context_mp, pangenome_name, pangenome_info, ppanggolin_context_args,
                                     output=output, tmpdir=tmpdir,
                                     threads_per_task=threads_per_task)
            futures[pangenome_name] = future

        all_gene_contexts = set()
        for pangenome_name, future in tqdm(futures.items(), unit="pangenome", disable=disable_bar):
            gene_contexts, graph_outfile = future.result()
            all_gene_contexts |= gene_contexts
            pangenome_name_to_graph_outfile[pangenome_name] = graph_outfile

    # Write tsv listing graph file to be able to rerun compare context without reruning ppanggolin  
    graph_list_file = output / "context_graph_files.tsv"
    logging.info(f'Writting list of context graph files in {graph_list_file}')

    with open(graph_list_file, "w") as fl:
        file_content = ''.join((f"{pangenome_name}\t{graph_outfile}\n" for pangenome_name, graph_outfile in
                                pangenome_name_to_graph_outfile.items()))
        fl.write(file_content)

    return all_gene_contexts


def search_context_mp(pangenome_name: str, pangenome_info: Dict[str, str], ppanggolin_context_args: dict, output: Path,
                      tmpdir: Path,
                      threads_per_task: int = 1, ):
    # output has been already created using the force flag.
    output_for_current_pan = output / pangenome_name
    mkdir(output_for_current_pan, force=True)

    pangenome = load_pangenome(pangenome_name, pangenome_info["path"], pangenome_info["taxid"], {
        "need_families": True, 'need_annotations': True, }, disable_bar=True)

    gene_context_graph, graph_outfile = search_gene_context_in_pangenome(pangenome=pangenome,
                                                                         output=str(output_for_current_pan),
                                                                         tmpdir=tmpdir,
                                                                         cpu=threads_per_task, disable_bar=True,
                                                                         **ppanggolin_context_args)

    gene_contexts = make_gene_context_from_context_graph(pangenome, gene_context_graph)

    return gene_contexts, graph_outfile


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
        {"Gene family name": set, "Sequence ID": set})

    for gene_context_id in df_context_grp.index:
        gene_family_names = df_context_grp.loc[gene_context_id,
        "Gene family name"]

        gene_families = [pangenome.get_gene_family(
            f_name) for f_name in gene_family_names]

        context_obj = GeneContext(
            pangenome, gc_id=gene_context_id, families=gene_families)
        context_objs.add(context_obj)

    return context_objs


def clean_context_objects(contexts):
    """
    Remove some of the object contain in context to make it pickable. 
    """

    families = {gf for gc in contexts for gf in gc.families}
    for gf in families:
        gf.genes = set()
        gf._genePerOrg = dict()


def make_gene_context_from_context_graph(pangenome: Pangenome, contexts_graph: nx.Graph) -> Set[GeneContext]:
    """
    Create gene contexts from a context graph.

    :param pangenome: The Pangenome object.
    :param context_graph: The context graph
    :return: A set of GeneContext objects.
    """
    context_objs = set()
    families_of_interest = {pangenome.get_gene_family(f_name) for f_name, d in contexts_graph.nodes(data=True) if
                            d["families_of_interest"]}

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

        context_obj = GeneContext(pangenome, gc_id=i, families=gene_families, families_of_interest=families_of_interest)

        context_obj.add_context_graph(context_graph)

        context_objs.add(context_obj)

    clean_context_objects(context_objs)

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


def load_pangenome_and_get_contexts_from_result(context_result_path: Path, pangenome_name: str, pangenome_path: Path,
                                                taxid: str) -> List[GeneContext]:
    """
    Retrieve gene contexts from a table and create GeneContext objects.

    :param context_result_path: The path to the context table file or graph.
    :param pangenome_name: The name of the pangenome.
    :param pangenome_path: The path to the pangenome file.
    :param taxid: The taxonomic ID associated with the pangenome.
    """

    pangenome = load_pangenome(pangenome_name, pangenome_path, taxid, {
        "need_families": True, })

    # if context_result_path.suffix == ".tsv":
    #     gene_contexts = make_gene_context_from_context_table(pangenome, context_result_path)
    if context_result_path.suffix in [".graphml", ".gexf"]:

        if context_result_path.suffix == ".graphml":
            contexts_graph = nx.read_graphml(context_result_path)

        elif context_result_path.suffix == ".gexf":
            contexts_graph = nx.read_gexf(context_result_path)

        gene_contexts = make_gene_context_from_context_graph(pangenome, contexts_graph)
    else:
        # TODO File extension should be checked when parsing file list. 
        raise ValueError(f'The gene context result has not the correct extension. {context_result_path}'
                         'Panorama expects "graphml" or "gexf".')
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

            future = executor.submit(load_pangenome_and_get_contexts_from_result, context_result, pangenome_name,
                                     pangenome_path_info["path"],
                                     pangenome_path_info["taxid"])
            futures.append(future)

        gene_contexts = [gene_context for future in tqdm(
            futures, unit="pangenome", disable=disable_bar) for gene_context in future.result()]

    return gene_contexts


def compare_pair_of_contexts(context_pair: Tuple[GeneContext, GeneContext], min_jaccard: float) -> Tuple[
    GeneContext, GeneContext, float]:
    """
    Compares a pair of gene contexts and calculates the Jaccard similarity between their family clusters.

    :param context_pair: A tuple containing two GeneContext objects to be compared.
    :param min_jaccard: min jaccard cutoff to report a pair of contexts.
    :return: A tuple containing the two GeneContext objects and the Jaccard similarity between their family clusters.
    """

    contextA, contextB = context_pair
    contextA_clst_family = {gf.akin for gf in contextA.families}
    contextB_clst_family = {gf.akin for gf in contextB.families}
    shared_family = len(contextA_clst_family & contextB_clst_family)

    clst_family_jaccard = shared_family / len(contextA_clst_family | contextB_clst_family)
    if clst_family_jaccard >= min_jaccard:
        return contextA.ID, contextB.ID, {'jaccard': clst_family_jaccard,
                                          "shared family": shared_family,
                                          "jaccard_edge": True}
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

            meta_node = f"{cluster}" if i == 0 else f"{cluster}{'´' * i}"

            meta_node_attr = {"cluster": cluster, "node_gA": n_gA, "node_gB": n_gB, 'full_name': full_name_meta_node}
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
    g_multigraph_edges_2_attr = {}
    for n1, n2, d in g.edges(data=True):
        if n1 in g_node_2_meta_nodes and n2 in g_node_2_meta_nodes:
            n1_meta_nodes = g_node_2_meta_nodes[n1]
            n2_meta_nodes = g_node_2_meta_nodes[n2]

            g_multigraph_edges = list(product(n1_meta_nodes, n2_meta_nodes))

            g_multigraph_edges_2_attr.update({(n, v): dict(d) for n, v in g_multigraph_edges})

    return g_multigraph_edges_2_attr


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


def get_shortest_path_edges_cc_strategy(g, weight="mean_transitivity"):
    """
    TODO : This is an aproximation of the shortest path algo. 

    A more robust implementation is required. 
    https://www.baeldung.com/cs/shortest-path-visiting-all-nodes

    """
    G = g.copy()

    intial_number_cc = nx.number_connected_components(G)

    sum_edges = 0
    for u, v, d in sorted(G.edges(data=True), key=lambda n_v_d: n_v_d[2][weight], reverse=True):

        G.remove_edge(u, v)

        cc_count = nx.number_connected_components(G)

        if cc_count != intial_number_cc:
            G.add_edge(u, v, **d)
            sum_edges += d[weight]

    return G.edges, sum_edges


def get_conserved_genomics_contexts(gcA_graph: nx.Graph, gcB_graph: nx.Graph,
                                    gene_fam_2_cluster_fam: Dict[str, str],
                                    min_cgc_size: int = 2,
                                    return_multigraph: bool = False) -> List[Tuple[Set[str], Set[str]]]:
    """
    Get the conserved genomics contexts between two gene context graphs.

    :param gcA_graph: The gene context graph A.
    :param gcB_graph: The gene context graph B.
    :param gene_fam_2_cluster_fam: Dictionary mapping gene families to cluster families.
    :param min_cgc_size: Minimum size of a conserved genomic context to report.
    :param return_multigraph: Flag indicating whether to return the multigraph representation.

    :return: Tuple containing a list of tuples representing the conserved genomics contexts and
             the multigraph representation (if return_multigraph is True).
             Each tuple in the list contains two sets: the gene nodes from graph A and the gene nodes from graph B.
    """

    # TODO clean this copy 
    # this uselful when investigate intermediate graph to not have noise from other context comparison in the context graph attributes... 
    gcA_graph = gcA_graph.copy()
    gcB_graph = gcB_graph.copy()

    multigraph = None

    gfA_to_cf = {gf: cf for gf, cf in gene_fam_2_cluster_fam.items() if gf in gcA_graph}
    gfB_to_cf = {gf: cf for gf, cf in gene_fam_2_cluster_fam.items() if gf in gcB_graph}

    # add cluster family in A and B graphs. Only useful when return_multigraph is True
    nx.set_node_attributes(gcA_graph, gfA_to_cf, name="cluster")
    nx.set_node_attributes(gcB_graph, gfB_to_cf, name="cluster")

    if len(set(gfA_to_cf.values()) & set(gfB_to_cf.values())) < min_cgc_size:
        # in case the two graph share not enough cluster to reach minimum context size
        return [], None

    meta_nodes_2_attributes, gA_node_2_meta_nodes, gB_node_2_meta_nodes = create_metanodes(gfA_to_cf, gfB_to_cf)

    # add metanode family in A and B graphs. Only useful when return_multigraph is True
    nx.set_node_attributes(gcA_graph, {n: len(mn) > 0 for n, mn in gA_node_2_meta_nodes.items()}, name="is_metanode")
    nx.set_node_attributes(gcB_graph, {n: len(mn) > 0 for n, mn in gB_node_2_meta_nodes.items()}, name="is_metanode")

    gA_multigraph_edges_2_attr = get_multigraph_edges(gcA_graph, gA_node_2_meta_nodes)
    gB_multigraph_edges_2_attr = get_multigraph_edges(gcB_graph, gB_node_2_meta_nodes)

    conserved_genomic_contexts = compute_CCC(meta_nodes_2_attributes.keys(), gA_multigraph_edges_2_attr.keys(),
                                             gB_multigraph_edges_2_attr.keys())

    cgc_infos = []
    cgc_graphs = []

    # filter cgc that does not match required size
    conserved_genomic_contexts_filtered = (cgc_nodes for cgc_nodes in conserved_genomic_contexts if
                                           len(cgc_nodes) >= min_cgc_size)

    for i, meta_nodes in enumerate(conserved_genomic_contexts_filtered):

        gB_nodes = {meta_nodes_2_attributes[meta_node]['node_gB'] for meta_node in meta_nodes}
        gA_nodes = {meta_nodes_2_attributes[meta_node]['node_gA'] for meta_node in meta_nodes}

        # filter graph A and graph B to include only nodes included in the current conserved_genomic_context 
        graphA_cgc = nx.subgraph_view(gcA_graph, filter_node=lambda x: x in gA_nodes)
        graphB_cgc = nx.subgraph_view(gcB_graph, filter_node=lambda x: x in gB_nodes)

        shortest_path_edges_A, sum_transitivity_edges_A = get_shortest_path_edges_cc_strategy(graphA_cgc,
                                                                                              weight="mean_transitivity")
        shortest_path_edges_B, sum_transitivity_edges_B = get_shortest_path_edges_cc_strategy(graphB_cgc,
                                                                                              weight="mean_transitivity")

        shortest_path_attributes_A = {(n, v): (n, v) in shortest_path_edges_A for n, v in graphA_cgc.edges()}
        shortest_path_attributes_B = {(n, v): (n, v) in shortest_path_edges_B for n, v in graphB_cgc.edges()}

        nx.set_edge_attributes(gcA_graph, shortest_path_attributes_A, name="in shortest path")
        nx.set_edge_attributes(gcB_graph, shortest_path_attributes_B, name="in shortest path")

        nx.set_edge_attributes(graphA_cgc, shortest_path_attributes_A, name="in shortest path")
        nx.set_edge_attributes(graphB_cgc, shortest_path_attributes_B, name="in shortest path")

        cgc_graphs.append(nx.union(graphA_cgc, graphB_cgc))

        cgc_score = (len(gA_nodes) / (1 + sum_transitivity_edges_A) + len(gB_nodes) / (
                1 + sum_transitivity_edges_B)) / 2
        cgc_mean_size_in_both_graph = (len(gA_nodes) + len(gB_nodes)) / 2

        for meta_node in meta_nodes:
            meta_nodes_2_attributes[meta_node]["CGC"] = f"CGC_{i}"
            meta_nodes_2_attributes[meta_node]["CGC_score"] = cgc_score
            meta_nodes_2_attributes[meta_node]["cgc_mean_size_in_both_graph"] = cgc_mean_size_in_both_graph

        cgc_info = {"cgc_id": i,
                    'graphA_nodes': gA_nodes,
                    "graphB_nodes": gB_nodes,
                    "score": cgc_score,
                    "mean_size": (len(gA_nodes) + len(gB_nodes)) / 2
                    }

        cgc_infos.append(cgc_info)

    # Construct the multigraph if requested. 
    # This graph is used for visualization and verification.
    if return_multigraph:
        pass_graph_attribute_to_multigraph(meta_nodes_2_attributes, gcA_graph, node_mapper="node_gA")
        pass_graph_attribute_to_multigraph(meta_nodes_2_attributes, gcB_graph, node_mapper="node_gB")

        multigraph = nx.MultiGraph()
        multigraph.add_nodes_from([(mn, attr) for mn, attr in meta_nodes_2_attributes.items()])

        for (n, v), attributes in gB_multigraph_edges_2_attr.items():
            edge_att = {k: v for k, v in attributes.items() if 'transitivity' in k}
            multigraph.add_edge(n, v, origin="graphB", **edge_att)

        for (n, v), attributes in gA_multigraph_edges_2_attr.items():
            edge_att = {k: v for k, v in attributes.items() if 'transitivity' in k}
            multigraph.add_edge(n, v, origin="graphA", **edge_att)

        rename_cgc_graph = [f'{i}-' for i in range(len(cgc_graphs))]
        all_gc_graph = nx.union_all([gcA_graph, gcB_graph] + cgc_graphs, ['', ''] + rename_cgc_graph)
        multigraph_and_cgc_graphs = nx.union(nx.MultiGraph(all_gc_graph), multigraph, rename=['', 'multi'])

    return cgc_infos, multigraph_and_cgc_graphs


def pass_graph_attribute_to_multigraph(meta_nodes_2_attributes, gcA_graph, node_mapper):
    graph_node_2_attributes = {}
    for n, attribute in gcA_graph.nodes(data=True):
        data = {f"{node_mapper}_{k}": v for k, v in attribute.items()}
        graph_node_2_attributes[n] = data

    for meta_node, attributes in meta_nodes_2_attributes.items():
        graph_node = attributes[node_mapper]

        attributes.update(graph_node_2_attributes[graph_node])


def compare_pair_of_context_graphs(context_pair: Tuple[GeneContext, GeneContext],
                                   gene_fam_2_cluster_fam: Dict[str, str],
                                   return_multigraph: bool, outdir):
    contextA, contextB = sorted(context_pair)
    # Get conserved genomic context from the two context graph 
    cgc_infos, multigraph = get_conserved_genomics_contexts(contextA.graph, contextB.graph, gene_fam_2_cluster_fam,
                                                            return_multigraph=return_multigraph)

    if len(cgc_infos) > 1:
        print(f"MORE THAN ONE CONSERVED GENOMIC CONTEXTS FROM TWO CONTEXTS GRAPHS: {contextA.ID}_vs_{contextB.ID}")

    cgc_sizes = [info['mean_size'] for info in cgc_infos]
    cgc_scores = [info['score'] for info in cgc_infos]

    # Compute simple metrics 
    contextA_clst_family = {gf.akin for gf in contextA.families}
    contextB_clst_family = {gf.akin for gf in contextB.families}
    shared_family = len(contextA_clst_family & contextB_clst_family)
    clst_family_jaccard = shared_family / len(contextA_clst_family | contextB_clst_family)

    # Writting multigraph
    if return_multigraph and len(cgc_infos) > 0:
        logging.debug(f"Writting multigraph{(contextA, contextB)}, {multigraph}")

        multigraph_outdir = outdir / 'context_pair_graph_comparison'  # / f"{contextA.ID}_vs_{contextB.ID}"

        mkdir(multigraph_outdir, True)

        # multigraph, graphA, graphB = graphs
        nx.readwrite.graphml.write_graphml(multigraph, multigraph_outdir / f"{contextA.ID}_vs_{contextB.ID}.graphml")
        # nx.readwrite.graphml.write_graphml(graphA, multigraph_outdir / "graph_A.graphml")
        # nx.readwrite.graphml.write_graphml(graphB, multigraph_outdir / "graph_B.graphml")

    # return comparison if cgc are found..
    if len(cgc_infos) > 0:
        return contextA.ID, contextB.ID, {'family_compo_jaccard': clst_family_jaccard,
                                          "shared family": shared_family,
                                          "max_cgc_score": round(max(cgc_scores), 2),
                                          "cgc_scores": '|'.join([str(round(score, 2)) for score in cgc_scores]),
                                          "cgc_count": len(cgc_infos),
                                          "cgc_sizes": ':'.join([str(size) for size in cgc_sizes]), }
    else:
        return contextA.ID, contextB.ID, None


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


def compare_gene_contexts_on_cluster_families(gene_contexts: List[GeneContext], min_jaccard, max_workers: int,
                                              disable_bar: bool) -> List[GeneContext]:
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

    context_pairs = combinations(gene_contexts, 2)
    comparison_arguments = ((p, min_jaccard) for p in context_pairs)
    pair_count = (len(gene_contexts) ** 2 - len(gene_contexts)) / 2

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for gcA, gcB, metrics in tqdm(executor.map(launch_context_comparison, comparison_arguments, chunksize=5),
                                      total=pair_count,
                                      disable=disable_bar, unit="context pair"):
            if metrics:
                context_graph.add_edge(gcA, gcB, **metrics)

    logging.info(f'Context graph: {context_graph}')
    return context_graph


def compare_gene_contexts_graph_mp(gene_contexts: List[GeneContext],
                                   gene_fam_2_cluster_fam: Dict[str, str],
                                   max_workers: int, disable_bar: bool,
                                   return_multigraph: bool, outdir) -> List[GeneContext]:
    """
    Compares gene contexts by looking at their context graphs.

    :param gene_contexts: A list of GeneContext objects to be compared.
    :param gene_fam_2_cluster_fam: dict mapping gene family name to cluster family name
    :param max_workers: The maximum number of worker processes for parallel execution.
    :param disable_bar: A boolean flag indicating whether to disable the progress bar.
    :return: A list of GeneContext objects.
    """
    context_graph = nx.Graph()
    for gc in gene_contexts:
        context_graph.add_node(gc.ID, pangenome=gc.pangenome)

    context_pairs = combinations(gene_contexts, 2)
    comparison_arguments = ((p, gene_fam_2_cluster_fam, return_multigraph, outdir) for p in context_pairs)
    pair_count = (len(gene_contexts) ** 2 - len(gene_contexts)) / 2

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for gcA, gcB, metrics in tqdm(
                executor.map(launch_compare_pair_of_context_graphs, comparison_arguments, chunksize=5),
                total=pair_count,
                disable=disable_bar, unit="context pair"):
            if metrics:
                context_graph.add_edge(gcA, gcB, **metrics)

    logging.info(f'Context graph: {context_graph}')
    return context_graph


def context_comparison(pangenomes: Pangenomes, family_clusters: bool,
                       synteny_score: float, lock: Lock, output: Path, tmpdir: Path, task: int = 1,
                       threads_per_task: int = 1,
                       disable_bar: bool = False, force: bool = False, ppanggolin_context_args: dict = None):
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
    :param ppanggolin_context_args: Additional keyword arguments to search context using ppanggolin.
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
        # with Manager() as manager:
        #     lock = manager.Lock()

        #     manager_dict = manager.dict()

        check_run_context_arguments(ppanggolin_context_args)

        gene_contexts = launch_ppanggolin_context(pangenome_to_path, ppanggolin_context_args, output=output,
                                                  tmpdir=tmpdir,
                                                  max_workers=task, threads_per_task=threads_per_task,
                                                  disable_bar=disable_bar)

    # write gene context summary
    summary_out_table = output / "gene_context_summary.tsv"
    logging.info(f'Writting gene context summary: {summary_out_table} ')
    write_context_summary(gene_contexts, summary_out_table)

    if family_clusters:
        logging.info(f"Retrieving family clusters from existing clustering: {family_clusters}")

        # Parse the given cluster family results
        gene_family_to_family_cluster, cluster2family = read_clustering_table(family_clusters)
        # remove fragmentation info
        gene_family_to_family_cluster = {gf: fc for gf, (fc, is_fragmented) in gene_family_to_family_cluster.items()}

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

    # context_graph_clstr_families = compare_gene_contexts_on_cluster_families(gene_contexts, min_jaccard, task, disable_bar)

    # Compare gene contexts based on their synteny information
    context_synteny_graph = compare_gene_contexts_graph_mp(gene_contexts, gene_family_to_family_cluster, task,
                                                           disable_bar,
                                                           return_multigraph=True, outdir=output)

    # context_graph_merged = nx.compose(context_graph_clstr_families, context_synteny_graph)

    # add node attributes
    nx.set_node_attributes(context_synteny_graph, {gc.ID: gc.summarize() for gc in gene_contexts})

    context_graph_file = output / f"context.graphml"
    logging.info(f'Writting gene context graph: {context_graph_file}')
    nx.readwrite.graphml.write_graphml(context_synteny_graph, context_graph_file)


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    check_context_comparison(args)
    pangenomes, tmpdir = common_launch(args, check_pangenome_for_context_search,
                                       {"need_families": True, 'need_annotations': True})

    #
    #
    #     context_comparison(pan_to_path, args.context_results, args.family_clusters,
    #                        ppanggolin_context_args=run_context_args,
    #                        synteny_score=args.synteny_score, lock=lock, output=args.output,
    #                        tmpdir=args.tmpdir, task=args.task, threads_per_task=args.threads_per_task,
    #                        disable_bar=args.disable_prog_bar, force=args.force)
    #
    # elif args.compare_subcommand == "module":
    #     raise NotImplementedError(f"Module comparison is not yet implemented.")
    if not args.keep_tmp:
        rmtree(tmpdir, ignore_errors=True)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("compare_context",
                                   description='Comparison of modules and gene contexts among pangenomes')

    parser_comparison_context(parser)
    return parser


def parser_comparison_context(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required, compare_opt, _ = parser_comparison(parser)

    compare_opt.add_argument('--synteny_score', type=int, required=False,
                             help="minimum synteny score used to filter edges between genomic contexts.")

    ## PPANGGOLIN CONTEXT ARGUMENTS
    onereq = required.add_mutually_exclusive_group(required=True)

    onereq.add_argument('-R', '--context_results', type=Path,
                        help="Already computed contexts: Tsv file with two columns: name of pangenome and path to the corresponding context results."
                             "Results can be a table (tsv) or a graph (graphml or gexf)")

    onereq.add_argument('-S', '--sequences', type=Path,
                        help="Fasta file with the sequences of interest")

    onereq.add_argument('-F', '--family', type=Path,
                        help="List of family IDs of interest from the pan")

    ppanggolin_context = parser.add_argument_group(title="PPanGGOLiN search context arguments",
                                                   description="Following arguments are used to search context "
                                                               "with PPanGGOLiN API:")
    ppanggolin_context.add_argument('--no_defrag', required=False, action="store_true",
                                    help="DO NOT Realign gene families to link fragments with"
                                         "their non-fragmented gene family.")

    ppanggolin_context.add_argument('--identity', required=False, type=float, default=0.5,
                                    help="Minimum identity threshold")
    ppanggolin_context.add_argument('--coverage', required=False, type=float, default=0.8,
                                    help="Minimum coverage threshold")
    ppanggolin_context.add_argument("-t", "--transitive", required=False, type=int, default=4,
                                    help="Size of the transitive closure used to build the graph. This indicates "
                                         "the number of non-related genes allowed in-between two related genes. "
                                         "Increasing it will improve precision but lower sensitivity a little.")
    ppanggolin_context.add_argument("-w", "--window", required=False, type=int, default=5,
                                    help="Number of neighboring genes that are considered on each side of "
                                         "a gene of interest when searching for conserved genomic contexts.")

    ppanggolin_context.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                                    help="Minimum Jaccard similarity used to filter edges between gene families. "
                                         "Increasing it will improve precision but lower sensitivity a lot.")
    ppanggolin_context.add_argument('--graph_format', help="Format of the context graph. Can be gexf or graphml.",
                                    default='graphml', choices=['graphml', 'gexf'])
