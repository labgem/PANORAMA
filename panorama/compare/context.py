#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from typing import Dict, Tuple
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Lock
import logging
from typing import Dict, Union, List, Set
from itertools import combinations
import networkx as nx


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

        context_graph = nx.subgraph_view(contexts_graph, filter_node=lambda n: n in families_in_context)

        
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


def compare_pair_of_contexts(context_pair: Tuple[GeneContext, GeneContext], min_jaccard) -> Tuple[GeneContext, GeneContext, float]:
    """
    Compares a pair of gene contexts and calculates the Jaccard similarity between their family clusters.

    :param context_pair: A tuple containing two GeneContext objects to be compared.
    :return: A tuple containing the two GeneContext objects and the Jaccard similarity between their family clusters.
    """

    contextA, contextB = context_pair
    contextA_clst_family =  {gf.family_cluster for gf in contextA.families}
    contextB_clst_family = {gf.family_cluster for gf in contextB.families}
    shared_family = len(contextA_clst_family & contextB_clst_family)

    clst_family_jaccard = shared_family / len(contextA_clst_family | contextB_clst_family)
    if clst_family_jaccard >= min_jaccard: 
        return contextA.ID, contextB.ID, {'jaccard':clst_family_jaccard, "shared family":shared_family}
    else: 
        return contextA.ID, contextB.ID, None
    
def launch_context_comparison(pack: tuple) -> tuple:
    """ 
    Allow to launch in multiprocessing the context comparison

    :param pack: Pack of argument for context comparison

    :return: edge metrics 
    """
    return compare_pair_of_contexts(*pack)


def compare_gene_contexts(gene_contexts: List[GeneContext], min_jaccard, max_workers: int, disable_bar: bool) -> List[GeneContext]:
    """
    Compares gene contexts by calculating the Jaccard similarity between their family clusters.

    :param gene_contexts: A list of GeneContext objects to be compared.
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

def compare_gene_contexts_on_synteny(gene_contexts: List[GeneContext], min_score, max_workers: int, disable_bar: bool) -> List[GeneContext]:
    """
    Compares gene contexts by calculating the Jaccard similarity between their family clusters.

    :param gene_contexts: A list of GeneContext objects to be compared.
    :param max_workers: The maximum number of worker processes for parallel execution.
    :param disable_bar: A boolean flag indicating whether to disable the progress bar.
    :return: A list of GeneContext objects.
    """
    context_graph = nx.Graph()
    for gc in gene_contexts:
        context_graph.add_node(gc.ID, pangenome=gc.pangenome)

    context_pairs =  combinations(gene_contexts, 2)
    comparison_arguments = ((p, min_score) for p in context_pairs)
    pair_count = (len(gene_contexts) **2 - len(gene_contexts)) / 2

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for gcA, gcB, metrics in tqdm(executor.map(launch_context_comparison, comparison_arguments, chunksize=5), total=pair_count, 
                                      disable=disable_bar, unit="context pair"):
            if metrics:
                context_graph.add_edge( gcA, gcB, **metrics)

    logging.info(f'Context graph: {context_graph}')
    return context_graph


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
        
    else:
        # run cluster familly
        # family_clusters_file = panorama cluster function
        raise NotImplementedError

    # TODO Check that all gene family in context have a family cluster 

    # add family cluster info in gene contexts 
    for gene_context in gene_contexts:
        for gene_family in gene_context.families:
            family_cluster, is_fragmented = gene_family_to_family_cluster[gene_family.name]
            gene_family.add_family_cluster(family_cluster)

    # Compare gene contexts based on their family clusters  

    context_graph = compare_gene_contexts(gene_contexts, min_jaccard, task, disable_bar)

    context_synteny_graph = compare_gene_contexts_on_synteny(gene_contexts, min_jaccard, task, disable_bar)
    
    context_graph_file = output / f"context.graphml"
    logging.info(f'Writting gene context graph: {context_graph_file}')
    nx.readwrite.graphml.write_graphml(context_graph, context_graph_file)

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
