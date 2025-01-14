#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from pathlib import Path
import tempfile
from typing import Any, List, Set, Tuple
from multiprocessing import Manager, Lock

# installed libraries
import networkx as nx

# local libraries
from panorama.pangenomes import Pangenomes
from panorama.geneFamily import GeneFamily
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.cluster import cluster_gene_families, write_clustering, parser_mmseqs2_cluster


def compute_frr(queries: Set[GeneFamily], targets: Set[GeneFamily]) -> Tuple[float, float, int]:
    akins = {query_gf.akin.ID for query_gf in queries for target_gf in targets if query_gf.akin == target_gf.akin}
    min_frr = len(akins) / min(len(queries), len(targets))
    max_frr = len(akins) / max(len(queries), len(targets))

    return min_frr, max_frr, len(akins)


def cluster_on_frr(graph: nx.Graph, frr_metrics: str) -> List[Set[Any]]:
    partitions = nx.algorithms.community.louvain_communities(graph, weight=frr_metrics)

    # Add partition index in node attributes
    for i, cluster_nodes in enumerate(partitions):
        nx.set_node_attributes(graph, {node: f"cluster_{i}" for node in cluster_nodes}, name=f"{frr_metrics}_cluster")

    logging.info(f"Graph has {len(partitions)} clusters using {frr_metrics}")
    return partitions


def common_launch(args, check_func, need_info: dict, **kwargs) -> Tuple[Pangenomes, Path, Manager, Lock]:
    manager = Manager()
    lock = manager.Lock()
    if args.cluster is None:
        need_info["need_families_sequences"] = True
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_func, max_workers=args.cpus, lock=lock,
                                 disable_bar=args.disable_prog_bar, **kwargs)
    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    if args.cluster is None:
        mmseqs2_opt = {"max_seqs": args.max_seqs, "min_ungapped": args.min_ungapped,
                       "comp_bias_corr": args.comp_bias_corr, "sensitivity": args.sensitivity,
                       "kmer_per_seq": args.kmer_per_seq, "identity": args.clust_identity,
                       "coverage": args.clust_coverage, "cov_mode": args.clust_cov_mode, "eval": args.eval,
                       "max_seq_len": args.max_seq_len, "max_reject": args.max_reject, "align_mode": args.align_mode,
                       "clust_mode": args.clust_mode, "reassign": args.reassign}

        clust_res = cluster_gene_families(pangenomes=pangenomes, method=args.method, mmseqs2_opt=mmseqs2_opt,
                                          tmpdir=tmpdir, keep_tmp=args.keep_tmp, threads=args.cpus, lock=lock,
                                          disable_bar=args.disable_prog_bar)

        cluster = tmpdir / "pangenome_gf_clustering.tsv"
        write_clustering(clust_res, cluster)

    else:
        cluster = args.cluster
    pangenomes.read_clustering(cluster, args.disable_prog_bar)
    return pangenomes, tmpdir, manager, lock


def parser_comparison(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')

    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")

    compare_opt = parser.add_argument_group(title="Comparison optional arguments")
    compare_opt.add_argument('--cluster', required=False, type=Path, nargs='?', default=None,
                             help="A tab-separated file listing the cluster names, the family IDs")
    compare_opt.add_argument('--frr_cutoff', required=False, type=tuple, default=(0.5, 0.8), nargs=2,
                             help="The frr (Families Repertoire Relatedness) is used to assess the similarity between two "
                                  "elements based on their gene families.\n"
                                  "\tThe 'min_frr': Computes the number of gene families shared between the two elements "
                                  "and divides it by the smaller number of gene families among the two elements.\n"
                                  "\tThe 'max_frr': Computes the number of gene families shared between the two elements "
                                  "and divides it by the larger number of gene families among the two elements."
                             )
    cluster = parser_mmseqs2_cluster(parser)
    cluster.description = "MMSeqs2 arguments to cluster gene families if no cluster result given"
    cluster.add_argument("--method", required=False, type=str, choices=["linclust", "cluster"], default="linclust",
                         help="Choose MMSeqs2 clustering methods:"
                              "\t-linclust fast but less sensitive clustering"
                              "\t-cluster slower but more sensitive clustering")

    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument('--graph_formats', required=False, type=str, choices=['gexf', "graphml"], nargs="+",
                          default=None, help="Format of the output graph.")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                          help="Keeping temporary files (useful for debugging).")
    optional.add_argument("-c", "--cpus", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    return required, compare_opt, optional
