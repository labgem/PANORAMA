#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from pathlib import Path
import tempfile
from multiprocessing import Manager

from sqlalchemy.testing.plugin.plugin_base import logging

# installed libraries

# local libraries
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.cluster import cluster_gene_families, write_clustering, parser_mmseqs2_cluster


def common_launch(args, check_func, need_info: dict):
    manager = Manager()
    lock = manager.Lock()
    if args.cluster is None:
        need_info["need_families_sequences"] = True
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_func, max_workers=args.cpus, lock=lock,
                                 disable_bar=args.disable_prog_bar)
    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    if args.cluster is None:
        mmseqs2_opt = {"max_seqs": args.max_seqs, "min_ungapped": args.min_ungapped,
                       "comp_bias_corr": args.comp_bias_corr, "sensitivity": args.sensitivity,
                       "kmer_per_seq": args.kmer_per_seq, "identity": args.clust_identity,
                       "coverage": args.clust_coverage, "cov_mode": args.cov_mode, "eval": args.eval,
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

    parser_mmseqs2_cluster(parser)

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                          help="Keeping temporary files (useful for debugging).")
    optional.add_argument("-c", "--cpus", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    return required, compare_opt, optional
