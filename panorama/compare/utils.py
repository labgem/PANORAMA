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
from panorama.alignment.cluster import cluster_gene_families, write_clustering


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
                       "kmer_per_seq": args.kmer_per_seq, "coverage": args.coverage, "identity": args.identity,
                       "cov_mode": args.cov_mode, "eval": args.eval, "max_seq_len": args.max_seq_len,
                       "max_reject": args.max_reject, "align_mode": args.align_mode, "clust_mode": args.clust_mode,
                       "reassign": args.reassign}

        clust_res = cluster_gene_families(pangenomes=pangenomes, method=args.method, mmseqs2_opt=mmseqs2_opt,
                                          tmpdir=tmpdir, keep_tmp=args.keep_tmp, threads=args.cpus, lock=lock,
                                          disable_bar=args.disable_prog_bar)

        cluster = tmpdir / "pangenome_gf_clustering.tsv"
        write_clustering(clust_res, cluster)

    else:
        cluster = args.cluster
    pangenomes.read_clustering(cluster, args.disable_prog_bar)
    return pangenomes, tmpdir


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

    mmseqs2 = parser.add_argument_group(title="MMSeqs2 arguments",
                                        description="The following arguments are optional."
                                                    "Look at MMSeqs2 documentation for more information."
                                                    "If one MMSeqs2 is missing fell free to ask on our github or "
                                                    "to make a pull request")
    mmseqs2.add_argument("-m", "--method", required=False, type=str, choices=["linclust", "cluster"],
                         default="linclust", help="Choose MMSeqs2 clustering methods:"
                                                  "\t-linclust fast but less sensitive clustering"
                                                  "\t-cluster slower but more sensitive clustering")
    mmseqs2.add_argument("--sensitivity", type=int, required=False, nargs='?', default=4,
                         help='sensitivity used with MMSeqs2')
    mmseqs2.add_argument("--max_seqs", required=False, type=int, default=400, nargs='?',
                         help="Maximum results per query sequence allowed to pass the prefilter")
    mmseqs2.add_argument("--min_ungapped", required=False, type=int, default=1, nargs='?',
                         help='Accept only matches with ungapped alignment score above threshold')
    mmseqs2.add_argument("--comp_bias_corr", required=False, type=float, nargs='?', default=1,
                         help='Correct for locally biased amino acid composition')
    mmseqs2.add_argument("--kmer_per_seq", required=False, nargs='?', type=int, default=80,
                         help='k-mers per sequence')
    mmseqs2.add_argument("--cov_mode", required=False, nargs="?", type=int, default=0,
                         help="Coverage mode used by MMSeqs2 to cluster")
    mmseqs2.add_argument("--align_mode", required=False, nargs="?", type=int, default=2,
                         help="Alignment mode used by MMSeqs2 to cluster")
    mmseqs2.add_argument("--eval", required=False, nargs="?", default=0.001, type=float,
                         help="List matches below this E-value")
    mmseqs2.add_argument("--max_seq_len", required=False, nargs="?", default=32768, type=int,
                         help="Maximum sequence length")
    mmseqs2.add_argument("--max_reject", required=False, nargs="?", default=2147483647, type=int,
                         help="Maximum rejected alignments before alignment calculation for a query is stopped")
    mmseqs2.add_argument("--clust_mode", required=False, nargs="?", default=1,
                         help="Clustering mode used by MMSeqs2 to cluster")
    mmseqs2.add_argument("--reassign", required=False, action="store_false", default=True,
                         help="Correct errors from cascaded clustering")

    optional = parser.add_argument_group(title="Optional arguments")

    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                          help="Keeping temporary files (useful for debugging).")
    optional.add_argument("-c", "--cpus", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    return required, compare_opt, optional
