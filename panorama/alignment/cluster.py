#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from pathlib import Path
import tempfile
from typing import Dict, Union
from multiprocessing import Manager, Lock
import subprocess
from time import time

# installed libraries
import pandas as pd

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenomes
from panorama.format.read_binaries import load_multiple_pangenomes
from panorama.alignment.common import get_gf_pangenomes, createdb

clust_col_names = ["cluster_id", "referent", "in_clust"]


def write_clustering(clust_res: Path, outfile: Path):
    """Write the clustering results with clustering ID

    :param clust_res: Path to the clustering results
    :param outfile: Path to the final output file
    """

    logging.debug("Begin writing clustering...")
    clust_df = pd.read_csv(clust_res, sep="\t", names=clust_col_names[1:])
    clust_id = {index: ref_fam for index, ref_fam in enumerate(clust_df[clust_col_names[1]].unique().tolist())}
    clust_id_df = pd.DataFrame.from_dict(clust_id, orient="index").reset_index()
    clust_id_df.columns = clust_col_names[0:2]
    merge_clust = clust_id_df.merge(clust_df, on=clust_col_names[1])
    merge_clust.to_csv(outfile, sep="\t", header=True, index=False)
    logging.info(f"Pangenomes gene families similarities are saved here: {outfile.absolute().as_posix()}")


def create_tsv(db: Path, clust: Path, output: Path, threads: int = 1):
    """Create a TSV from clustering result thanks to MMSeqs2

    :param db: database of the sequences
    :param clust: clustering results
    :param output: output directory
    :param threads: number of cpu available
    """

    cmd = ["mmseqs", "createtsv", db.absolute().as_posix(), db.absolute().as_posix(), clust.absolute().as_posix(),
           output.absolute().as_posix(), "--threads", str(threads), "--full-header"]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)


def linclust_launcher(seq_db: Path, mmseqs2_opt: Dict[str, Union[int, float, str]], tmpdir: tempfile.TemporaryDirectory,
                      lclust_db: Path = None, threads: int = 1) -> Path:
    """Launch a linclust clustering provide by MMSeqs2 on pangenomes gene families

    :param seq_db: Path to the MMSeqs2 database with all gene families
    :param mmseqs2_opt: Dictionnary with all the option provide by MMSeqs2
    :param tmpdir: Temporary directory for MMSeqs2
    :param lclust_db: Path to resulting database if you prefer a specific file
    :param threads: Number of available threads

    :return: Path to resulting database
    """

    if lclust_db is None:
        lclust_db = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, delete=False).name)
    logging.debug("Clustering all gene families with linclust process...")
    cmd = list(map(str, ['mmseqs', 'cluster', seq_db.absolute().as_posix(), lclust_db.absolute().as_posix(),
                         tmpdir.name, '--threads', threads, '--comp-bias-corr', mmseqs2_opt["comp_bias_corr"],
                         "--kmer-per-seq", mmseqs2_opt["kmer_per_seq"], "--min-seq-id", mmseqs2_opt["identity"],
                         "-c", mmseqs2_opt["coverage"], "--cov-mode", mmseqs2_opt["cov_mode"], "-e",
                         mmseqs2_opt["eval"],
                         "--alignment-mode", mmseqs2_opt["align_mode"], "--max-seq-len", mmseqs2_opt["max_seq_len"],
                         "--max-rejected", mmseqs2_opt["max_reject"], "--cluster-mode", mmseqs2_opt["clust_mode"]]
                   )
               )
    logging.debug(" ".join(cmd))
    begin_time = time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    lclust_time = time() - begin_time
    logging.debug(f"Linclust done in {round(lclust_time, 2)} seconds")
    return lclust_db


def cluster_launcher(seq_db: Path, mmseqs2_opt: Dict[str, Union[int, float, str]], tmpdir: tempfile.TemporaryDirectory,
                     cluster_db: Path = None, threads: int = 1) -> Path:
    """Launch a cluster clustering provide by MMSeqs2 on pangenomes gene families

    :param seq_db: Path to the MMSeqs2 database with all gene families
    :param mmseqs2_opt: Dictionnary with all the option provide by MMSeqs2
    :param tmpdir: Temporary directory for MMSeqs2
    :param cluster_db: Path to resulting database if you prefer a specific file
    :param threads: Number of available threads

    :return: Path to resulting database
    """

    if cluster_db is None:
        cluster_db = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, delete=False).name)
    logging.debug("Clustering all gene families with cluster process...")
    cmd = list(map(str,
                   ['mmseqs', 'cluster', seq_db.absolute().as_posix(), cluster_db.absolute().as_posix(),
                    tmpdir.name, '--threads', threads, '--max-seqs', mmseqs2_opt["max_seqs"], '--min-ungapped-score',
                    mmseqs2_opt["min_ungapped"], '--comp-bias-corr', mmseqs2_opt["comp_bias_corr"], '-s',
                    mmseqs2_opt["sensitivity"], "--kmer-per-seq", mmseqs2_opt["kmer_per_seq"], "--min-seq-id",
                    mmseqs2_opt["identity"], "-c", mmseqs2_opt["coverage"], "--cov-mode", mmseqs2_opt["cov_mode"],
                    "-e", mmseqs2_opt["eval"], "--alignment-mode", mmseqs2_opt["align_mode"], "--max-seq-len",
                    mmseqs2_opt["max_seq_len"], "--max-rejected", mmseqs2_opt["max_reject"], "--cluster-mode",
                    mmseqs2_opt["clust_mode"]]
                   )
               )
    logging.debug(" ".join(cmd))
    begin_time = time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    lclust_time = time() - begin_time
    logging.debug(f"Linclust done in {round(lclust_time, 2)} seconds")
    return cluster_db


def cluster_gene_families(pangenomes: Pangenomes, method: str, mmseqs2_opt: Dict[str, Union[int, float, str]],
                          lock: Lock, tmpdir: tempfile.TemporaryDirectory, threads: int = 1,
                          disable_bar: bool = False) -> Path:
    """Cluster pangenomes gene families with MMSeqs2

    :param pangenomes: Pangenomes obejct containing all the pangenome to cluster
    :param method: Choosen method to cluster pangenomes gene families
    :param mmseqs2_opt: Dictionnary with all the option provide by MMSeqs2
    :param lock: Global lock for multiprocessing execution
    :param tmpdir: Temporary directory for MMSeqs2
    :param threads: Number of available threads
    :param disable_bar: Disable progressive bar

    :return: Clustering results from MMSeqs2
    """

    gf_seqs = get_gf_pangenomes(pangenomes=pangenomes, create_db=False, lock=lock, tmpdir=tmpdir, threads=threads,
                                disable_bar=disable_bar)
    merge_db = createdb(list(gf_seqs.values()), tmpdir)
    if method == "linclust":
        clust_db = linclust_launcher(seq_db=merge_db, mmseqs2_opt=mmseqs2_opt, tmpdir=tmpdir, threads=threads)
    elif method == "cluster":
        clust_db = cluster_launcher(seq_db=merge_db, mmseqs2_opt=mmseqs2_opt, tmpdir=tmpdir, threads=threads)
    else:
        raise ValueError("You must choose between linclut or cluster for clustering")
    clust_res = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, suffix=".tsv", delete=False).name)
    logging.getLogger().info("Writting clustering in tsv file")
    create_tsv(db=merge_db, clust=clust_db, output=clust_res, threads=threads)
    logging.debug("Clustering done")
    return clust_res


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """

    # check_parameter(args)
    mkdir(args.output, args.force)
    mmseqs2_opt = {"max_seqs": args.max_seqs, "min_ungapped": args.min_ungapped, "comp_bias_corr": args.comp_bias_corr,
                   "sensitivity": args.sensitivity, "kmer_per_seq": args.kmer_per_seq, "coverage": args.coverage,
                   "identity": args.identity, "cov_mode": args.cov_mode, "eval": args.eval,
                   "max_seq_len": args.max_seq_len, "max_reject": args.max_reject, "align_mode": args.align_mode,
                   "clust_mode": args.clust_mode, "reassign": args.reassign}
    manager = Manager()
    lock = manager.Lock()
    pangenomes = load_multiple_pangenomes(pangenome_list=args.pangenomes, need_info={"need_families": True}, lock=lock,
                                          max_workers=args.threads, disable_bar=args.disable_prog_bar)
    tmpdir = tempfile.TemporaryDirectory(dir=args.tmpdir)
    clust_res = cluster_gene_families(pangenomes=pangenomes, method=args.method, mmseqs2_opt=mmseqs2_opt, lock=lock,
                                      tmpdir=tmpdir, threads=args.threads, disable_bar=args.disable_prog_bar)

    outfile = args.output / "pangenome_gf_clustering.tsv"
    write_clustering(clust_res, outfile)
    tmpdir.cleanup()


def subparser(sub_parser) -> argparse.ArgumentParser:
    """Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return: parser arguments for align command
    """

    parser = sub_parser.add_parser("cluster", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_clust(parser)
    return parser


def parser_clust(parser):
    """Parser for specific argument of annot command

    :param parser: parser for annot argument
    """

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")
    required.add_argument("-m", "--method", required=True, type=str, choices=["linclust", "cluster"],
                          help="Choose MMSeqs2 clustering methods:"
                               "\t-linclust fast but less sensitive clustering"
                               "\t-cluster slower but more sensitive clustering")
    mmseqs2 = parser.add_argument_group(title="MMSeqs2 arguments",
                                        description="The following arguments are optional."
                                                    "Look at MMSeqs2 documentation for more information."
                                                    "If one MMSesq2 is missing fell free to ask on our github or "
                                                    "to make a pull request")
    mmseqs2.add_argument("-s", "--sensitivity", type=int, required=False, nargs='?', default=4,
                         help='sensitivity used with MMSeqs2')
    mmseqs2.add_argument("--max_seqs", required=False, type=int, default=400, nargs='?',
                         help="Maximum results per query sequence allowed to pass the prefilter")
    mmseqs2.add_argument("--min_ungapped", required=False, type=int, default=1, nargs='?',
                         help='Accept only matches with ungapped alignment score above threshold')
    mmseqs2.add_argument("--comp_bias_corr", required=False, type=float, nargs='?', default=1,
                         help='Correct for locally biased amino acid composition')
    mmseqs2.add_argument("--kmer_per_seq", required=False, nargs='?', type=int, default=80,
                         help='k-mers per sequence')
    mmseqs2.add_argument("--identity", required=False, nargs="?", default=0.8, type=float,
                         help="Set the identity use to construct clustering [0-1]")
    mmseqs2.add_argument("--coverage", required=False, nargs="?", type=float, default=0.8,
                         help="Set the coverage use to construct clustering [0-1]")
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
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
