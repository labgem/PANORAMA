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
from shutil import rmtree

# installed libraries
import pandas as pd

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenome, Pangenomes
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.common import write_pangenomes_families_sequences, createdb

clust_col_names = ["cluster_id", "referent", "in_clust"]


def check_cluster_parameters(args: argparse.Namespace):
    """
    Checks if given command argument are legit if so return a dictionary with information needed to load pangenomes.

    Args:
        args: Argument in the command line
    """
    if not args.tmpdir.exists() or not args.tmpdir.is_dir():
        raise NotADirectoryError("The given path for temporary directory is not found or not a directory."
                                 "Please check your path and try again.")


def check_pangenome_cluster(pangenome: Pangenome) -> None:
    """
    Function to check and load one pangenome in the loading process
    Args:
        pangenome: Pangenome object
    """
    if pangenome.status["genesClustered"] == "No":
        raise AttributeError("You did not cluster the genes. "
                             "See the 'ppanggolin cluster' command if you want to do that.")
    if pangenome.status["geneFamilySequences"] == "No":
        raise AttributeError("There is no sequences associated to the gene families.")


def write_clustering(clust_res: Path, outfile: Path):
    """
    Write the clustering results with clustering ID

    Args:
        clust_res: Path to the clustering results
        outfile: Path to the final output file
    """

    logging.getLogger("PANORAMA").debug("Begin writing clustering...")
    clust_df = pd.read_csv(clust_res, sep="\t", names=clust_col_names[1:])
    clust_id = {index: ref_fam for index, ref_fam in enumerate(clust_df[clust_col_names[1]].unique().tolist())}
    clust_id_df = pd.DataFrame.from_dict(clust_id, orient="index").reset_index()
    clust_id_df.columns = clust_col_names[0:2]
    merge_clust = clust_id_df.merge(clust_df, on=clust_col_names[1], how="inner", validate="one_to_many")
    merge_clust.to_csv(outfile, sep="\t", header=True, index=False)
    logging.getLogger("PANORAMA").info(f"Pangenomes gene families similarities are saved here: {outfile.as_posix()}")


def create_tsv(db: Path, clust: Path, output: Path, threads: int = 1):
    """
    Create a TSV from clustering result thanks to MMSeqs2

    Args:
        db: database of the sequences
        clust: clustering results
        output: output directory
        threads: Number of available threads. (Defaults 1).
    """

    cmd = ["mmseqs", "createtsv", db.absolute().as_posix(), db.absolute().as_posix(), clust.absolute().as_posix(),
           output.absolute().as_posix(), "--threads", str(threads), "--full-header"]
    logging.getLogger().debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)


def linclust_launcher(seq_db: Path, mmseqs2_opt: Dict[str, Union[int, float, str]], lclust_db: Path = None,
                      tmpdir: Path = None, keep_tmp: bool = False, threads: int = 1) -> Path:
    """Launch a linclust clustering provide by MMSeqs2 on pangenomes gene families

    Args:
        seq_db: List of path to gene families sequences
        mmseqs2_opt: Dictionary with all the option provide by MMSeqs2
        lclust_db: Path to resulting database if you prefer a specific file (Defaults None)
        tmpdir: Temporary directory for MMSeqs2 (Defaults user temporary directory)
        threads: Number of available threads. (Defaults 1).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).

    Returns:
        Dataframe with clustering results
    """
    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    if lclust_db is None:
        lclust_db = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=keep_tmp).name)
    logging.getLogger("PANORAMA").debug("Clustering all gene families with linclust process...")
    cmd = list(map(str, ['mmseqs', 'cluster', seq_db.absolute().as_posix(), lclust_db.absolute().as_posix(),
                         tmpdir.name, '--threads', threads, '--comp-bias-corr', mmseqs2_opt["comp_bias_corr"],
                         "--kmer-per-seq", mmseqs2_opt["kmer_per_seq"], "--min-seq-id", mmseqs2_opt["identity"],
                         "-c", mmseqs2_opt["coverage"], "--cov-mode", mmseqs2_opt["cov_mode"], "-e",
                         mmseqs2_opt["eval"],
                         "--alignment-mode", mmseqs2_opt["align_mode"], "--max-seq-len", mmseqs2_opt["max_seq_len"],
                         "--max-rejected", mmseqs2_opt["max_reject"], "--cluster-mode", mmseqs2_opt["clust_mode"]]
                   )
               )
    logging.getLogger("PANORAMA").debug(" ".join(cmd))
    begin_time = time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    lclust_time = time() - begin_time
    logging.getLogger("PANORAMA").debug(f"Linclust done in {round(lclust_time, 2)} seconds")
    return lclust_db


def cluster_launcher(seq_db: Path, mmseqs2_opt: Dict[str, Union[int, float, str]], cluster_db: Path = None,
                     tmpdir: Path = None, keep_tmp: bool = False, threads: int = 1) -> Path:
    """Launch a cluster clustering provide by MMSeqs2 on pangenomes gene families

    Args:
        seq_db: List of path to gene families sequences
        mmseqs2_opt: Dictionary with all the option provide by MMSeqs2
        cluster_db: Path to resulting database if you prefer a specific file (Defaults None)
        tmpdir: Temporary directory for MMSeqs2 (Defaults user temporary directory)
        threads: Number of available threads. (Defaults 1).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).

    Returns:
        Dataframe with clustering results
    """
    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    if cluster_db is None:
        cluster_db = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=keep_tmp).name)
    logging.getLogger("PANORAMA").debug("Clustering all gene families with cluster process...")
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
    logging.getLogger("PANORAMA").debug(" ".join(cmd))
    begin_time = time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    lclust_time = time() - begin_time
    logging.getLogger("PANORAMA").debug(f"Linclust done in {round(lclust_time, 2)} seconds")
    return cluster_db


def cluster_gene_families(pangenomes: Pangenomes, method: str, mmseqs2_opt: Dict[str, Union[int, float, str]],
                          tmpdir: Path = None, keep_tmp: bool = False, threads: int = 1, lock: Lock = None,
                          disable_bar: bool = False) -> Path:
    """Cluster pangenomes gene families with MMSeqs2

    Args:
        pangenomes: Pangenomes objects containing pangenome
        method:  Chosen method to cluster pangenomes gene families
        mmseqs2_opt: Dictionary with all the option provide by MMSeqs2
        tmpdir: Temporary directory for MMSeqs2 (Defaults user temporary directory)
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).
        threads: Number of available threads. (Defaults 1).
        lock: Lock object to write sequences in multithreading (Defaults None).
        disable_bar: Disable progressive bar (Defaults False).

    Returns:
        Dataframe with clustering results
    """
    tmpdir = Path(tempfile.gettempdir()) if tmpdir is None else tmpdir
    pangenome2families_seq = write_pangenomes_families_sequences(pangenomes=pangenomes, tmpdir=tmpdir, threads=threads,
                                                                 lock=lock, disable_bar=disable_bar)
    logging.getLogger("PANORAMA").info("Begin to cluster pangenomes gene families")
    merge_db = createdb(list(pangenome2families_seq.values()), tmpdir)
    if method == "linclust":
        clust_db = linclust_launcher(seq_db=merge_db, mmseqs2_opt=mmseqs2_opt,
                                     tmpdir=tmpdir, keep_tmp=keep_tmp, threads=threads)
    elif method == "cluster":
        clust_db = cluster_launcher(seq_db=merge_db, mmseqs2_opt=mmseqs2_opt,
                                    tmpdir=tmpdir, keep_tmp=keep_tmp, threads=threads)
    else:
        raise ValueError("You must choose between linclust or cluster methods for clustering. "
                         "Look at MMSeqs2 documentation for more information.")
    clust_res = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, suffix=".tsv", delete=keep_tmp).name)
    logging.getLogger().debug("Writing clustering in tsv file")
    create_tsv(db=merge_db, clust=clust_db, output=clust_res, threads=threads)
    logging.getLogger("PANORAMA").info("Clustering done")
    return clust_res


def launch(args):
    """
    Launch functions to align gene families from pangenomes

    Args:
        args: argument given in CLI
    """

    check_cluster_parameters(args)
    mkdir(args.output, args.force)
    mmseqs2_opt = {"max_seqs": args.max_seqs, "min_ungapped": args.min_ungapped, "comp_bias_corr": args.comp_bias_corr,
                   "sensitivity": args.sensitivity, "kmer_per_seq": args.kmer_per_seq, "coverage": args.coverage,
                   "identity": args.identity, "cov_mode": args.cov_mode, "eval": args.eval,
                   "max_seq_len": args.max_seq_len, "max_reject": args.max_reject, "align_mode": args.align_mode,
                   "clust_mode": args.clust_mode, "reassign": args.reassign}
    manager = Manager()
    lock = manager.Lock()
    need_info = {"need_families": True, 'need_families_sequences': True}
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_cluster, max_workers=args.threads,
                                 lock=lock, disable_bar=args.disable_prog_bar)

    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    clust_res = cluster_gene_families(pangenomes=pangenomes, method=args.method, mmseqs2_opt=mmseqs2_opt, tmpdir=tmpdir,
                                      keep_tmp=args.keep_tmp, threads=args.threads, lock=lock,
                                      disable_bar=args.disable_prog_bar)

    outfile = args.output / "pangenome_gf_clustering.tsv"
    write_clustering(clust_res, outfile)
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

    parser = sub_parser.add_parser("cluster")
    parser_clust(parser)
    return parser


def parser_clust(parser):
    """
    Add argument to parser for cluster command

    Args:
        parser: parser for cluster argument
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
                                                    "If one MMSeqs2 is missing fell free to ask on our github or "
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
    mmseqs2.add_argument("--identity", required=False, nargs="?", default=0.5, type=float,
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
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    optional.add_argument("--tmpdir", required=False, type=Path, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                          help="Keeping temporary files (useful for debugging).")
