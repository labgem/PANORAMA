#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import argparse
import logging
from pathlib import Path
import tempfile
from typing import Dict, List, Tuple
from itertools import combinations
from shutil import rmtree
from multiprocessing import Manager
import subprocess
from time import time

# installed libraries
from tqdm import tqdm
import pandas as pd
from pandas import DataFrame

# local libraries
from panorama.utils import mkdir
from panorama.pangenomes import Pangenome
from panorama.format.read_binaries import load_pangenomes
from panorama.alignment.common import createdb, write_pangenomes_families_sequences

align_format = ["query", "target", "fident", "qlen", "tlen", "alnlen", "evalue", "bits"]
align_column = ["query", "target", "identity", "qlength", "tlength", "alnlength", "e_value", "bits"]  # noqa: E501


def check_align_parameters(args: argparse.Namespace):
    """
    Checks if given command argument are legit if so return a dictionary with information needed to load pangenomes.

    Args:
        args: Argument in the command line
    """
    if not args.tmpdir.exists() or not args.tmpdir.is_dir():
        raise NotADirectoryError("The given path for temporary directory is not found or not a directory."
                                 "Please check your path and try again.")


def check_pangenome_align(pangenome: Pangenome) -> None:
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


def write_alignment(query_db: Path, target_db: Path, aln_db: Path, outfile: Path, threads: int = 1):
    """
    Write alignment result provide by MMSeqs2

    Args:
        query_db: MMSeqs2 database for the query sequences
        target_db: MMSeqs2 database for the target sequences
        aln_db: MMSeqs2 database for the alignment results
        outfile: Path to the output file
        threads: Number of available threads
    """
    cmd = ["mmseqs", "convertalis", query_db.as_posix(), target_db.as_posix(), aln_db.as_posix(),
           outfile.as_posix(), "--format-output", ",".join(align_format), "--threads", str(threads)]
    logging.getLogger("PANORAMA").debug(" ".join(cmd))
    logging.getLogger("PANORAMA").info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)


def align_db(query_db: Path, target_db: Path, aln_db: Path = None, tmpdir: Path = Path(tempfile.gettempdir()),
             identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0, threads: int = 1,
             keep_tmp: bool = False) -> Path:
    """
    Align with MMSeqs2 query and target sequences

    Args:
        query_db: MMSeqs2 database for the query sequences
        target_db: MMSeqs2 database for the target sequences
        aln_db: MMSeqs2 alignment database results
        identity: Set the identity use to construct clustering [0-1]. (Defaults 0.8).
        coverage: Coverage used to construct clustering [0-1]. (Defaults 0.8).
        cov_mode: Coverage mode used by MMSeqs2 to cluster. (Defaults 0).
        tmpdir: Temporary directory for MMSeqs2. (Defaults tmp user directory).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).
        threads: Number of available threads. (Defaults 1).
    """
    aln_db = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=keep_tmp).name) if aln_db is None else aln_db
    cmd = ["mmseqs", "search", query_db.absolute().as_posix(), target_db.absolute().as_posix(),
           aln_db.absolute().as_posix(), tmpdir.name, "-a", "--min-seq-id", str(identity), "-c", str(coverage),
           "--cov-mode", str(cov_mode), "--threads", str(threads)]
    logging.getLogger("PANORAMA").debug(" ".join(cmd))
    begin_time = time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    align_time = time() - begin_time
    logging.getLogger("PANORAMA").debug(f"Aligning done in {round(align_time, 2)} seconds")
    return aln_db


def align_pangenomes_pair(pangenomes_pair: Tuple[str, str], db_pair: Tuple[Path, Path], identity: float = 0.8,
                          coverage: float = 0.8, cov_mode: int = 0, tmpdir: Path = Path(tempfile.gettempdir()),
                          keep_tmp: bool = False, threads: int = 1) -> Path:
    """
    Align with MMSeqs2 gene families from 2 different pangenomes

    Args:
        pangenomes_pair: Pangenomes pair to align pangenomes in multiprocessing
        db_pair: MMSeqs2 database for the pangenomes pair
        identity: Set the identity use to construct clustering [0-1]. (Defaults 0.8).
        coverage: Coverage used to construct clustering [0-1]. (Defaults 0.8).
        cov_mode: Coverage mode used by MMSeqs2 to cluster. (Defaults 0).
        tmpdir: Temporary directory for MMSeqs2. (Defaults tmp user directory).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).
        threads: Number of available threads. (Defaults 1).

    Returns:
        Path to the alignment results
    """
    logging.getLogger("PANORAMA").debug(f"Aligning gene families between {pangenomes_pair[0]} and {pangenomes_pair[1]}")
    aln_db = align_db(query_db=db_pair[0], target_db=db_pair[1], tmpdir=tmpdir, identity=identity,
                      coverage=coverage, cov_mode=cov_mode, threads=threads)
    logging.getLogger("PANORAMA").debug(f"Write alignment results between {pangenomes_pair[0]} "
                                        f"and {pangenomes_pair[1]}")
    aln_res = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, suffix=".tsv", delete=keep_tmp).name)
    write_alignment(db_pair[0], db_pair[1], aln_db, aln_res, threads)
    logging.getLogger("PANORAMA").debug("Write alignment done")
    return aln_res


def align_pangenomes(pangenome2db: Dict[str, Path], identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0,
                     tmpdir: Path = Path(tempfile.gettempdir()), keep_tmp: bool = False,
                     threads: int = 1, disable_bar: bool = False) -> List[Path]:
    """
    Align all gene families between pangenomes

    Args:
        pangenome2db: Dictionary with for each pangenome (key) a MMSeqs2 database of his gene families
        identity: Set the identity use to construct clustering [0-1]. (Defaults 0.8).
        coverage: Coverage used to construct clustering [0-1]. (Defaults 0.8).
        cov_mode: Coverage mode used by MMSeqs2 to cluster. (Defaults 0).
        tmpdir: Temporary directory for MMSeqs2. (Defaults tmp user directory).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).
        threads: Number of available threads. (Defaults 1).
        disable_bar: Disable progressive bar. (Defaults False).

    Returns:
        List of alignment results between pangenomes
    """
    pangenomes_pairs = list(combinations(pangenome2db.keys(), 2))
    with tqdm(total=len(pangenomes_pairs), unit='pangenomes pair', disable=disable_bar) as progress:
        logging.getLogger("PANORAMA").info("Aligning gene families between pangenomes...")
        results = []
        for pangenomes_pair in pangenomes_pairs:
            db_pair = (pangenome2db[pangenomes_pair[0]], pangenome2db[pangenomes_pair[1]])
            res = align_pangenomes_pair(pangenomes_pair, db_pair, identity, coverage, cov_mode,
                                        tmpdir, keep_tmp, threads)
            results.append(res)
            progress.update()
    return results


def merge_aln_res(align_results: List[Path], outfile: Path):
    """
    Merge pangenome pair alignments in one file

    Args:
        align_results: list of pair alignments files
        outfile: Path to the final output file
    """
    merge_res = pd.read_csv(align_results[0], sep="\t", names=align_column)
    for aln_res in align_results[1:]:
        merge_res = pd.concat([merge_res, pd.read_csv(aln_res, sep="\t", names=align_column)],
                              ignore_index=True, copy=False)
    merge_res.to_csv(outfile, sep="\t", header=True, index=False)
    logging.getLogger("PANORAMA").info(f"Pangenomes gene families similarities are saved here:"
                                       f"{outfile.absolute().as_posix()}")


def inter_pangenome_align(pangenome2families_seq: Dict[str, Path], output: Path,
                          identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0,
                          tmpdir: Path = Path(tempfile.gettempdir()), keep_tmp: bool = False, threads: int = 1,
                          disable_bar: bool = False):
    """
    Main function to align gene families between pangenomes without inside pangenome alignment

    Args:
        pangenome2families_seq: Dictionary mapping pangenome names to their respective sequence files
        output: Path to the output directory with the alignment results
        identity: Set the identity use to construct clustering [0-1]. (Defaults 0.8).
        coverage: Coverage used to construct clustering [0-1]. (Defaults 0.8).
        cov_mode: Coverage mode used by MMSeqs2 to cluster. (Defaults 0).
        tmpdir: Temporary directory for MMSeqs2. (Defaults tmp user directory).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).
        threads: Number of available threads. (Defaults 1).
        disable_bar: Disable progressive bar. (Defaults False).
    """
    pangenome2db = {}
    for name, sequences in pangenome2families_seq.items():
        logging.getLogger("PANORAMA").debug(f"Constructing sequences database for pangenome {name}")
        pangenome2db[name] = createdb([sequences], tmpdir, keep_tmp=keep_tmp)
    align_results = align_pangenomes(pangenome2db, tmpdir=tmpdir, identity=identity, coverage=coverage,
                                     cov_mode=cov_mode, threads=threads, disable_bar=disable_bar)
    logging.getLogger("PANORAMA").info("Merging pangenomes gene families alignment...")
    outfile = output / "pangenome_gf_similarities.tsv"
    merge_aln_res(align_results, outfile)


def all_against_all(families_seq: List[Path], output: Path, identity: float = 0.8, coverage: float = 0.8,
                    cov_mode: int = 0, tmpdir: Path = Path(tempfile.gettempdir()), keep_tmp: bool = False,
                    threads: int = 1) -> DataFrame:
    """
    Main function to align all gene families from all pangenomes with inside alignment

    Args:
        families_seq: List of path to gene families sequences
        output: Path to the output directory with the alignment results
        tmpdir: Temporary directory for MMSeqs2
        identity: Set the identity use to construct clustering [0-1]. (Defaults 0.9)
        coverage: Coverage used to construct clustering [0-1]. (Defaults 0.8)
        cov_mode: Coverage mode used by MMSeqs2 to cluster. (Defaults 1).
        threads: Number of available threads. (Defaults 1).
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).

    Returns:
        Dataframe with alignment results
    """
    merge_db = createdb(families_seq, tmpdir)
    aln_db = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=keep_tmp).name)
    logging.getLogger("PANORAMA").debug("Aligning all gene families...")
    align_db(query_db=merge_db, target_db=merge_db, aln_db=aln_db, tmpdir=tmpdir,
             identity=identity, coverage=coverage, cov_mode=cov_mode, threads=threads)
    aln_res = Path(tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, suffix=".tsv", delete=keep_tmp).name)
    write_alignment(query_db=merge_db, target_db=merge_db, aln_db=aln_db, outfile=aln_res, threads=threads)
    logging.getLogger("PANORAMA").debug("Write alignment done")
    align_df = pd.read_csv(aln_res, sep="\t", names=align_column)
    outfile = output / "pangenome_gf_similarities.tsv"
    align_df.to_csv(outfile, sep="\t", header=True, index=False)
    logging.getLogger("PANORAMA").info(f"Pangenomes gene families similarities are saved here: "
                                       f"{outfile.absolute().as_posix()}")
    return align_df


def launch(args):
    """
    Launch functions to align gene families from pangenomes

    Args:
        args: argument given in CLI
    """
    check_align_parameters(args)
    mkdir(args.output, args.force)
    manager = Manager()
    lock = manager.Lock()
    need_info = {"need_families": True, 'need_families_sequences': True}
    pangenomes = load_pangenomes(pangenome_list=args.pangenomes, need_info=need_info,
                                 check_function=check_pangenome_align, max_workers=args.threads,
                                 lock=lock, disable_bar=args.disable_prog_bar)
    tmpdir = Path(tempfile.mkdtemp(dir=args.tmpdir))
    pangenome2families_seq = write_pangenomes_families_sequences(pangenomes=pangenomes, tmpdir=tmpdir, lock=lock,
                                                                 threads=args.threads,
                                                                 disable_bar=args.disable_prog_bar)
    if args.inter_pangenomes:
        inter_pangenome_align(pangenome2families_seq=pangenome2families_seq, output=args.output, identity=args.identity,
                              coverage=args.coverage, cov_mode=args.cov_mode, tmpdir=tmpdir, keep_tmp=args.keep_tmp,
                              threads=args.threads, disable_bar=args.disable_prog_bar)
    elif args.all_against_all:
        all_against_all(families_seq=list(pangenome2families_seq.values()), output=args.output,
                        identity=args.identity, coverage=args.coverage, cov_mode=args.cov_mode,
                        tmpdir=tmpdir, threads=args.threads, keep_tmp=args.keep_tmp)
    else:
        raise argparse.ArgumentError(argument=args, message="You must choose between inter_pangenome alignment or "
                                                            "all_against_all alignment")
    if not args.keep_tmp:
        rmtree(tmpdir, ignore_errors=True)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    Args:
        sub_parser: sub_parser for align command

    Returns:
        argparse.ArgumentParser: parser arguments for align command
    """
    parser = sub_parser.add_parser("align", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_align(parser)
    return parser


def parser_align(parser):
    """
    Add argument to parser for align command

    Args:
        parser: parser for align argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument('-o', '--output', required=True, type=Path,
                          help="Output directory where the file(s) will be written")
    exclusive = required.add_mutually_exclusive_group(required=True)
    exclusive.add_argument("--inter_pangenomes", action="store_true",
                           help="Align only gene families between pangenomes and not inside pangenome. "
                                "Not compatible with --all_against_all option")
    exclusive.add_argument("--all_against_all", action="store_true",
                           help="Align gene families between pangenomes and intra-pangenome."
                                "Not compatible with --inter_pangenome option")
    mmseqs = parser.add_argument_group(title="MMSeqs2 arguments",
                                       description="The following arguments are optional."
                                                   "Look at MMSeqs2 documentation for more information.")
    mmseqs.add_argument('--identity', required=False, type=float, default=0.5,
                        help="min identity percentage threshold")
    mmseqs.add_argument('--coverage', required=False, type=float, default=0.8,
                        help="min coverage percentage threshold")
    mmseqs.add_argument('--cov_mode', required=False, type=int, default=0,
                        help="covery_mode")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--tmpdir", required=False, type=Path, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--keep_tmp", required=False, default=False, action="store_true",
                          help="Keeping temporary files (useful for debugging).")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads.")
