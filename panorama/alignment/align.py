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
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Lock
import subprocess

# installed libraries
import tables
from tqdm import tqdm
from ppanggolin.align.alignOnPang import write_gene_fam_sequences
import pandas as pd

# local libraries
from panorama.utils import check_tsv_sanity, init_lock, mkdir
from panorama.format.read_binaries import check_pangenome_info
from panorama.pangenomes import Pangenome

align_format = ["query", "target", "fident", "qlen", "tlen", "alnlen", "evalue", "bits"]
align_column = ["query", "target", "identity", "qlength", "tlength", "alnlength", "e_value", "bits"]  # noqa: E501


def check_align(pangenomes_list: dict):
    for pangenome, info in pangenomes_list.items():
        h5f = tables.open_file(info["path"], "r")
        if not h5f.root.status._v_attrs.geneFamilySequences:
            raise Exception(f"Cannot use this function as your pangenome {pangenome} does not have gene families "
                            f"representatives associated sequences. Please use the clustering option of PPanGGOLiN.")
        h5f.close()


def createdb(seq_file: Path, tmpdir: tempfile.TemporaryDirectory, other_seq: List[Path] = None) -> Path:
    """
    Create a MMseqs2 sequence database with the given fasta file

    :param file_obj: Fasta file
    :param tmpdir: temporary directory

    :return: DB file
    """
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, delete=False)
    cmd = ["mmseqs", "createdb", seq_file.absolute().as_posix()]
    if other_seq is not None:
        cmd += list(map(Path.as_posix, map(Path.absolute, other_seq)))
    cmd += [seqdb.name, "--dbtype", "0"]
    logging.debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return Path(seqdb.name)


def get_gf_pangenome_db(pangenome_name: str, pangenome_info: dict,
                        tmpdir: tempfile.TemporaryDirectory) -> Tuple[Pangenome, Path]:
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    check_pangenome_info(pangenome, need_families=True, disable_bar=True)
    logging.debug(f"Begin create {pangenome_name} gene families database...")
    gf_fasta_path = Path(f"{tmpdir.name}/{pangenome_name}_gf.fna")
    with open(gf_fasta_path, "w") as gf_fasta:
        write_gene_fam_sequences(pangenome, gf_fasta)
        pan_db = createdb(gf_fasta_path, tmpdir)
    logging.debug(f"{pangenome_name} gene families database created")
    return pangenome, pan_db


def get_gf_pangenomes_db(pangenomes_list: dict, lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                         task: int = 1, disable_bar: bool = False) -> Dict[Pangenome, Path]:
    logging.debug("Begin create pangenomes gene families database...")
    with ProcessPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes_list), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome_name, pangenome_info in pangenomes_list.items():
                future = executor.submit(get_gf_pangenome_db, pangenome_name, pangenome_info, tmpdir)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            pangenomes_db = {}
            for future in futures:
                results = future.result()
                pangenomes_db[results[0]] = results[1]
    logging.debug("All pangenomes gene families database created")
    return pangenomes_db


def write_alignment(query_db: Path, target_db: Path, aln_db: Path, outfile: Path, threads: int = 1):
    cmd = ["mmseqs", "convertalis", query_db.as_posix(), target_db.as_posix(), aln_db.as_posix(),
           outfile.as_posix(), "--format-output", ",".join(align_format), "--threads", str(threads)]
    logging.debug(" ".join(cmd))
    logging.info("Extracting alignments...")
    subprocess.run(cmd, stdout=subprocess.DEVNULL)


def align_db(query_db: Path, target_db: Path, aln_db: Path, tmpdir: tempfile.TemporaryDirectory,
             identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0, threads: int = 1):
    cmd = ["mmseqs", "search", query_db.absolute().as_posix(), target_db.absolute().as_posix(),
           aln_db.absolute().as_posix(), tmpdir.name, "-a", "--min-seq-id", str(identity), "-c", str(coverage),
           "--cov-mode", str(cov_mode), "--threads", str(threads)]
    logging.debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    logging.debug("Aligning done")


def align_pangenomes_pair(pangenomes_pair: Tuple[Pangenome, Pangenome], db_pair: Tuple[Path, Path],
                          tmpdir: tempfile.TemporaryDirectory, identity: float = 0.8,
                          coverage: float = 0.8, cov_mode: int = 0, threads: int = 1) -> Path:
    aln_db = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, delete=False)
    logging.debug(f"Aligning gene families between {pangenomes_pair[0]} and {pangenomes_pair[1]}")
    align_db(db_pair[0], db_pair[1], Path(aln_db.name), tmpdir, identity, coverage, cov_mode, threads)
    logging.debug(f"Write alignment results between {pangenomes_pair[0]} and {pangenomes_pair[1]}")
    aln_res = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, suffix=".tsv", delete=False)
    write_alignment(db_pair[0], db_pair[1], Path(aln_db.name), Path(aln_res.name), threads)
    logging.debug("Write alignment done")
    return Path(aln_res.name)


def align_pangenomes(pangenomes_db: Dict[Pangenome, Path], lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                     identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0,
                     task: int = 1, threads_per_task: int = 1, disable_bar: bool = False) -> List[Path]:
    with ProcessPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
        pangenomes_pairs = list(combinations(pangenomes_db.keys(), 2))
        with tqdm(total=len(pangenomes_pairs), unit='pangenomes pair', disable=disable_bar) as progress:
            futures = []
            logging.info("Aligning gene families between pangenomes...")
            for pangenomes_pair in pangenomes_pairs:
                db_pair = (pangenomes_db[pangenomes_pair[0]], pangenomes_db[pangenomes_pair[1]])
                future = executor.submit(align_pangenomes_pair, pangenomes_pair, db_pair, tmpdir,
                                         identity, coverage, cov_mode, threads_per_task)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            results = []
            for future in futures:
                results.append(future.result())
    return results


def merge_aln_res(align_results: List[Path]) -> pd.DataFrame:
    merge_res = pd.read_csv(align_results[0], sep="\t", names=align_column)
    for aln_res in align_results[1:]:
        merge_res = pd.concat([merge_res, pd.read_csv(aln_res, sep="\t", names=align_column)],
                              ignore_index=True, copy=False)
    return merge_res


def inter_pangenome_align(pangenomes_list: dict, output: Path, lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                          identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0,
                          task: int = 1, threads_per_task: int = 1, disable_bar: bool = False):
    pangenomes_db = get_gf_pangenomes_db(pangenomes_list=pangenomes_list, lock=lock, tmpdir=tmpdir, task=task,
                                         disable_bar=disable_bar)
    align_results = align_pangenomes(pangenomes_db, lock, tmpdir=tmpdir,
                                     identity=identity, coverage=coverage, cov_mode=cov_mode,
                                     task=task, threads_per_task=threads_per_task, disable_bar=disable_bar)
    logging.info("Merging pangenomes gene families aligment...")
    merge_res = merge_aln_res(align_results)
    outfile = output / "pangenome_gf_similarities.tsv"
    merge_res.to_csv(outfile, sep="\t", header=True, index=False)
    logging.info(f"Pangenomes gene families similarities are saved here: {outfile.absolute().as_posix()}")


def get_gf_pangenome_sequences(pangenome_name: str, pangenome_info: dict,
                               tmpdir: tempfile.TemporaryDirectory) -> Path:
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    check_pangenome_info(pangenome, need_families=True, disable_bar=True)
    logging.debug(f"Begin create {pangenome_name} gene families database...")
    gf_seq = Path(f"{tmpdir.name}/{pangenome_name}_gf.fna")
    with open(gf_seq, "w") as gf_fasta:
        write_gene_fam_sequences(pangenome, gf_fasta)
    logging.debug(f"{pangenome_name} gene families fasta created")
    return gf_seq


def get_gf_pangenomes_sequences(pangenomes_list: dict, lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                                task: int = 1, disable_bar: bool = False) -> List[Path]:
    logging.debug("Begin create pangenomes gene families database...")
    with ProcessPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes_list), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome_name, pangenome_info in pangenomes_list.items():
                future = executor.submit(get_gf_pangenome_sequences, pangenome_name, pangenome_info, tmpdir)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            pangenomes_seq = []
            for future in futures:
                pangenomes_seq.append(future.result())
    logging.debug("All pangenomes gene families database created")
    return pangenomes_seq


def all_against_all(pangenomes_list: dict, output: Path, lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                    identity: float = 0.8, coverage: float = 0.8, cov_mode: int = 0,
                    task: int = 1, threads_per_task: int = 1, disable_bar: bool = False):
    gf_seqs = get_gf_pangenomes_sequences(pangenomes_list=pangenomes_list, lock=lock, tmpdir=tmpdir, task=task,
                                          disable_bar=disable_bar)
    merge_db = createdb(gf_seqs[0], tmpdir, gf_seqs[1:])
    aln_db = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, delete=False)
    logging.debug("Aligning all gene families...")
    align_db(merge_db, merge_db, Path(aln_db.name), tmpdir,
             identity, coverage, cov_mode, task * threads_per_task)
    aln_res = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, suffix=".tsv", delete=False)
    write_alignment(merge_db, merge_db, Path(aln_db.name),
                    Path(aln_res.name), task * threads_per_task)
    logging.debug("Write alignment done")
    align_df = pd.read_csv(Path(aln_res.name), sep="\t", names=align_column)
    outfile = output / "pangenome_gf_similarities.tsv"
    align_df.to_csv(outfile, sep="\t", header=True, index=False)


def launch(args):
    """
    Launch functions to annotate pangenomes

    :param args: Argument given
    """
    # check_parameter(args)
    mkdir(args.output, args.force)
    pan_to_path = check_tsv_sanity(args.pangenomes)
    check_align(pan_to_path)
    manager = Manager()
    lock = manager.Lock()
    tmpdir = tempfile.TemporaryDirectory(dir=args.tmpdir)
    if args.inter_pangenomes:
        inter_pangenome_align(pan_to_path, args.output, lock, tmpdir, args.identity, args.coverage, args.cov_mode,
                              args.task, args.threads_per_task, args.disable_prog_bar)
    elif args.all_against_all:
        all_against_all(pan_to_path, args.output, lock, tmpdir, args.identity, args.coverage, args.cov_mode,
                        args.task, args.threads_per_task, args.disable_prog_bar)
    else:
        raise argparse.ArgumentError(argument=args, message="You must choose between inter_pangenome alignment or "
                                                            "all_against_all alignment")
    # tmpdir.cleanup()


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("align", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_align(parser)
    return parser


def parser_align(parser):
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
    exclusive = parser.add_mutually_exclusive_group(required=True)
    exclusive.add_argument("--inter_pangenomes", action="store_true",
                           help="Align only gene families between pangenomes and not inside pangenome")
    exclusive.add_argument("--all_against_all", action="store_true",
                           help="Align gene families between pangenomes and intra-pangenome")
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
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--task", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")
    optional.add_argument("--threads_per_task", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_align(main_parser)
    common = main_parser._action_groups.pop(1)  # get the 'optional arguments' action group.
    common.title = "Common arguments"
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    main_parser._action_groups.append(common)
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
