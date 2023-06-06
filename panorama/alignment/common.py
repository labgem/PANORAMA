#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from pathlib import Path
import tempfile
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import Manager, Lock
import subprocess

# installed libraries
from tqdm import tqdm
from ppanggolin.align.alignOnPang import write_gene_fam_sequences

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import init_lock


def createdb(seq_file: Path, tmpdir: tempfile.TemporaryDirectory, other_seq: List[Path] = None) -> Path:
    """
    Create a MMseqs2 sequence database with the given fasta file

    :param seq_file: Fasta file
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


def get_gf_pangenome(pangenome: Pangenome, tmpdir: tempfile.TemporaryDirectory, create_db: bool = False) -> Tuple[
    str, Path]:
    logging.debug(f"Creating {pangenome.name} gene families sequence file...")
    gf_fasta_path = Path(f"{tmpdir.name}/{pangenome.name}_gf.fna")
    with open(gf_fasta_path, "w") as gf_fasta:
        write_gene_fam_sequences(pangenome, gf_fasta)
        if create_db:
            logging.debug(f"Creating {pangenome.name} gene families sequence file...")
            pan_db = createdb(gf_fasta_path, tmpdir)
            logging.debug(f"{pangenome.name} gene families database created")
            return pangenome.name, pan_db
        else:
            logging.debug(f"{pangenome.name} gene families fasta created")
            return pangenome.name, gf_fasta_path


def get_gf_pangenomes(pangenomes: Pangenomes, create_db: bool, lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                      threads: int = 1, disable_bar: bool = False) -> Dict[Pangenome, Path]:
    logging.debug("Begin create pangenomes gene families database...")
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome in pangenomes:
                future = executor.submit(get_gf_pangenome, pangenome, tmpdir, create_db)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            pangenomes_db = {}
            for future in futures:
                results = future.result()
                pangenomes_db[results[0]] = results[1]
    logging.debug("All pangenomes gene families database created")
    return pangenomes_db
