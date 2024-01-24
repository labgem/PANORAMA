#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
import logging
from pathlib import Path
import tempfile
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Lock
import subprocess

# installed libraries
from tqdm import tqdm
from ppanggolin.align.alignOnPang import write_gene_fam_sequences

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import init_lock


def createdb(seq_files: List[Path], tmpdir: tempfile.TemporaryDirectory, db_type: int = 0) -> Path:
    """
    Create a MMseqs2 sequence database with the given fasta file

    :param seq_files: List of fasta file
    :param tmpdir: temporary directory
    :param db_type: type of MMSeqs2 database

    :return: DB file
    """
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir.name, delete=False)
    cmd = ["mmseqs", "createdb"] + list(map(Path.as_posix, map(Path.absolute, seq_files))) + \
          [seqdb.name, "--dbtype", str(db_type)]
    logging.getLogger("PANORAMA").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    return Path(seqdb.name)


def get_gf_pangenome(pangenome: Pangenome, tmpdir: tempfile.TemporaryDirectory,
                     create_db: bool = False) -> Tuple[str, Path]:
    """Get gene families sequences in a pangenome and if ask create an MMSeqs2 database

    :param pangenome: Pangenome object containing gene families
    :param tmpdir: temporary directory
    :param create_db: boolean to create a MMSeqs2 database

    :return: Pangenome name and path to the sequences or the database if asked
    """
    logging.getLogger("PANORAMA").debug(f"Creating {pangenome.name} gene families sequence file...")
    gf_fasta_path = Path(f"{tmpdir.name}/{pangenome.name}_gf.fna")
    with open(gf_fasta_path, "w") as gf_fasta:
        write_gene_fam_sequences(pangenome, gf_fasta)
        if create_db:
            logging.getLogger("PANORAMA").debug(f"Creating {pangenome.name} gene families sequence file...")
            pan_db = createdb([gf_fasta_path], tmpdir)
            logging.getLogger("PANORAMA").debug(f"{pangenome.name} gene families database created")
            return pangenome.name, pan_db
        else:
            logging.getLogger("PANORAMA").debug(f"{pangenome.name} gene families fasta created")
            return pangenome.name, gf_fasta_path


def get_gf_pangenomes(pangenomes: Pangenomes, create_db: bool, lock: Lock, tmpdir: tempfile.TemporaryDirectory,
                      threads: int = 1, disable_bar: bool = False) -> Dict[str, Path]:
    """Get gene families sequences from pangenomes and if ask create an MMSeqs2 database for each pangenome gene families

    :param pangenomes: Pangenomes with gene families
    :param create_db: boolean to create or not database of gene families
    :param lock: Global lock for multiprocessing execution
    :param tmpdir: Temporary directory for MMSeqs2
    :param threads: Number of available threads
    :param disable_bar: Disable progressive bar

    :return: Dictionary with pangenome name and path to the sequences or database
    """
    logging.getLogger("PANORAMA").debug("Begin create pangenomes gene families database...")
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome in pangenomes:
                future = executor.submit(get_gf_pangenome, pangenome, tmpdir, create_db)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            pangenomes_gf = {}
            for future in futures:
                results = future.result()
                pangenomes_gf[results[0]] = results[1]
    logging.getLogger("PANORAMA").debug("All pangenomes gene families database created")
    return pangenomes_gf
