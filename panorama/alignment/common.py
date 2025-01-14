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
from ppanggolin.formats.writeSequences import write_fasta_prot_fam_from_pangenome_file

# local libraries
from panorama.pangenomes import Pangenomes, Pangenome
from panorama.utils import init_lock, mkdir


def createdb(seq_files: List[Path], tmpdir: Path, db_type: int = 0, keep_tmp: bool = False) -> Path:
    """
    Create a MMseqs2 sequence database with the given fasta file

    Args:
        seq_files: List of fasta file
        tmpdir: temporary directory
        db_type: type of MMSeqs2 database (Default 0)
        keep_tmp: Whether to keep the temporary directory after execution. (Defaults False).

    Returns:
        DB file
    """
    logging.getLogger('PANORAMA').debug('Create MMSeqs2 database')
    seqdb = tempfile.NamedTemporaryFile(mode="w", dir=tmpdir, delete=keep_tmp)
    cmd = ["mmseqs", "createdb"] + list(map(Path.as_posix, map(Path.absolute, seq_files))) + \
          [seqdb.name, "--dbtype", str(db_type)]
    logging.getLogger("PANORAMA").debug(" ".join(cmd))
    subprocess.run(cmd, stdout=subprocess.DEVNULL)
    logging.getLogger('PANORAMA').debug(f'MMSeqs2 database created here: {seqdb.name}')
    return Path(seqdb.name)


def write_pangenomes_families_sequences(pangenomes: Pangenomes, tmpdir: Path, threads: int = 1,
                                        lock: Lock = None, disable_bar: bool = False) -> Dict[str, Path]:
    """
    Write pangenomes families sequences

    Args:
        pangenomes: Pangenomes objects containing pangenome
        tmpdir: Temporary directory to write sequences
        threads: Number of available threads for each worker. (Defaults 1).
        lock: Lock object to write sequences in multithreading
        disable_bar: Disable progressive bar. (Defaults False).

    Returns:
        Dictionary with for each pangenome the path to gene families sequences
    """

    def write_protein_families_sequences(pan: Pangenome, **kwargs) -> Tuple[str, Path]:
        """Wrapper to write protein families sequences in multithreading

        Args:
            pan: Pangenome object to get gene families sequences
            **kwargs: Additional arguments to pass to called function

        Returns:
            Name of the pangenome with the path to the gene families sequences

        todo add in ppanggolin the return of the output path and also
        add in the begin of sequence name the pangenome_name to be sure that there is no duplicate families
        between pangenomes
        """
        write_fasta_prot_fam_from_pangenome_file(pan.file, **kwargs)
        return pan.name, kwargs["output"] / "all_protein_families.faa.gz"

    logging.getLogger("PANORAMA").info("Writing pangenomes families sequences...")
    with ThreadPoolExecutor(max_workers=threads, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenomes), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome in pangenomes:
                args = {"output": mkdir(tmpdir / f"{pangenome.name}"), "family_filter": "all",
                        "compress": True, "disable_bar": True}
                future = executor.submit(write_protein_families_sequences, pangenome, **args)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            pangenomes2families_sequences = {}
            for future in futures:
                result = future.result()
                pangenomes2families_sequences[result[0]] = result[1]
    return pangenomes2families_sequences
