#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from typing import Dict
from pathlib import Path
import tempfile
from concurrent.futures import ProcessPoolExecutor, CancelledError, Executor, Future
from multiprocessing import Manager, Lock

# installed libraries
from tqdm import tqdm

# local libraries
from panorama.format.read_binaries import check_pangenome_info
from panorama.pangenomes import Pangenome
from panorama.compare import init_db_lock


def context_mp(pangenome_name: str, pangenome_info: Dict[str, str]):
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, disable_bar=True)
    return True


def context_comparison(pangenome_path, lock: Lock, cpu: int = 1, disable_bar: bool = False):
    with ProcessPoolExecutor(max_workers=cpu, initializer=init_db_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenome_path), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome_name, pangenome_info in pangenome_path.items():
                future = executor.submit(context_mp, pangenome_name, pangenome_info)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            for future in futures:
                results = future.result()
    print(results)
