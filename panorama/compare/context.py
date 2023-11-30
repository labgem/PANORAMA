#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from typing import Dict, Union
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Lock

# installed libraries
from tqdm import tqdm
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import search_gene_context_in_pangenome

# local libraries
from panorama.utils import mkdir, init_lock
from panorama.pangenomes import Pangenome


def check_context_comparison(kwargs):
    if "sequences" not in kwargs and "family" not in kwargs:
        raise Exception("At least one of --sequences or --family option must be given")


def search_context_mp(pangenome_name: str, pangenome_info: Dict[str, Union[int, str]], output: Path, tmpdir: Path,
                      threads_per_task: int = 1, **kwargs):
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    search_gene_context_in_pangenome(pangenome=pangenome, output=output, tmpdir=tmpdir,
                                     cpu=threads_per_task, disable_bar=True, **kwargs)
    return True


def context_comparison(pangenome_path, lock: Lock, output: Path, tmpdir: Path, task: int = 1, threads_per_task: int = 1,
                       disable_bar: bool = False, force: bool = False, **kwargs):
    check_context_comparison(kwargs)
    mkdir(output, force)
    with ProcessPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
        with tqdm(total=len(pangenome_path), unit='pangenome', disable=disable_bar) as progress:
            futures = []
            for pangenome_name, pangenome_info in pangenome_path.items():
                mkdir(output/pangenome_name, force)
                future = executor.submit(search_context_mp, pangenome_name, pangenome_info, output, tmpdir,
                                         threads_per_task, **kwargs)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)
            for future in futures:
                results = future.result()
    print(results)


def context_comparison_parser(parser):
    required = parser._action_groups[2]
    required.add_argument('-S', '--sequences', type=Path, required=False,
                          help="Fasta file with the sequences of interest")
    required.add_argument('-F', '--family', type=Path, required=False,
                          help="List of family IDs of interest from the pan")
    optional = parser._action_groups[4]
    optional.add_argument('--no_defrag', required=False, action="store_true",
                          help="DO NOT Realign gene families to link fragments with"
                               "their non-fragmented gene family.")
    optional.add_argument('--identity', required=False, type=float, default=0.5,
                          help="min identity percentage threshold")
    optional.add_argument('--coverage', required=False, type=float, default=0.8,
                          help="min coverage percentage threshold")
    optional.add_argument("-t", "--transitive", required=False, type=int, default=4,
                          help="Size of the transitive closure used to build the graph. This indicates the number of "
                               "non related genes allowed in-between two related genes. Increasing it will improve "
                               "precision but lower sensitivity a little.")
    optional.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
                          help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
                               "will improve precision but lower sensitivity a lot.")
    parser._action_groups.pop(3)
