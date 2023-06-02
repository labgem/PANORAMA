#!/usr/bin/env python3
# coding:utf-8

# default libraries
from __future__ import annotations
from typing import Dict
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Lock

# installed libraries
from tqdm import tqdm
import pandas as pd
from ppanggolin.utils import restricted_float
from ppanggolin.context.searchGeneContext import search_gene_context_in_pangenome

# local libraries
from panorama.utils import mkdir, init_lock, load_multiple_pangenomes
from panorama.pangenomes import Pangenome


def check_run_context_arguments(kwargs):
    if "sequences" not in kwargs and "family" not in kwargs:
        raise Exception("At least one of --sequences or --family option must be given")


def search_context_mp(pangenome_name: str, pangenome_info: Dict[str, str], output: Path, tmpdir: Path,
                      threads_per_task: int = 1, **kwargs):
    pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
    pangenome.add_file(pangenome_info["path"])
    search_gene_context_in_pangenome(pangenome=pangenome, output=output, tmpdir=tmpdir,
                                     cpu=threads_per_task, disable_bar=True, **kwargs)
    return True

def parse_cluster_family(cluster_family_file, family_of_interest=None):

    family_to_cluster = {}
    with open(cluster_family_file) as fh:
        for line in fh:
            cluster_id, familly = line.rstrip().split('\t')

            family_to_cluster[familly] = cluster_id
    
    return family_to_cluster

def parse_context_results(contexts_result_file_list):
    
    pangenome_to_context_file = {}
    
    with open(contexts_result_file_list) as fh:
        for line in fh:
            pan_name, context_file = line.rstrip().split('\t')
            pangenome_to_context_file[pan_name] = context_file

    return pangenome_to_context_file

def parse_context_table(context_table):

    return pd.read_csv(context_table, sep='\t')

    

def context_comparison(pangenome_to_path, contexts_results, familly_clusters, lock: Lock, output: Path, 
                       tmpdir: Path, task: int = 1, threads_per_task: int = 1,
                       disable_bar: bool = False, force: bool = False, **kwargs):
    
    mkdir(output, force)

    need_info = {"need_families":True,}

    load_multiple_pangenomes(pangenome_to_path, task, disable_bar, lock, need_info)



    if contexts_results:
        # Parse the given context results
        pangenome_to_context_file =  parse_context_results(contexts_results)

        pang_to_context_df = {pan:pd.read_csv(context, sep='\t')  for pan, context in pangenome_to_context_file.items()}


        

    else:
        # run ppanggolin context in parallel
        raise NotImplementedError
    
        check_run_context_arguments(kwargs)

        with ProcessPoolExecutor(max_workers=task, initializer=init_lock, initargs=(lock,)) as executor:
            with tqdm(total=len(pangenome_to_path), unit='pangenome', disable=disable_bar) as progress:
                futures = []
                for pangenome_name, pangenome_info in pangenome_to_path.items():
                    mkdir(output/pangenome_name, force)
                    future = executor.submit(search_context_mp, pangenome_name, pangenome_info, output, tmpdir,
                                            threads_per_task, **kwargs)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)
                for future in futures:
                    results = future.result()
        
    if familly_clusters:
        # Parse the given cluster family results 
        pass
       
    else:
        # run cluster familly
        # familly_clusters_file = panorama cluster function
        raise NotImplementedError
         
    
    


def context_comparison_parser(parser):

    use_context_arg = parser.add_argument_group("Use already computed contexts arguments")

    use_context_arg.add_argument('--context_results', type=Path, required=False,
                        help="Tsv file with two columns: name of pangenome and path to the corresponding context results.")
    
    use_context_arg.add_argument('--familly_clusters', type=Path, required=False,
                        help="A tab-separated file listing the cluster names, the family IDs,")

    # TODO: use context parser of ppanggolin to have correct args and prevent duplication

    # run_context_arg = parser.add_argument_group("PPanGGOLiN context arguments")

    # # required = parser._action_groups[2]
    # run_context_arg.add_argument('-S', '--sequences', type=Path, required=False,
    #                       help="Fasta file with the sequences of interest")
    
    # run_context_arg.add_argument('-F', '--family', type=Path, required=False,
    #                       help="List of family IDs of interest from the pan")
    
    # # optional = parser._action_groups[4]

    # run_context_arg.add_argument('--no_defrag', required=False, action="store_true",
    #                       help="DO NOT Realign gene families to link fragments with"
    #                            "their non-fragmented gene family.")
    
    # run_context_arg.add_argument('--identity', required=False, type=float, default=0.5,
    #                       help="min identity percentage threshold")
    
    # run_context_arg.add_argument('--coverage', required=False, type=float, default=0.8,
    #                       help="min coverage percentage threshold")
    
    # run_context_arg.add_argument("-t", "--transitive", required=False, type=int, default=4,
    #                       help="Size of the transitive closure used to build the graph. This indicates the number of "
    #                            "non related genes allowed in-between two related genes. Increasing it will improve "
    #                            "precision but lower sensitivity a little.")
    
    # run_context_arg.add_argument("-s", "--jaccard", required=False, type=restricted_float, default=0.85,
    #                       help="minimum jaccard similarity used to filter edges between gene families. Increasing it "
    #                            "will improve precision but lower sensitivity a lot.")


