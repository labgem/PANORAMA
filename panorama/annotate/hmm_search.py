#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import collections
import logging
from pathlib import Path
from typing import List, Generator, Tuple, Union
from tqdm import tqdm
import tempfile
from concurrent.futures import ThreadPoolExecutor

# installed libraries
import pyhmmer
import pandas as pd

# local libraries
from geneFamily import GeneFamily
from pangenomes import Pangenome
from utils import path_is_dir, path_is_file

res_col_names = ['Gene_family', 'Accession', 'protein_name', 'e_value',
                 'score', 'bias', 'secondary_name', 'Description']
meta_col_names = ["accession", "hmm_name", "protein_name", "secondary_name", "score_threshold",
                  "eval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
meta_dtype = {"accession": "string", "hmm_name": "string", "protein_name": "string", "secondary_name": "string",
              "score_threshold": "float", "eval_threshold": "float", "hmm_cov_threshold": "float",
              "target_cov_threshold": "float", "description": "string"}


def digit_gf_seq(pangenome: Pangenome, disable_bar: bool = False) -> List[pyhmmer.easel.DigitalSequence]:
    """Digitalised pangenome gene families sequences for hmmsearch

    :param pangenome: Pangenome object with gene families

    :return: list of digitalised gene family sequences
    """
    digit_gf_sequence = []
    logging.getLogger().info("Digitalized gene families sequences")
    for family in tqdm(pangenome.gene_families, unit="gene families", disable=disable_bar):
        bit_name = family.name.encode('UTF-8')
        seq = pyhmmer.easel.TextSequence(name=bit_name,
                                         sequence=family.sequence if family.hmm is None
                                         else family.hmm.consensus.upper() + '*')
        digit_gf_sequence.append(seq.digitize(pyhmmer.easel.Alphabet.amino()))
    return digit_gf_sequence


def get_msa():
    pass

def profile_gf(gf: GeneFamily, msa_file: Path,  msa_format: str = "afa", ):
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    if os.stat(msa_file).st_size == 0:
        logging.getLogger().warning(f"{msa_file.absolute().as_posix()} is empty, so it's not readable."
                                    f"Pass to next file")
    else:
        try:
            with pyhmmer.easel.MSAFile(msa_file.absolute().as_posix(), format=msa_format,
                                       digital=True, alphabet=alphabet) as easel_msa:
                msa = next(easel_msa)
                msa.name = gf.name.encode('UTF-8')
                msa.accession = f"PAN{gf.ID}".encode('UTF-8')
                gf.hmm, gf.profile, gf.optimized_profile = builder.build_msa(msa, background)
        except Exception as error:
            raise Exception(f"The following error happened with file {msa_file} : {error}")


def profile_gfs(pangenome: Pangenome, msa_path: Path = None, msa_format: str = "afa", threads: int = 1, disable_bar: bool = False):
    """Create an HMM profile for each gene families"""
    if msa_path is None:
        get_msa()
    else:
        msa_file_list = list(msa_path.iterdir())
        with ThreadPoolExecutor(max_workers=threads) as executor:
            logging.getLogger().info("Compute gene families HMM and profile")
            with tqdm(total=len(msa_file_list), unit='file', disable=disable_bar) as progress:
                futures = []
                for msa_file in msa_file_list:
                    gf = pangenome.get_gene_family(msa_file.stem)
                    future = executor.submit(profile_gf, gf, msa_file, msa_format)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)

                for future in futures:
                    future.result()


def read_metadata(metadata: Path):
    metadata_df = pd.read_csv(metadata, delimiter="\t", names=meta_col_names,
                              dtype=meta_dtype, header=0).set_index('accession')
    metadata_df['description'] = metadata_df["description"].fillna('unknown')
    hmm_full_name = collections.namedtuple("HMM_full_name", ['protein_name', 'secondary_name'])
    hmm_to_metaname = {}
    for index, row in metadata_df.iterrows():
        hmm_to_metaname[row[0]] = hmm_full_name(row[1], row[2])
    return metadata_df, hmm_to_metaname


def read_hmm(hmm_dir: Generator, disable_bar: bool = False):
    hmms = []
    hmm_list_path = list(hmm_dir)
    logging.getLogger().info("Read HMM")
    for hmm_file in tqdm(hmm_list_path, total=len(hmm_list_path), unit='HMM', disable=disable_bar):
        try:
            hmm = next(pyhmmer.plan7.HMMFile(hmm_file))
        except ValueError as val_error:
            raise ValueError(f'{hmm_file} {val_error}')
        except Exception as error:
            raise Exception(f'Unexpected error on HMM file {hmm_file}, caused by {error}')
        else:
            hmms.append(hmm)
    return hmms


def split_hmm_file(hmm_path: Path, tmpdir: tempfile.TemporaryDirectory) -> List[Path]:
    """ Split a big HMM file into multiple temporary files with one HMM by file for multithreading

    :param hmm_path: Path to HMM file
    :param tmpdir: Temporary directory to store HMM file

    :return: list of path to temporary HMM file
    """
    hmm_list = []
    with open(hmm_path.absolute().as_posix(), 'r') as hmm_file:
        counter = 0
        lines = [hmm_file.readline()]
        for line in hmm_file.readlines():
            if line.startswith('HMMER'):
                hmm_list.append(Path(f"{tmpdir.name}/HMM_{counter}.hmm"))
                new_hmm_file = open(hmm_list[-1].as_posix(), 'w')
                new_hmm_file.write("\n".join(lines))
                counter += 1
                lines = [line]
            else:
                lines.append(line)
        hmm_list.append(Path(f"{tmpdir.name}/HMM_{counter}.hmm"))
        new_hmm_file = open(hmm_list[-1].as_posix(), 'w')
        new_hmm_file.write("\n".join(lines))
    return hmm_list


def annot_with_hmmsearch(hmm_list: List[pyhmmer.plan7.HMM], gf_sequences: List[pyhmmer.easel.DigitalSequence],
                         hmm_meta: pd.DataFrame = None, threads: int = 1,
                         disable_bar: bool = False) -> List[Tuple[str, str, str, float, float, float, str, str]]:
    res = []
    result = collections.namedtuple("Result", res_col_names)
    logging.getLogger().info("Align gene families to HMM")
    bar = tqdm(range(len(hmm_list)), unit="hmm", disable=disable_bar)
    for top_hits in pyhmmer.hmmsearch(hmm_list, gf_sequences, cpus=threads):
        for hit in top_hits:
            cog = hit.best_domain.alignment
            target_covery = (max(cog.target_to, cog.target_from) - min(cog.target_to, cog.target_from))/len(cog.target_sequence)
            hmm_covery = (max(cog.hmm_to, cog.hmm_from) - min(cog.hmm_to, cog.hmm_from))/len(cog.hmm_sequence)
            hmm_info = hmm_meta.loc[cog.hmm_accession.decode('UTF-8')]
            add_res = False
            if target_covery >= hmm_info["target_cov_threshold"] and hmm_covery >= hmm_info["hmm_cov_threshold"]:
                if pd.isna(hmm_info['score_threshold']):
                    if hit.evalue < hmm_info['eval_threshold']:
                        add_res = True
                else:
                    if hit.score > hmm_info['score_threshold']:
                        add_res = True
            if add_res:
                res.append(result(hit.name.decode('UTF-8'), cog.hmm_accession.decode('UTF-8'), hmm_info.protein_name,
                                  hit.evalue, hit.score, hit.bias, hmm_info.secondary_name, hmm_info.description))
        bar.update()
    bar.close()
    return res


def annot_with_hmm(pangenome: Pangenome, hmm_path: Path, meta: Path = None, mode: str = 'fast', msa: Path = None,
                   tmpdir: Path = Path(tempfile.gettempdir()), threads: int = 1, disable_bar: bool = False):
    """ Launch hmm search against pangenome gene families and HMM

    :param pangenome: Pangenome with gene families
    :param hmm_path: Path to one file with multiple HMM or a directory with one HMM by file
    :param meta:
    :param mode:
    :param msa:
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param disable_bar: allow to disable progress bar

    :return: dataframe with best result for each pangenome families
    """
    logging.getLogger().info("Begin HMM searching")
    if mode == 'profile':
        profile_gfs(pangenome, msa, threads=threads, disable_bar=disable_bar)
        logging.getLogger().debug("Gene families HMM and profile computation... Done")
    else:
        if mode != 'fast':
            raise ValueError(f"{mode} unrecognized mode")
    gf_sequences = digit_gf_seq(pangenome, disable_bar=disable_bar)
    tmp_dir = tempfile.TemporaryDirectory(prefix="hmm_panorama", dir=tmpdir)
    metadata, hmm_to_metaname = read_metadata(meta)
    # Get list of HMM with Plan7 data model
    if path_is_file(hmm_path):
        split_hmm_file(hmm_path, tmp_dir)
        hmms = read_hmm(hmm_dir=Path(tmp_dir.name).iterdir(), disable_bar=disable_bar)

    elif path_is_dir(hmm_path):
        hmms = read_hmm(hmm_dir=hmm_path.iterdir(), disable_bar=disable_bar)
    else:
        raise Exception("Unexpected error")
    res = annot_with_hmmsearch(hmms, gf_sequences, metadata, threads, disable_bar)
    return pd.DataFrame(res)


if __name__ == "__main__":
    # default libraries
    import argparse
    from time import time
    # installed libraries
    from ppanggolin.formats.readBinaries import check_pangenome_info
    # local libraries
    from panorama.utils import check_log, set_verbosity_level


    def plot_res(pred_size_max, nb_pred, run_time):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        fig.subplots_adjust(right=0.75)

        twin1 = ax.twinx()

        p1, = ax.plot(range(1, pred_size_max + 1), nb_pred, 'r.', label="Nombre de prédiction total")
        p2, = twin1.plot(range(1, pred_size_max + 1), run_time, 'bx', label="Temps de calcul")

        ax.set_xlabel("Nombre maximum de prédiction par familles")
        ax.set_ylabel("Nombre de prédiction total")
        twin1.set_ylabel("Temps de calcul (s)")

        tkw = dict(size=4, width=1.5)
        ax.tick_params(axis='y', colors=p1.get_color(), **tkw)
        twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        ax.tick_params(axis='x', **tkw)

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        # Put a legend to the right of the current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), handles=[p1, p2],
                  fancybox=True, shadow=True, ncol=2)

        plt.show()


    def plot_stat(pred_size_max, res_stat):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        p1 = ax.plot(range(1, pred_size_max + 1), [res.mean for res in res_stat],
                     'ro', label="Moyenne")
        p2 = ax.plot(range(1, pred_size_max + 1), [res.min for res in res_stat],
                     'bo', label="Minimum")
        p3 = ax.plot(range(1, pred_size_max + 1), [res.max for res in res_stat],
                     'go', label="Maximum")

        ax.errorbar(range(1, pred_size_max + 1), [res.mean for res in res_stat],
                    yerr=[res.std for res in res_stat], fmt="+", capsize=6)

        ax.set_xlabel("Nombre maximum de prédiction par familles")
        ax.set_ylabel("Nombre de prédiction par famille")

        # Put a legend to the right of the current axis
        ax.legend(loc='upper center', fancybox=True, shadow=True, ncol=2)

        plt.show()


    def test_hmm_search(args: argparse.Namespace):
        pangenome = Pangenome(name="pangenome")
        pangenome.add_file(pangenome_file=args.pangenome)
        check_pangenome_info(pangenome, need_families=True)
        res_df = annot_with_hmm(pangenome=pangenome, hmm_path=args.hmm, method=args.method, tmpdir=args.tmpdir,
                                meta=args.meta, max_prediction=args.max_prediction,
                                threads=args.threads, disable_bar=args.disable_prog_bar)
        res_df.to_csv(path_or_buf="hmm_res.tsv", sep="\t")
        return res_df


    main_parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
        formatter_class=argparse.RawTextHelpFormatter)

    req = main_parser.add_argument_group(title="Required for test")
    req.add_argument("-p", "--pangenome", required=True, type=Path,
                     help="Pangenome.h5 file to test annotation with hmmer")
    req.add_argument("--hmm", required=True, type=Path,
                     help="HMM file to test annotation with hmmer")
    req.add_argument("--meta", required=True, type=Path,
                     help="Metadata link to HMM with protein name, description and cutoff")
    opt = main_parser.add_argument_group(title="Optional argument")
    opt.add_argument("--tmpdir", required=False, type=str, default=Path(tempfile.gettempdir()),
                     help="directory for storing temporary files")
    opt.add_argument("--threads", required=False, type=int, default=1)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    main_args = main_parser.parse_args()
    set_verbosity_level(main_args)

    res_df = test_hmm_search(main_args)

    # run_time = []
    # nb_pred = []
    # stat = collections.namedtuple("Stat", ['min', 'max', 'mean', "std"])
    # res_stat = []
    # max_pred = main_args.max_prediction
    # for i in range(1, max_pred + 1):
    #     print(i)
    #     b_time = time()
    #     main_args.max_prediction = i
    #     res_df = test_hmm_search(main_args)
    #     res_vc = res_df['Gene_family'].value_counts()
    #     res_stat.append(stat(res_vc.min(), res_vc.max(), res_vc.mean(), res_vc.std()))
    #     nb_pred.append(res_vc.sum())
    #     run_time.append(time() - b_time)
    # print(nb_pred, run_time)
    # plot_res(max_pred, nb_pred, run_time)
    # plot_stat(max_pred, res_stat)
