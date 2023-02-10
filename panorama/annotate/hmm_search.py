#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import collections
import logging
import re
from pathlib import Path
from typing import List, Generator, Tuple
from tqdm import tqdm
import tempfile
from concurrent.futures import ThreadPoolExecutor

# installed libraries
import pyhmmer
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily

res_col_names = ['Gene_family', 'Accession', 'protein_name', 'e_value',
                 'score', 'bias', 'secondary_name', 'Description']
meta_col_names = ["accession", "hmm_name", "protein_name", "secondary_name", "score_threshold",
                  "eval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
meta_dtype = {"accession": "string", "hmm_name": "string", "protein_name": "string", "secondary_name": "string",
              "score_threshold": "float", "eval_threshold": "float", "hmm_cov_threshold": "float",
              "target_cov_threshold": "float", "description": "string"}
seq_length = {}
hmm_length = {}


def digit_gf_seq(pangenome: Pangenome, disable_bar: bool = False) -> List[pyhmmer.easel.DigitalSequence]:
    """Digitalised pangenome gene families sequences for hmmsearch

    :param pangenome: Pangenome object with gene families
    :param disable_bar: Disable bar

    :return: list of digitalised gene family sequences
    """
    global seq_length
    digit_gf_sequence = []
    logging.getLogger().info("Digitalized gene families sequences")
    for family in tqdm(pangenome.gene_families, unit="gene families", disable=disable_bar):
        bit_name = family.name.encode('UTF-8')
        seq = pyhmmer.easel.TextSequence(name=bit_name,
                                         sequence=family.sequence if family.HMM is None
                                         else family.HMM.consensus.upper() + '*')
        seq_length[seq.name] = len(seq.sequence)
        digit_gf_sequence.append(seq.digitize(pyhmmer.easel.Alphabet.amino()))
    return digit_gf_sequence


def get_msa(tmpdir: Path):
    _ = tempfile.TemporaryDirectory(prefix="msa_panorama", dir=tmpdir)
    logging.getLogger().warning("In Dev")
    raise NotImplementedError


def profile_gf(gf: GeneFamily, msa_file: Path, msa_format: str = "afa", ):
    """ Compute a profile for a gene family

    :param gf: Gene family to profile
    :param msa_file: path to file containing msa
    :param msa_format: format used to write msa

    :raise Exception: Problem to compute profile
    """
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
                gf._hmm, gf.profile, gf.optimized_profile = builder.build_msa(msa, background)
        except Exception as error:
            raise Exception(f"The following error happened with file {msa_file} : {error}")


def profile_gfs(pangenome: Pangenome, msa_path: Path = None, msa_format: str = "afa",
                tmpdir: Path = Path(tempfile.gettempdir()), threads: int = 1, disable_bar: bool = False):
    """Create an HMM profile for each gene families

    :param pangenome: Pangenome contaning gene families to profile
    :param msa_path: path to file containing msa
    :param msa_format: format used to write msa
    :param threads: Number of available threads
    :param disable_bar: Disable progress bar
    """
    if msa_path is None:
        get_msa(tmpdir)
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


def read_metadata(metadata: Path) -> pd.DataFrame:
    """ Read metadata associate with HMM

    :param metadata: PAth to metadata

    :return: metadata dataframe and dictionnary associate hmm name with protein and secondary name
    """
    metadata_df = pd.read_csv(metadata, delimiter="\t", names=meta_col_names,
                              dtype=meta_dtype, header=0).set_index('accession')
    metadata_df['description'] = metadata_df["description"].fillna('unknown')
    return metadata_df  # , hmm_to_metaname


def read_hmm(hmm_dir: Generator[Path, None, None], disable_bar: bool = False) -> List[pyhmmer.plan7.HMM]:
    """Read HMM file to create HMM object

    :param hmm_dir: Directory with the HMM
    :param disable_bar: Disable progress bar

    :return: List of HMM object
    """
    hmms = []
    hmm_list_path = list(hmm_dir)
    logging.getLogger().info("Read HMM")
    for hmm_path in tqdm(hmm_list_path, total=len(hmm_list_path), unit='HMM', disable=disable_bar):
        # if hmm_path.stem == "RM__Type_I_MTases":
        #     print("pika")
        end = False
        hmm_file = pyhmmer.plan7.HMMFile(hmm_path)
        while not end:
            try:
                hmm = next(hmm_file)
            except StopIteration:
                end = True
            except Exception as error:
                raise Exception(f'Unexpected error on HMM file {hmm_path}, caused by {error}')
            else:
                try:
                    with open(hmm_path, 'r') as file:
                        line = file.readline()
                        while not re.search('LENG', line):
                            line = file.readline()
                        hmm_length[hmm.name] = int(line.split(' ')[-1])
                except ValueError as val_error:
                    raise ValueError(f'{hmm_path} {val_error}')
                except Exception as error:
                    raise Exception(f'Unexpected error to get HMM length on file {hmm_path}.'
                                    f'Error caused by {error}')
                hmms.append(hmm)
    return hmms


def annot_with_hmmsearch(hmm_list: List[pyhmmer.plan7.HMM], gf_sequences: List[pyhmmer.easel.DigitalSequence],
                         hmm_meta: pd.DataFrame = None, threads: int = 1,
                         disable_bar: bool = False) -> List[Tuple[str, str, str, float, float, float, str, str]]:
    """Compute HMMer alignment between gene families sequences and HMM

    :param hmm_list: List of HMM
    :param gf_sequences: List of gene families sequences
    :param hmm_meta: Metadata associate with HMM
    :param threads: Number of available threads
    :param disable_bar: Disable progress bar

    :return: Alignment results
    """
    res = []
    result = collections.namedtuple("Result", res_col_names)
    logging.getLogger().info("Align gene families to HMM")
    bar = tqdm(range(len(hmm_list)), unit="hmm", disable=disable_bar)
    bit_cutoffs = "gathering" if hmm_meta is None else None
    for top_hits in pyhmmer.hmmsearch(hmm_list, gf_sequences, cpus=threads, bit_cutoffs=bit_cutoffs):
        for hit in top_hits:
            cog = hit.best_domain.alignment
            target_covery = (max(cog.target_to, cog.target_from) - min(cog.target_to, cog.target_from)) / seq_length[
                cog.target_name]
            hmm_covery = (max(cog.hmm_to, cog.hmm_from) - min(cog.hmm_to, cog.hmm_from)) / hmm_length[cog.hmm_name]
            if hmm_meta is not None:
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
                    res.append(
                        result(hit.name.decode('UTF-8'), cog.hmm_accession.decode('UTF-8'), hmm_info.protein_name,
                               hit.evalue, hit.score, hit.bias, hmm_info.secondary_name, hmm_info.description))
            else:
                res.append(result(hit.name.decode('UTF-8'), cog.hmm_accession.decode('UTF-8'),
                                  cog.hmm_name.decode('UTF-8'), hit.evalue, hit.score, hit.bias, "",
                                  "unknow" if hit.description is None else hit.description.decode('UTF-8')))
        bar.update()
    bar.close()
    return res


def annot_with_hmm(pangenome: Pangenome, hmm_path: Path, meta: Path = None,
                   mode: str = 'fast', msa: Path = None, tmpdir: Path = Path(tempfile.gettempdir()),
                   threads: int = 1, disable_bar: bool = False) -> pd.DataFrame:
    """ Launch hmm search against pangenome gene families and HMM

    :param pangenome: Pangenome with gene families
    :param hmm_path: Path to one file with multiple HMM or a directory with one HMM by file
    :param meta: Path to metadata associate with HMM
    :param mode: alignment method used
    :param msa: Path to msa results
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param disable_bar: allow to disable progress bar

    :return: dataframe with best result for each pangenome families
    """
    logging.getLogger().info("Begin HMM searching")
    if mode == 'profile':
        profile_gfs(pangenome, msa, tmpdir=tmpdir, threads=threads, disable_bar=disable_bar)
        logging.getLogger().debug("Gene families HMM and profile computation... Done")
    else:
        if mode != 'fast':
            raise ValueError(f"{mode} unrecognized mode")
    gf_sequences = digit_gf_seq(pangenome, disable_bar=disable_bar)
    if meta is not None:
        metadata = read_metadata(meta)
    else:
        metadata = None
    # Get list of HMM with Plan7 data model
    hmms = read_hmm(hmm_dir=hmm_path.iterdir(), disable_bar=disable_bar)
    res = annot_with_hmmsearch(hmms, gf_sequences, metadata, threads, disable_bar)
    return pd.DataFrame(res)
