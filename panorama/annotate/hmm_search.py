#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import collections
import logging
from pathlib import Path
from typing import List, Tuple
from tqdm import tqdm
import tempfile
from concurrent.futures import ThreadPoolExecutor

# installed libraries
import pyhmmer
import pandas as pd
from numpy import nan

# local libraries
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily

res_col_names = ['families', 'Accession', 'protein_name', 'e_value',
                 'score', 'bias', 'secondary_name', 'Description']
meta_col_names = ["name", "accession", "path", "length", "protein_name", "secondary_name", "score_threshold",
                  "eval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
meta_dtype = {"accession": "string", "name": "string", "path": "string", "length": "int", "description": "string",
              "protein_name": "string", "secondary_name": "string", "score_threshold": "float",
              "eval_threshold": "float", "hmm_cov_threshold": "float", "target_cov_threshold": "float"}
seq_length = {}


def digit_gf_seq(pangenome: Pangenome, disable_bar: bool = False) -> List[pyhmmer.easel.DigitalSequence]:
    """
    Digitalised pangenome gene families sequences for hmmsearch

    Args:
        pangenome: Pangenome object with gene families
        disable_bar: Flag to disable progress bar

    Returns:
        List[pyhmmer.easel.DigitalSequence]: list of digitalised gene family sequences
    """
    global seq_length
    digit_gf_sequence = []
    logging.getLogger("PANORAMA").info("Begin to digitalized gene families sequences...")
    for family in tqdm(pangenome.gene_families, total=pangenome.number_of_gene_families, unit="gene families",
                       desc="Digitalized gene families sequences", disable=disable_bar):
        bit_name = family.name.encode('UTF-8')
        seq = pyhmmer.easel.TextSequence(name=bit_name,
                                         sequence=family.sequence if family.HMM is None
                                         else family.HMM.consensus.upper() + '*')
        seq_length[seq.name] = len(seq.sequence)
        digit_gf_sequence.append(seq.digitize(pyhmmer.easel.Alphabet.amino()))
    return digit_gf_sequence


def get_msa(tmpdir: Path):
    """
    Get the MSA profile

    Args:
        tmpdir: Path to the temporary directory

    Raises:
        NotImplementedError: Function is not implemented yet
    """
    _ = tempfile.TemporaryDirectory(prefix="msa_panorama", dir=tmpdir)
    logging.getLogger("PANORAMA").warning("In Dev")
    raise NotImplementedError


def profile_gf(gf: GeneFamily, msa_file: Path, msa_format: str = "afa", ):
    """
    Compute a profile for a gene family

    Args:
        gf: Gene family to profile
        msa_file: path to file containing msa
        msa_format: format used to write msa

    Raises:
        Exception: Problem to compute profile
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    if os.stat(msa_file).st_size == 0:
        logging.getLogger("PANORAMA").warning(f"{msa_file.absolute().as_posix()} is empty, so it's not readable. Pass to next file")
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
    """
    Create an HMM profile for each gene families

    Args:
        pangenome: Pangenome contaning gene families to profile
        msa_path: Path to file containing msa
        msa_format: format used to write msa
        tmpdir: Temporary directory for profiling
        threads: Number of available threads
        disable_bar: Flag to disable progress bar
    """
    if msa_path is None:
        get_msa(tmpdir)
    else:
        msa_file_list = list(msa_path.iterdir())
        with ThreadPoolExecutor(max_workers=threads) as executor:
            logging.getLogger("PANORAMA").info("Compute gene families HMM and profile")
            with tqdm(total=len(msa_file_list), unit='file', disable=disable_bar) as progress:
                futures = []
                for msa_file in msa_file_list:
                    gf = pangenome.get_gene_family(msa_file.stem)
                    future = executor.submit(profile_gf, gf, msa_file, msa_format)
                    future.add_done_callback(lambda p: progress.update())
                    futures.append(future)

                for future in futures:
                    future.result()


def read_hmms(hmm: Path, disable_bar: bool = False) -> Tuple[List[pyhmmer.plan7.HMM], pd.DataFrame]:
    """
    Read HMM file to create HMM object

    Args:
        hmm: Path to the HMM list file
        disable_bar: Flag to disable the progress bar

    Returns:
        Tuple[List[pyhmmer.plan7.HMM], pd.DataFrame]: List of HMM object

    Raises:
        Exception: Unexpected error occurred while reading HMM
    """
    hmms = []
    hmm_df = pd.read_csv(hmm, delimiter="\t", names=meta_col_names,
                         dtype=meta_dtype, header=0).set_index('accession')
    hmm_df['description'] = hmm_df["description"].fillna('unknown')
    logging.getLogger("PANORAMA").info("Begin to read HMM...")
    for hmm_path in tqdm(map(lambda x: Path(x), hmm_df["path"]), total=hmm_df.shape[0],
                         desc="Reading HMM", unit='HMM', disable=disable_bar):
        end = False
        try:
            hmm_file = pyhmmer.plan7.HMMFile(hmm_path)
        except Exception as error:
            logging.getLogger("PANORAMA").error(f"Problem reading HMM: {hmm_path}")
            raise error
        while not end:
            try:
                hmm = next(hmm_file)
            except StopIteration:
                end = True
            except Exception as error:
                raise Exception(f'Unexpected error on HMM file {hmm_path}, caused by {error}')
            else:
                hmms.append(hmm)
    return hmms, hmm_df


def annot_with_hmmsearch(hmm_list: List[pyhmmer.plan7.HMM], gf_sequences: List[pyhmmer.easel.DigitalSequence],
                         meta: pd.DataFrame = None, bit_cutoffs: str = None, threads: int = 1,
                         disable_bar: bool = False) -> List[Tuple[str, str, str, float, float, float, str, str]]:
    """
    Compute HMMer alignment between gene families sequences and HMM

    Args:
        hmm_list: List of HMM
        gf_sequences: List of digitalized gene families sequences
        meta: Metadata associate with HMM
        bit_cutoffs:
        threads: Number of available threads
        disable_bar:  Disable progress bar

    Returns:
         List[Tuple[str, str, str, float, float, float, str, str]]: Alignment results
    """
    res = []
    result = collections.namedtuple("Result", res_col_names)
    logging.getLogger("PANORAMA").info("Begin alignment of gene families to HMM")
    bar = tqdm(range(len(hmm_list)), unit="hmm", desc="Align gene families to HMM", disable=disable_bar)
    options = {"bit_cutoffs": bit_cutoffs}
    for top_hits in pyhmmer.hmmsearch(hmm_list, gf_sequences, cpus=threads, **options):
        for hit in top_hits:
            cog = hit.best_domain.alignment
            hmm_info = meta.loc[cog.hmm_accession.decode('UTF-8')]
            target_covery = ((max(cog.target_to, cog.target_from) - min(cog.target_to, cog.target_from)) /
                             seq_length[cog.target_name])
            hmm_covery = (max(cog.hmm_to, cog.hmm_from) - min(cog.hmm_to, cog.hmm_from)) / hmm_info["length"]
            add_res = False
            if ((target_covery >= hmm_info["target_cov_threshold"] or pd.isna(hmm_info["target_cov_threshold"])) and
                    (hmm_covery >= hmm_info["hmm_cov_threshold"] or pd.isna(hmm_info["hmm_cov_threshold"]))):
                if pd.isna(hmm_info['score_threshold']):
                    if hit.evalue < hmm_info['eval_threshold'] or pd.isna(hmm_info['eval_threshold']):
                        add_res = True
                else:
                    if hit.score > hmm_info['score_threshold']:
                        add_res = True
            if add_res:
                secondary_name = "" if pd.isna(hmm_info.secondary_name) else hmm_info.secondary_name
                res.append(result(hit.name.decode('UTF-8'), cog.hmm_accession.decode('UTF-8'), hmm_info.protein_name,
                                  hit.evalue, hit.score, hit.bias, secondary_name, hmm_info.description))
        bar.update()
    bar.close()
    return res


def annot_with_hmm(pangenome: Pangenome, hmms: List[pyhmmer.plan7.HMM], meta: pd.DataFrame = None, mode: str = "fast",
                   bit_cutoffs: str = None, threads: int = 1, disable_bar: bool = False) -> pd.DataFrame:
    """
    Takes a pangenome and a list of HMMs as input, and returns the best hit for each gene family in the pangenome.

    Args:
        pangenome: Pangenome with gene families
        hmms: Specify the hmm profiles to be used for annotation
        meta: Store the metadata of the hmms
        mode: Specify the number of threads to use for the hmm search
        bit_cutoffs: Specify the model-specific thresholding
        threads: Number of available threads
        disable_bar: bool: Disable the progress bar

    Returns:
        pd.DataFrame: A dataframe with the best results for each pangenome families

    Raises:
        ValueError: if the mode is not recognized. Possible value are fast and profile
        NotImplementedError: If the mode is profile. This option is not implemented yet

    Todo:
        Make the possibility to use the profile with the profile mode
    """
    if mode == "fast":
        gf_sequences = digit_gf_seq(pangenome, disable_bar=disable_bar)
    elif mode == "profile":
        raise NotImplementedError("Really sorry. It's not implemented yet, but be sure it's on schedule")
    else:
        raise ValueError("Unrecognized mode: {}".format(mode))
    logging.getLogger("PANORAMA").debug("Launch pyHMMer-HMMsearch")
    res = annot_with_hmmsearch(hmms, gf_sequences, meta, bit_cutoffs, threads, disable_bar)
    metadata_df = pd.DataFrame(res).fillna(nan)
    metadata_df.replace(to_replace='-', value=nan, inplace=True)
    metadata_df.replace(to_replace='', value=nan, inplace=True)
    logging.getLogger("PANORAMA").info("HMM search done.")
    return metadata_df
