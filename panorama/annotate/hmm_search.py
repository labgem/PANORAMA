#!/usr/bin/env python3
# coding:utf-8

# default libraries
import os
import sys
from collections import defaultdict, namedtuple
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Union

from tqdm import tqdm
import tempfile
from concurrent.futures import ThreadPoolExecutor

# installed libraries
import psutil
from pyhmmer.easel import (TextSequence, DigitalSequence, DigitalSequenceBlock,
                           SequenceBlock, Alphabet, MSAFile, SequenceFile)
from pyhmmer.plan7 import Builder, Background, HMM, HMMFile, Hit, TopHits
from pyhmmer import hmmsearch, hmmscan, hmmpress
import pandas as pd
from numpy import nan
from ppanggolin.formats.writeSequences import write_gene_protein_sequences
from ppanggolin.formats.writeMSA import write_msa_files

# local libraries
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily

res_col_names = ['families', 'Accession', 'protein_name', 'e_value',
                 'score', 'bias', 'i_e_value', 'secondary_name', 'Description']
meta_col_names = ["name", "accession", "path", "length", "protein_name", "secondary_name", "score_threshold",
                  "eval_threshold", "ieval_threshold", "hmm_cov_threshold", "target_cov_threshold", "description"]
meta_dtype = {"accession": "string", "name": "string", "path": "string", "length": "int", "description": "string",
              "protein_name": "string", "secondary_name": "string", "score_threshold": "float",
              "eval_threshold": "float", "ieval_threshold": "float", "hmm_cov_threshold": "float",
              "target_cov_threshold": "float"}


def digit_gene_sequences(pangenome: Pangenome, threads: int = 1, tmp: Path = None, keep_tmp: bool = False,
                         disable_bar: bool = False) -> Tuple[SequenceFile, bool]:
    """
    Digitalised pangenome genes sequences for hmmsearch

    Args:
        pangenome: Pangenome object with genes
        threads: Number of threads to use
        tmp: Temporary directory to save the gene protein sequences
        keep_tmp: Keep temporary files
        disable_bar: Flag to disable progress bar

    Returns:
        List[pyhmmer.easel.Sequence]: list of digitalised gene family sequences
    """

    write_gene_protein_sequences(pangenome.file, tmp, "all", cpu=threads, keep_tmp=keep_tmp,
                                 tmp=tmp, disable_bar=disable_bar)
    available_memory = psutil.virtual_memory().available
    target_size = os.stat(tmp / "all_protein_genes.fna").st_size
    seq_file = SequenceFile(tmp / "all_protein_genes.fna", digital=True)
    return seq_file, True if target_size < available_memory * 0.1 else False


def digit_family_sequences(pangenome: Pangenome,
                           disable_bar: bool = False) -> tuple[List[DigitalSequence], bool]:
    """
    Digitalised pangenome gene families sequences for HMM alignment

    Args:
        pangenome: Pangenome object with gene families
        disable_bar: Flag to disable progress bar

    Returns:
        List[pyhmmer.easel.Sequence]: list of digitalised gene family sequences
    """
    sequences = []
    logging.getLogger("PANORAMA").info("Begin to digitalized gene families sequences...")
    for family in tqdm(pangenome.gene_families, total=pangenome.number_of_gene_families, unit="gene families",
                       desc="Digitalized gene families sequences", disable=disable_bar):
        bit_name = family.name.encode('UTF-8')
        bit_acc = str(family.ID).encode('UTF-8')
        sequence = family.sequence if family.HMM is None else family.HMM.consensus.upper()
        sequence = sequence.replace('*', '')
        text_seq = TextSequence(name=bit_name, accession=bit_acc, sequence=sequence)
        sequences.append(text_seq.digitize(Alphabet.amino()))
    available_memory = psutil.virtual_memory().available
    return sequences, True if sum(map(sys.getsizeof, sequences)) < available_memory else False


def get_msa(pangenome: Pangenome, tmpdir: Path, threads: int = 1, disable_bar: bool = False) -> pd.DataFrame:
    """
    Get the MSA for each gene families of the pangenome

    Args:
        pangenome: pangenome object with genes
        tmpdir: Temporary directory to save the MSA
        threads: Number of threads to use
        disable_bar: Flag to disable progress bar

    Returns:
        A pandas dataframe containing the MSA for each gene families of the pangenome
    """

    if "translation_table" in pangenome.parameters["annotate"]:
        code = pangenome.parameters["annotate"]["translation_table"]
    else:
        if "translation_table" in pangenome.parameters["cluster"]:
            code = pangenome.parameters["cluster"]["translation_table"]
        else:
            code = "11"
    write_msa_files(pangenome, output=tmpdir, tmpdir=tmpdir, cpu=threads, partition="all", source="protein",
                    use_gene_id=True, translation_table=code, force=True, disable_bar=disable_bar)

    family2msa = {"ID": [], "Path": []}
    for msa_file in Path(tmpdir / 'msa_all_protein').glob(pattern="*.aln"):
        family_name = msa_file.stem
        family2msa["ID"].append(family_name)
        family2msa["Path"].append(msa_file.absolute().as_posix())
    return pd.DataFrame.from_dict(family2msa)


def profile_gf(gf: GeneFamily, msa_path: Path, msa_format: str = "afa", ):
    """
    Compute a profile for a gene family

    Args:
        gf: Gene family to profile
        msa_path: path to file containing msa
        msa_format: format used to write msa

    Raises:
        Exception: Problem to compute profile
    """
    alphabet = Alphabet.amino()
    builder = Builder(alphabet)
    background = Background(alphabet)
    if os.stat(msa_path).st_size == 0:
        logging.getLogger("PANORAMA").debug(f"{msa_path.absolute().as_posix()} is empty, so it's not readable."
                                            f"Pass to next file")
    else:
        try:
            with MSAFile(msa_path.absolute().as_posix(), format=msa_format,
                         digital=True, alphabet=alphabet) as msa_file:
                msa = msa_file.read()
        except Exception as error:
            raise Exception(f"The following error happened while reading file {msa_path} : {error}")
        else:
            try:
                msa.name = gf.name.encode('UTF-8')
                msa.accession = f"PAN{gf.ID}".encode('UTF-8')
                gf._hmm, gf.profile, gf.optimized_profile = builder.build_msa(msa, background)
            except Exception as error:
                raise Exception(f"The following error happened while building HMM from file {msa_path} : {error}")


def profile_gfs(pangenome: Pangenome, msa_df: pd.DataFrame, msa_format: str = "afa",
                threads: int = 1, disable_bar: bool = False):
    """
    Create an HMM profile for each gene families

    Args:
        pangenome: Pangenome containing gene families to profile
        msa_df: Dataframe linking gene families to msa
        msa_format: format used to write msa
        threads: Number of available threads
        disable_bar: Flag to disable progress bar
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.getLogger("PANORAMA").info("Compute gene families HMM and profile")
        with tqdm(total=pangenome.number_of_gene_families, unit='family', disable=disable_bar) as progress:
            futures = []
            msa_df["ID"] = msa_df["ID"].apply(str)
            for family in pangenome.gene_families:
                msa_file = Path(msa_df.loc[msa_df["ID"] == family.name]["Path"].values[0])
                future = executor.submit(profile_gf, family, msa_file, msa_format)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def read_hmms(hmm_db: Path, disable_bar: bool = False) -> Tuple[Dict[str, List[HMM]], pd.DataFrame]:
    """
    Read HMM file to create HMM object

    Args:
        hmm_db: Path to the HMM list file
        disable_bar: Flag to disable the progress bar

    Returns:
        Tuple[Dict[str, List[pyhmmer.plan7.HMM], pd.DataFrame]: A dictionary to identify which cutoff use to align HMM and a pandas dataframe with hmm information

    Raises:
        Exception: Unexpected error occurred while reading HMM
    """
    hmms = defaultdict(list)
    hmm_df = pd.read_csv(hmm_db, delimiter="\t", names=meta_col_names,
                         dtype=meta_dtype, header=0).set_index('accession')
    hmm_df['description'] = hmm_df["description"].fillna('')
    logging.getLogger("PANORAMA").info("Begin to read HMM...")
    for hmm_path in tqdm(map(lambda x: Path(x), hmm_df["path"]), total=hmm_df.shape[0],
                         desc="Reading HMM", unit='HMM', disable=disable_bar):
        end = False
        try:
            hmm_file = HMMFile(hmm_path)
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
                if hmm.cutoffs.gathering_available():
                    hmms["gathering"].append(hmm)
                elif hmm.cutoffs.noise_available():
                    hmms["noise"].append(hmm)
                elif hmm.cutoffs.trusted_available():
                    hmms["trusted"].append(hmm)
                else:
                    hmms["other"].append(hmm)
    return hmms, hmm_df


def assign_hit(hit: Hit, meta: pd.DataFrame) -> Union[Tuple[str, str, str, float, float, float, float, str, str], None]:
    """
    Check if a hit can be assigned to a target

    Args:
        hit: Hit object found by alignment
        meta: metadata information to check target assignment

    Returns:
        informative values to save if alignment is checked
    """
    cog = hit.best_domain.alignment
    hmm_info = meta.loc[cog.hmm_accession.decode('UTF-8')]
    target_coverage = (cog.target_to - cog.target_from) / cog.target_length
    hmm_coverage = (cog.hmm_to - cog.hmm_from) / cog.hmm_length

    check_target_cov = target_coverage >= hmm_info["target_cov_threshold"] or pd.isna(hmm_info["target_cov_threshold"])
    check_hmm_cov = hmm_coverage >= hmm_info["hmm_cov_threshold"] or pd.isna(hmm_info["hmm_cov_threshold"])
    check_score = hit.score >= hmm_info['score_threshold'] or pd.isna(hmm_info['score_threshold'])
    check_e_value = hit.evalue <= hmm_info['eval_threshold'] or pd.isna(hmm_info['eval_threshold'])
    check_ie_value = hit.best_domain.i_evalue <= hmm_info['ieval_threshold'] or pd.isna(hmm_info['ieval_threshold'])

    if check_target_cov and check_hmm_cov and check_score and check_e_value and check_ie_value:
        secondary_name = "" if pd.isna(hmm_info.secondary_name) else hmm_info.secondary_name
        return (hit.name.decode('UTF-8'), cog.hmm_accession.decode('UTF-8'), hmm_info.protein_name, hit.evalue,
                hit.score, hit.bias, hit.best_domain.i_evalue, secondary_name, hmm_info.description)


def annot_with_hmmscan(hmms: Dict[str, List[HMM]], gf_sequences: Union[SequenceFile, List[DigitalSequence]],
                       meta: pd.DataFrame = None, Z: int = 4000, threads: int = 1, tmp: Path = None,
                       disable_bar: bool = False
                       ) -> Tuple[List[Tuple[str, str, str, float, float, float, float, str, str]], List[TopHits]]:
    """
    Compute HMMer alignment between gene families sequences and HMM

    Args:
        hmms: List of HMM classified by bit_cutoffs
        gf_sequences: List of digitalized gene families sequences
        meta: Metadata associate with HMM
        threads: Number of available threads
        tmp: Temporary directory to store pressed HMM database
        disable_bar:  Disable progress bar

    Returns:
         Alignment results
    """

    def hmmscan_callback(seq, _):
        """HMMScan callback function for debugging
        """
        logging.getLogger('PANORAMA').debug(f"Finished annotation with target {seq.name.decode()}")
        bar.update()

    tmp = Path(tempfile.gettempdir()) if tmp is None else tmp
    res = []
    result = namedtuple("Result", res_col_names)
    all_top_hits = []
    hmmpress([hmm for hmm_list in hmms.values() for hmm in hmm_list], tmp / 'hmm_db')
    with HMMFile(tmp / "hmm_db") as hmm_db:
        models = hmm_db.optimized_profiles()
        logging.getLogger("PANORAMA").info("Begin alignment to HMM with HMMScan")
        with tqdm(total=len(gf_sequences), unit="target", desc="Align target to HMM", disable=disable_bar) as bar:
            options = {"Z": Z}
            for cutoff, hmm_list in hmms.items():
                if cutoff != "other":
                    options["bit_cutoffs"] = cutoff
            for top_hits in hmmscan(gf_sequences, models, cpus=threads, callback=hmmscan_callback, **options):
                all_top_hits.append(top_hits)
                for hit in top_hits:
                    assign = assign_hit(hit, meta)
                    if assign is not None:
                        res.append(result(*assign))
    return res, all_top_hits


def annot_with_hmmsearch(hmms: Dict[str, List[HMM]], gf_sequences: SequenceBlock, meta: pd.DataFrame = None,
                         Z: int = 4000, threads: int = 1, disable_bar: bool = False) -> Tuple[List[Tuple[str, str, str, float, float, float, float, str, str]], List[TopHits]]:
    """
    Compute HMMer alignment between gene families sequences and HMM

    Args:
        hmms: List of HMM classified by bit_cutoffs
        gf_sequences: List of digitalized gene families sequences
        meta: Metadata associate with HMM
        threads: Number of available threads
        disable_bar:  Disable progress bar

    Returns:
         Alignment results
    """

    def hmmsearch_callback(hmm, _):
        """HMMSearch callback function for debugging
        """
        logging.getLogger('PANORAMA').debug(f"Finished annotation with HMM {hmm.name.decode()}")
        bar.update()

    res = []
    result = namedtuple("Result", res_col_names)
    logging.getLogger("PANORAMA").info("Begin alignment to HMM with HMMSearch")
    with tqdm(total=sum(len(hmm_list) for hmm_list in hmms.values()), unit="hmm",
              desc="Align target to HMM", disable=disable_bar) as bar:
        all_top_hits = []
        for cutoff, hmm_list in hmms.items():
            options = {"Z": Z}
            if cutoff != "other":
                options["bit_cutoffs"] = cutoff
            for top_hits in hmmsearch(hmm_list, gf_sequences, cpus=threads, callback=hmmsearch_callback, **options):
                all_top_hits.append(top_hits)
                for hit in top_hits:
                    assign = assign_hit(hit, meta)
                    if assign is not None:
                        res.append(result(*assign))

    return res, all_top_hits


def write_top_hits(all_top_hits: List[TopHits], output: Path, source: str, tblout: bool = False,
                   domtblout: bool = False, pfamtblout: bool = False, name: str = 'panorama', mode: str = "fast"):
    """
    Write the pyhmmer hits in the designated format

    Args:
        all_top_hits: List of all pyhmmer hits
        output: Path to the output directory to write hits results
        source: name of the annotation source
        mode: Specify which methods use to align families to HMM
        tblout: Flag to write pyhmmer results in tabular format
        pfamtblout: Flag to write pyhmmer results in Pfam tabular format
        domtblout: Flag to write pyhmmer results with domain tabular format
        name: Name of the pangenome
    """

    header = True
    tbl, domtbl, pfamtbl = None, None, None
    output_path = output / f'{name}' / f'{source}'
    output_path.mkdir(parents=True, exist_ok=True)
    if tblout:
        tbl = open(output_path / f"hmmsearch_{mode}.tbl", "wb")
    if domtblout:
        domtbl = open(output_path / f"hmmsearch_{mode}.domtbl", "wb")
    if pfamtblout:
        pfamtbl = open(output_path / f"hmmsearch_{mode}.pfamtbl", "wb")
    for top_hits in all_top_hits:
        if tblout:
            top_hits.write(tbl, format="targets", header=header)
        if domtblout:
            top_hits.write(domtbl, format="domains", header=header)
        if pfamtblout:
            top_hits.write(pfamtbl, format="pfam", header=header)
        header = False
    if tblout:
        logging.getLogger("PANORAMA").info(f"Per-sequence hits save to file: {tbl.name}")
        tbl.close()
    if domtblout:
        logging.getLogger("PANORAMA").info(f"Per-domain hits save to file: {domtbl.name}")
        domtbl.close()
    if pfamtblout:
        logging.getLogger("PANORAMA").info(f"hits and domains save to file: {pfamtbl.name}")
        pfamtbl.close()


def get_metadata_df(result: List[Tuple[str, str, str, float, float, float, str, str]], mode: str = "fast",
                    gene2family: Dict[str, str] = None) -> pd.DataFrame:
    """
    Parse and refactor the HMM alignment results for filtering and writing step

    Args:
        result: Alignment results with HMM
        mode: Specify which methods use to align families to HMM
        gene2family: Dictionary that link gene to gene_families

    Returns:
        Parsed metadata dataframe
    """
    metadata_df = pd.DataFrame(result).fillna(nan)
    metadata_df.replace(to_replace='-', value=nan, inplace=True)
    metadata_df.replace(to_replace='', value=nan, inplace=True)
    if mode == "sensitive":
        assert gene2family is not None, "Gene and families must be linked in a dictionary"
        gene2family_df = pd.DataFrame.from_dict(gene2family, orient="index").reset_index()
        gene2family_df.columns = ["genes", "families"]
        merged_df = pd.merge(gene2family_df, metadata_df, left_on='genes', right_on='families', how='inner',
                             validate="one_to_many").drop(["genes", "families_y"], axis=1)
        merged_df = merged_df.rename(columns={"families_x": "families"})
        metadata_df = merged_df.drop_duplicates(subset=['families', 'Accession', 'protein_name',
                                                        'e_value', 'score', 'bias'])
        # Keep the best score, e_value, bias for each protein_name by families
        metadata_df = metadata_df.sort_values(by=['score', 'e_value', 'bias'], ascending=[False, True, False])
        group = metadata_df.groupby(["families", "protein_name"])
        metadata_df = group.first().assign(
            secondary_name=group.agg({"secondary_name": lambda x: ",".join(set(x.dropna()))}).replace("", nan)).reset_index()
    return metadata_df


def annot_with_hmm(pangenome: Pangenome, hmms: Dict[str, List[HMM]], meta: pd.DataFrame = None, source: str = "",
                   mode: str = "fast", msa: Path = None, msa_format: str = "afa", tblout: bool = False, Z: int = 4000,
                   domtblout: bool = False, pfamtblout: bool = False, output: Path = None,
                   threads: int = 1, tmp: Path = None, disable_bar: bool = False) -> pd.DataFrame:
    """
    Takes a pangenome and a list of HMMs as input, and returns the best hit for each gene family in the pangenome.

    Args:
        pangenome: Pangenome with gene families
        hmms: Specify the hmm profiles to be used for annotation
        meta: Store the metadata of the hmms
        source: name of the annotation source
        mode: Specify which methods use to align families to HMM
        msa: Path to a msa file listing the gene families and their msa
        msa_format: Specify the format of the msa file
        tblout: Flag to write pyhmmer results in tabular format
        pfamtblout: Flag to write pyhmmer results in Pfam tabular format
        domtblout: Flag to write pyhmmer results with domain tabular format
        output: Path to the output directory to write hits results
        threads: Number of available threads
        tmp: Temporary directory to store pressed HMM database
        disable_bar: bool: Disable the progress bar

    Returns:
        pd.DataFrame: A dataframe with the best results for each pangenome families

    Raises:
        ValueError: if the mode is not recognized. Possible value are fast and profile
        NotImplementedError: If the mode is profile. This option is not implemented yet

    Todo:
        Make the possibility to use the profile with the profile mode
    """
    assert mode in ["sensitive", "fast", "profile"], f"Unrecognized mode: {mode}"
    if (tblout or domtblout or pfamtblout) and output is None:
        raise AssertionError("Output path must be specified to save hits")
    gene2family = None
    if mode == "sensitive":
        sequences, fit_memory = digit_gene_sequences(pangenome, threads, tmp, disable_bar)
        gene2family = {gene.ID: family.name for family in pangenome.gene_families for gene in family.genes}
        if fit_memory:
            logging.getLogger("PANORAMA").debug("Launch pyHMMer-HMMSearch")
            res, all_top_hits = annot_with_hmmsearch(hmms, sequences.read_block(), meta, Z, threads, disable_bar)
        else:
            logging.getLogger("PANORAMA").debug("Launch pyHMMer-HMMScan")
            res, all_top_hits = annot_with_hmmscan(hmms, sequences, meta, Z, threads, tmp, disable_bar)

    else:
        if mode == "profile":
            if msa is not None:
                msa_df = pd.read_csv(msa, sep="\t", names=["ID", "Path"])
            else:
                msa_format = "afa"
                msa_df = get_msa(pangenome, tmp, threads, disable_bar)
            if msa_df.shape[0] != pangenome.number_of_gene_families:
                raise ValueError("The number of msa files does not match the number of gene families")
            profile_gfs(pangenome, msa_df, msa_format, threads, disable_bar)

        # Here either we have profiled our gene family or we will use the referent sequences
        sequences, fit_memory = digit_family_sequences(pangenome, disable_bar=disable_bar)
        if fit_memory:
            logging.getLogger("PANORAMA").debug("Launch pyHMMer-HMMSearch")
            sequence_block = DigitalSequenceBlock(alphabet=Alphabet.amino(), iterable=sequences)
            res, all_top_hits = annot_with_hmmsearch(hmms, sequence_block, meta, Z, threads, disable_bar)
        else:
            logging.getLogger("PANORAMA").debug("Launch pyHMMer-HMMScan")
            res, all_top_hits = annot_with_hmmscan(hmms, sequences, meta, Z, threads, tmp, disable_bar)

    if tblout or domtblout or pfamtblout:
        write_top_hits(all_top_hits, output, source, tblout, domtblout, pfamtblout, pangenome.name, mode)
    metadata_df = get_metadata_df(res, mode, gene2family)
    return metadata_df
