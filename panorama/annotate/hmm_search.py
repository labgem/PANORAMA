#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
import os
import sys
import tempfile
import time
from collections import defaultdict, namedtuple
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas as pd

# installed libraries
import psutil
from numpy import nan
from ppanggolin.formats.writeMSA import write_msa_files
from ppanggolin.formats.writeSequences import write_gene_protein_sequences
from pyhmmer import hmmpress, hmmscan, hmmsearch
from pyhmmer.easel import (
    Alphabet,
    DigitalSequence,
    DigitalSequenceBlock,
    MSAFile,
    SequenceBlock,
    SequenceFile,
    TextSequence,
)
from pyhmmer.plan7 import HMM, Background, Builder, Hit, HMMFile, TopHits
from tqdm import tqdm

from panorama.geneFamily import GeneFamily

# local libraries
from panorama.pangenomes import Pangenome

RESULT_COLUMNS = [
    "families",
    "Accession",
    "protein_name",
    "e_value",
    "score",
    "bias",
    "i_e_value",
    "secondary_names",
    "Description",
]
METADATA_COLUMNS = [
    "name",
    "accession",
    "path",
    "length",
    "protein_name",
    "secondary_name",
    "score_threshold",
    "eval_threshold",
    "ieval_threshold",
    "hmm_cov_threshold",
    "target_cov_threshold",
    "description",
]
METADATA_TYPES = {
    "accession": "string",
    "name": "string",
    "path": "string",
    "length": "int",
    "description": "string",
    "protein_name": "string",
    "secondary_name": "string",
    "score_threshold": "float",
    "eval_threshold": "float",
    "ieval_threshold": "float",
    "hmm_cov_threshold": "float",
    "target_cov_threshold": "float",
}


def digit_gene_sequences(
    pangenome: Pangenome,
    threads: int = 1,
    tmp: Path = None,
    keep_tmp: bool = False,
    disable_bar: bool = False,
) -> Tuple[SequenceFile, bool]:
    """
    Convert all gene sequences from the pangenome into a digital format for HMMER alignment.

    Args:
        pangenome (Pangenome): Pangenome object containing annotated gene data.
        threads (int, optional): Number of threads to use when exporting sequences. Default is 1.
        tmp (Path, optional): Temporary directory where intermediate files will be written.
        If None, uses system temp directory.
        keep_tmp (bool, optional): Whether to keep temporary files after execution. Default is False.
        disable_bar (bool, optional): If True, disables progress bars. Default is False.

    Returns:
        Tuple[SequenceFile, bool]: A tuple containing:
            - The digitalized SequenceFile object for downstream HMMER processing.
            - A boolean indicating whether the size of the sequence file is below 10% of available system memory.
    """

    write_gene_protein_sequences(
        pangenome.file,
        tmp,
        "all",
        cpu=threads,
        keep_tmp=keep_tmp,
        tmp=tmp,
        disable_bar=disable_bar,
    )
    available_memory = psutil.virtual_memory().available
    target_size = os.stat(tmp / "all_protein_genes.faa").st_size
    seq_file = SequenceFile(tmp / "all_protein_genes.faa", digital=True)
    return seq_file, True if target_size < available_memory * 0.1 else False


def digit_family_sequences(pangenome: Pangenome, disable_bar: bool = False) -> Tuple[List[DigitalSequence], bool]:
    """
    Convert each gene family’s consensus or sequence into a digital format for HMM profile creation.

    Args:
        pangenome (Pangenome): Pangenome object containing gene families.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.

    Returns:
        Tuple[List[DigitalSequence], bool]: A list of digitalized sequences for each gene family,
        and a boolean indicating whether the total size is below available system memory.
    """
    sequences = []
    logging.getLogger("PANORAMA").info("Begin to digitalized gene families sequences...")
    for family in tqdm(
        pangenome.gene_families,
        total=pangenome.number_of_gene_families,
        unit="gene families",
        desc="Digitalized gene families sequences",
        disable=disable_bar,
    ):
        bit_name = family.name.encode("UTF-8")
        bit_acc = str(family.ID).encode("UTF-8")
        sequence = family.sequence if family.HMM is None else family.HMM.consensus.upper()
        sequence = sequence.replace("*", "")
        text_seq = TextSequence(name=bit_name, accession=bit_acc, sequence=sequence)
        sequences.append(text_seq.digitize(Alphabet.amino()))
    available_memory = psutil.virtual_memory().available
    return sequences, (True if sum(map(sys.getsizeof, sequences)) < available_memory else False)


def get_msa(pangenome: Pangenome, tmpdir: Path, threads: int = 1, disable_bar: bool = False) -> pd.DataFrame:
    """
    Compute multiple sequence alignments (MSA) for all gene families in the pangenome.

    Args:
        pangenome (Pangenome): Pangenome object containing gene family information.
        tmpdir (Path): Directory to store temporary MSA output files.
        threads (int, optional): Number of threads to use for parallel execution. Default is 1.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.

    Returns:
        pd.DataFrame: A DataFrame mapping each gene family ID to its corresponding MSA file path.
    """

    if "translation_table" in pangenome.parameters["annotate"]:
        code = pangenome.parameters["annotate"]["translation_table"]
    else:
        if "translation_table" in pangenome.parameters["cluster"]:
            code = pangenome.parameters["cluster"]["translation_table"]
        else:
            code = "11"
    write_msa_files(
        pangenome,
        output=tmpdir,
        tmpdir=tmpdir,
        cpu=threads,
        partition="all",
        source="protein",
        use_gene_id=True,
        translation_table=code,
        force=True,
        disable_bar=disable_bar,
    )

    family2msa = {"ID": [], "Path": []}
    for msa_file in Path(tmpdir / "msa_all_protein").glob(pattern="*.aln"):
        family_name = msa_file.stem
        family2msa["ID"].append(family_name)
        family2msa["Path"].append(msa_file.absolute().as_posix())
    return pd.DataFrame.from_dict(family2msa)


def profile_gf(
    gf: GeneFamily,
    msa_path: Path,
    msa_format: str = "afa",
):
    """
    Build an HMM profile for a single gene family using its MSA.

    Args:
        gf (GeneFamily): Gene family object to be profiled.
        msa_path (Path): Path to the MSA file.
        msa_format (str, optional): Format of the MSA file (e.g., "afa"). Default is "afa".

    Raises:
        Exception: If the MSA file is unreadable or if building the HMM profile fails.
    """
    alphabet = Alphabet.amino()
    builder = Builder(alphabet)
    background = Background(alphabet)
    if os.stat(msa_path).st_size == 0:
        logging.getLogger("PANORAMA").debug(
            f"{msa_path.absolute().as_posix()} is empty, so it's not readable.Pass to next file"
        )
    else:
        try:
            with MSAFile(
                msa_path.absolute().as_posix(),
                format=msa_format,
                digital=True,
                alphabet=alphabet,
            ) as msa_file:
                msa = msa_file.read()
        except Exception as error:
            raise Exception(f"The following error happened while reading file {msa_path}") from error
        else:
            try:
                msa.name = gf.name.encode("UTF-8")
                msa.accession = f"PAN{gf.ID}".encode("UTF-8")
                gf._hmm, gf.profile, gf.optimized_profile = builder.build_msa(msa, background)
            except Exception as error:
                raise Exception(f"The following error happened while building HMM from file {msa_path}") from error


def profile_gfs(
    pangenome: Pangenome,
    msa_df: pd.DataFrame,
    msa_format: str = "afa",
    threads: int = 1,
    disable_bar: bool = False,
):
    """
    Generate HMM profiles for all gene families in the pangenome.

    Args:
        pangenome (Pangenome): Pangenome object containing gene families.
        msa_df (pd.DataFrame): DataFrame mapping gene family IDs to MSA file paths.
        msa_format (str, optional): Format used to read MSA files. Default is "afa".
        threads (int, optional): Number of threads for parallel processing. Default is 1.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.getLogger("PANORAMA").info("Compute gene families HMM and profile")
        with tqdm(total=pangenome.number_of_gene_families, unit="family", disable=disable_bar) as progress:
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
    Read a set of HMM files and categorize them based on available cutoffs.

    Args:
        hmm_db (Path): Path to the tab-delimited file listing HMM metadata.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.

    Returns:
        Tuple[Dict[str, List[HMM]], pd.DataFrame]: A dictionary categorizing HMMs by cutoff type
        (gathering, trusted, noise, or None), and a DataFrame with metadata.

    Raises:
        Exception: If reading an HMM file fails unexpectedly.
    """
    hmms = defaultdict(list)
    hmm_df = pd.read_csv(hmm_db, delimiter="\t", names=METADATA_COLUMNS, dtype=METADATA_TYPES, header=0).set_index(
        "accession"
    )
    hmm_df["description"] = hmm_df["description"].fillna("")
    logging.getLogger("PANORAMA").info("Begin to read HMM...")
    for hmm_path in tqdm(
        map(lambda x: Path(x), hmm_df["path"]),
        total=hmm_df.shape[0],
        desc="Reading HMM",
        unit="HMM",
        disable=disable_bar,
    ):
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
                raise Exception(f"Unexpected error on HMM file {hmm_path}") from error
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
    Evaluate whether a hit from HMMER alignment meets filtering criteria.

    Args:
        hit (Hit): HMMER hit object representing a match between a query and a model.
        meta (pd.DataFrame): DataFrame containing annotation thresholds and metadata.

    Returns:
        Union[Tuple[str, str, str, float, float, float, float, str, str], None]: A tuple of annotation
        information if the hit passes filters, otherwise None.
    """
    cog = hit.best_domain.alignment
    hmm_info = meta.loc[cog.hmm_accession.decode("UTF-8")]
    target_coverage = (cog.target_to - cog.target_from) / cog.target_length
    hmm_coverage = (cog.hmm_to - cog.hmm_from) / cog.hmm_length

    check_target_cov = target_coverage >= hmm_info["target_cov_threshold"] or pd.isna(hmm_info["target_cov_threshold"])
    check_hmm_cov = hmm_coverage >= hmm_info["hmm_cov_threshold"] or pd.isna(hmm_info["hmm_cov_threshold"])
    check_score = hit.score >= hmm_info["score_threshold"] or pd.isna(hmm_info["score_threshold"])
    check_e_value = hit.evalue <= hmm_info["eval_threshold"] or pd.isna(hmm_info["eval_threshold"])
    check_ie_value = hit.best_domain.i_evalue <= hmm_info["ieval_threshold"] or pd.isna(hmm_info["ieval_threshold"])

    if check_target_cov and check_hmm_cov and check_score and check_e_value and check_ie_value:
        secondary_name = "" if pd.isna(hmm_info.secondary_name) else hmm_info.secondary_name
        return (
            hit.name.decode("UTF-8"),
            cog.hmm_accession.decode("UTF-8"),
            hmm_info.protein_name,
            hit.evalue,
            hit.score,
            hit.bias,
            hit.best_domain.i_evalue,
            secondary_name,
            hmm_info.description,
        )
    return None


def annot_with_hmmscan(
    hmms: Dict[str, List[HMM]],
    gf_sequences: Union[SequenceFile, List[DigitalSequence]],
    meta: pd.DataFrame = None,
    Z: int = 4000,
    threads: int = 1,
    tmp: Path = None,
    disable_bar: bool = False,
) -> Tuple[List[Tuple[str, str, str, float, float, float, float, str, str]], List[TopHits]]:
    """
    Annotate sequences by scanning them against HMM profiles using HMMER's `hmmscan`.

    Args:
        hmms (Dict[str, List[HMM]]): Dictionary of HMMs grouped by cutoff type.
        gf_sequences (Union[SequenceFile, List[DigitalSequence]]): Digital sequences to annotate.
        meta (pd.DataFrame, optional): Optional metadata for evaluating hit criteria.
        Z (int, optional): Effective number of database comparisons. Default is 4000.
        threads (int, optional): Number of threads to use. Default is 1.
        tmp (Path, optional): Temporary directory for intermediate files. Default is system temp dir.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.

    Returns:
        Tuple[List[Tuple[str, str, str, float, float, float, float, str, str]], List[TopHits]]:
            - Annotation results that pass filters.
            - All raw TopHits results from HMMER.
    """

    def hmmscan_callback(seq, _):
        """
        Callback function triggered after each sequence is processed by `hmmscan`.

        Primarily used for debugging and tracking progress. Logs the completion
        of annotation for a given target sequence and updates the progress bar.

        Args:
            seq (Sequence): The annotated sequence object returned by `hmmscan`.
            _ (Any): Placeholder for an unused second argument (typically error object or result metadata).
        """
        logging.getLogger("PANORAMA").debug(f"Finished annotation with target {seq.name.decode()}")
        bar.update()

    tmp = Path(tempfile.gettempdir()) if tmp is None else tmp
    res = []
    result = namedtuple("Result", RESULT_COLUMNS)
    all_top_hits = []
    hmmpress([hmm for hmm_list in hmms.values() for hmm in hmm_list], tmp / "hmm_db")
    with HMMFile(tmp / "hmm_db") as hmm_db:
        models = hmm_db.optimized_profiles()
        logging.getLogger("PANORAMA").info("Begin alignment to HMM with HMMScan")
        with tqdm(
            total=len(gf_sequences),
            unit="target",
            desc="Align target to HMM",
            disable=disable_bar,
        ) as bar:
            options = {"Z": Z}
            for cutoff in hmms.keys():
                if cutoff != "other":
                    options["bit_cutoffs"] = cutoff
            for top_hits in hmmscan(gf_sequences, models, cpus=threads, callback=hmmscan_callback, **options):
                all_top_hits.append(top_hits)
                for hit in top_hits:
                    assign = assign_hit(hit, meta)
                    if assign is not None:
                        res.append(result(*assign))
    return res, all_top_hits


def annot_with_hmmsearch(
    hmms: Dict[str, List[HMM]],
    gf_sequences: SequenceBlock,
    meta: pd.DataFrame = None,
    Z: int = 4000,
    threads: int = 1,
    disable_bar: bool = False,
) -> Tuple[List[Tuple[str, str, str, float, float, float, float, str, str]], List[TopHits]]:
    """
    Annotate HMM profiles by searching them against a block of target sequences using HMMER's `hmmsearch`.

    Args:
        hmms (Dict[str, List[HMM]]): Dictionary of HMMs grouped by cutoff type.
        gf_sequences (SequenceBlock): Digital sequence block representing target sequences.
        meta (pd.DataFrame, optional): Optional metadata for evaluating hit criteria.
        Z (int, optional): Effective number of database comparisons. Default is 4000.
        threads (int, optional): Number of threads to use. Default is 1.
        disable_bar (bool, optional): If True, disables the progress bar. Default is False.

    Returns:
        Tuple[List[Tuple[str, str, str, float, float, float, float, str, str]], List[TopHits]]:
            - Filtered annotation results.
            - All TopHits results returned by HMMER.
    """

    def hmmsearch_callback(hmm, _):
        """
        Callback function triggered after each HMM is processed by `hmmsearch`.

        Used for debugging and progress tracking. Logs the completion of
        annotation for a given HMM profile and updates the progress bar.

        Args:
            hmm (HMMProfile): The HMM profile object that was used in the search.
            _ (Any): Placeholder for an unused second argument (typically error object or result metadata).
        """
        logging.getLogger("PANORAMA").debug(f"Finished annotation with HMM {hmm.name.decode()}")
        bar.update()

    res = []
    result = namedtuple("Result", RESULT_COLUMNS)
    logging.getLogger("PANORAMA").info("Begin alignment to HMM with HMMSearch")
    with tqdm(
        total=sum(len(hmm_list) for hmm_list in hmms.values()),
        unit="hmm",
        desc="Align target to HMM",
        disable=disable_bar,
    ) as bar:
        all_top_hits = []
        for cutoff, hmm_list in hmms.items():
            options = {"Z": Z}
            if cutoff != "other":
                options["bit_cutoffs"] = cutoff
            for top_hits in hmmsearch(
                hmm_list,
                gf_sequences,
                cpus=threads,
                callback=hmmsearch_callback,
                **options,
            ):
                all_top_hits.append(top_hits)
                for hit in top_hits:
                    assign = assign_hit(hit, meta)
                    if assign is not None:
                        res.append(result(*assign))

    return res, all_top_hits


def write_top_hits(
    all_top_hits: List[TopHits],
    output: Path,
    source: str,
    tblout: bool = False,
    domtblout: bool = False,
    pfamtblout: bool = False,
    name: str = "panorama",
    mode: str = "fast",
):
    """
    Write pyhmmer search hits to file in various tabular formats.

    Depending on the flags provided, writes per-sequence (`tbl`), per-domain (`domtbl`),
    and/or Pfam-style (`pfamtbl`) formatted results.

    Args:
        all_top_hits (List[TopHits]): List of pyhmmer TopHits objects.
        output (Path): Directory where output files will be written.
        source (str): Name of the annotation source (used in subfolder naming).
        tblout (bool): If True, write per-sequence hits (`*.tbl`).
        domtblout (bool): If True, write per-domain hits (`*.domtbl`).
        pfamtblout (bool): If True, write hits in Pfam format (`*.pfamtbl`).
        name (str): Name of the pangenome (used for folder structure).
        mode (str): Alignment mode used for the annotation (e.g., "fast", "sensitive").
    """

    header = True
    tbl, domtbl, pfamtbl = None, None, None
    output_path = output / f"{name}" / f"{source}"
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


def get_metadata_df(
    result: List[Tuple[str, str, str, float, float, float, float, str, str]],
    mode: str = "fast",
    gene2family: Dict[str, str] = None,
) -> pd.DataFrame:
    """
    Refactor HMM alignment results into a structured metadata DataFrame.

    Handles basic cleaning and optionally joins gene-family metadata in
    "sensitive" mode to allow grouping by family.

    Args:
        result (List[Tuple]):
            List of raw alignment results from HMM search.
        mode (str):
            Annotation mode used ("fast" or "sensitive").
        gene2family (Dict[str, str], optional):
            Required for "sensitive" mode. Maps gene IDs to family names.

    Returns:
        pd.DataFrame:
            Cleaned and optionally merged metadata DataFrame.
    """
    metadata_df = pd.DataFrame(result).fillna(nan)
    metadata_df = metadata_df.where(metadata_df != "-", nan)
    metadata_df = metadata_df.where(metadata_df != "", nan)
    if mode == "sensitive":
        assert gene2family is not None, "Gene and families must be linked in a dictionary"
        gene2family_df = pd.DataFrame.from_dict(gene2family, orient="index").reset_index()
        gene2family_df.columns = ["genes", "families"]
        merged_df = pd.merge(
            gene2family_df,
            metadata_df,
            left_on="genes",
            right_on="families",
            how="inner",
            validate="one_to_many",
        ).drop(["genes", "families_y"], axis=1)
        merged_df = merged_df.rename(columns={"families_x": "families"})
        metadata_df = merged_df.drop_duplicates(
            subset=["families", "Accession", "protein_name", "e_value", "score", "bias"]
        )
        # Keep the best score, e_value, bias for each protein_name by families
        metadata_df = metadata_df.sort_values(by=["score", "e_value", "bias"], ascending=[False, True, False])
        group = metadata_df.groupby(["families", "protein_name"])
        metadata_df = (
            group.first()
            .assign(secondary_name=group.agg({"secondary_names": lambda x: ",".join(set(x.dropna()))}).replace("", nan))
            .reset_index()
        )
    return metadata_df


def annot_with_hmm(
    pangenome: Pangenome,
    hmms: Dict[str, List[HMM]],
    meta: pd.DataFrame = None,
    source: str = "",
    mode: str = "fast",
    msa: Path = None,
    msa_format: str = "afa",
    tblout: bool = False,
    Z: int = 4000,
    domtblout: bool = False,
    pfamtblout: bool = False,
    output: Path = None,
    threads: int = 1,
    tmp: Path = None,
    disable_bar: bool = False,
) -> pd.DataFrame:
    """
    Annotate a pangenome using a collection of HMM profiles and return best hits per family.

    Supports three annotation modes:
    - "fast": Uses referent sequences or representative sequences.
    - "sensitive": Uses all genes in gene families with HMMScan/HMMSearch.
    - "profile": Not yet implemented; intended to use family alignments.

    Args:
        pangenome (Pangenome): Pangenome object with gene families.
        hmms (Dict[str, List[HMM]]): Dictionary of annotation source → list of HMMs.
        meta (pd.DataFrame, optional): Optional metadata for the HMMs.
        source (str): Name of the annotation source.
        mode (str): Annotation mode to use ("fast", "sensitive", "profile").
        msa (Path, optional): Path to a file listing MSAs (only used in "profile" mode).
        msa_format (str): Format of MSAs if provided (default: "afa").
        tblout (bool): If True, write per-sequence hits.
        pfamtblout (bool): If True, write Pfam-format hits.
        Z (int): Effective number of comparisons (default: 4000).
        domtblout (bool): If True, write per-domain hits.
        output (Path, optional): Directory to write annotation results.
        threads (int): Number of threads for parallel processing.
        tmp (Path, optional): Temporary directory for intermediate files.
        disable_bar (bool): If True, disable progress bars.

    Returns:
        pd.DataFrame: DataFrame containing best hit per gene or per family.

    Raises:
        ValueError: If the number of MSAs doesn't match the number of gene families.
        AssertionError: If output is required but not provided.
        NotImplementedError: If "profile" mode is selected (not yet available).

    TODO:
        Implement "profile" mode with family alignments.
    """
    assert mode in ["sensitive", "fast", "profile"], f"Unrecognized mode: {mode}"
    if (tblout or domtblout or pfamtblout) and output is None:
        raise AssertionError("Output path must be specified to save hits")
    gene2family = None
    t0 = time.time()
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
        write_top_hits(
            all_top_hits,
            output,
            source,
            tblout,
            domtblout,
            pfamtblout,
            pangenome.name,
            mode,
        )
    metadata_df = get_metadata_df(res, mode, gene2family)
    logging.getLogger("PANORAMA").info(
        f"Annotation with HMM done for {pangenome.name} in {time.time() - t0:2f} seconds"
    )
    return metadata_df
