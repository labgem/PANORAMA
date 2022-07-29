#!/usr/bin/env python3
# coding:utf-8

# default libraries
import collections
import logging
from pathlib import Path
from typing import List, Union, Generator
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile

# installed libraries
import pyhmmer
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenome
from panorama.utils import path_is_dir, path_is_file

res_col_names = ['Gene_family', 'Annotation', 'Accession', 'e_value', 'score', 'overlap', 'Description']
cutoffs_col_names = ["accession", "name", "eval", "hmm_cov", "target_cov", "description"]


def digit_seq(pangenome: Pangenome) -> List[pyhmmer.easel.DigitalSequence]:
    """Digitalised pangenome gene families sequences for hmmsearch

    :param pangenome: Pangenome object with gene families

    :return: list of digitalised gene family sequences
    """
    digit_gf_seq = []
    for family in pangenome.gene_families:
        bit_name = family.name.encode('UTF-8')
        seq = pyhmmer.easel.TextSequence(name=bit_name, sequence=family.sequence)
        digit_gf_seq.append(seq.digitize(pyhmmer.easel.Alphabet.amino()))
    return digit_gf_seq


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
                file = open(hmm_list[-1].as_posix(), 'w')
                file.write("\n".join(lines))
                counter += 1
                lines = [line]
            else:
                lines.append(line)
        hmm_list.append(Path(f"{tmpdir.name}/HMM_{counter}.hmm"))
        file = open(hmm_list[-1].as_posix(), 'w')
        file.write("\n".join(lines))
    return hmm_list


def hmm_search_hmmer(hmm_list: List[pyhmmer.plan7.HMM], gf_sequences: List[pyhmmer.easel.DigitalSequence],
                     hmm_meta: pd.DataFrame = None, threads: int = 1, disable_bar: bool = False):
    res = []
    result = collections.namedtuple("Result", res_col_names)
    bar = tqdm(range(len(hmm_list)), unit="hmm", disable=disable_bar)
    for top_hits in pyhmmer.hmmsearch(hmm_list, gf_sequences, cpus=threads, bit_cutoffs='trusted'):
        for hit in top_hits:
            cog = hit.best_domain.alignment
            if hmm_meta is None:
                res.append(result(hit.name.decode('UTF-8'), cog.hmm_name.decode('UTF-8'),
                                  cog.hmm_accession.decode('UTF-8'), hit.evalue, hit.score, None, None))
            else:
                hmm_info = hmm_meta.loc[cog.hmm_accession.decode('UTF-8')]
                res.append(result(hit.name.decode('UTF-8'), hmm_info['name'], cog.hmm_accession.decode('UTF-8'),
                                  hit.evalue, hit.score, None, hmm_info.description))
        bar.update()
    bar.close()
    return res


def annot_with_hmmsearch(hmm_list: List[pyhmmer.plan7.HMM], gf_sequences: List[pyhmmer.easel.DigitalSequence],
                         hmm_meta: pd.DataFrame = None, threads: int = 1, disable_bar: bool = False):
    results = hmm_search_hmmer(hmm_list, gf_sequences, hmm_meta, threads, disable_bar)
    best_results = {}
    keep_query = set()
    keep_size = 5
    for result in results:
        if result.Gene_family in best_results:
            # Search if Annotation already exist for families
            annot_in = False
            for index, res in enumerate(best_results[result.Gene_family]):
                if result.Annotation == res.Annotation:
                    if result.score > res.score:
                        best_results[result.Gene_family][index] = result
                        annot_in = True
                        break
            if not annot_in:
                if len(best_results[result.Gene_family]) < keep_size:
                    best_results[result.Gene_family].append(result)
                else:
                    previous_bitscore = best_results[result.Gene_family][-1].score
                    if result.score > previous_bitscore:
                        best_results[result.Gene_family][-1] = result
                    elif result.score == previous_bitscore:
                        if result.Annotation not in [res.Annotation for res in best_results[result.Gene_family]]:
                            fam_rm = best_results[result.Gene_family].pop(-1)
                            if len(best_results[result.Gene_family]) == 0:
                                keep_query.remove(fam_rm)
                best_results[result.Gene_family] = sorted(best_results[result.Gene_family],
                                                          key=lambda x: x.score, reverse=True)
        else:
            best_results[result.Gene_family] = [result]
            keep_query.add(result.Gene_family)
    return sorted([res for gf, res_list in best_results.items() for res in res_list if gf in keep_query],
                  key=lambda x: x.Gene_family)


def hmm_search_plan7(hmm: pyhmmer.plan7.HMM,
                     gf_sequences: List[pyhmmer.easel.DigitalSequence]) -> List[List[Union[str, str, int, str]]]:
    """ Compute a hmmsearch between gene families and one HMM

    :param hmm: HMM with Plan7 data model
    :param gf_sequences: Digitalized sequences

    :return: result of hmmsearch against gene family sequences
    """
    pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet, bit_cutoffs='trusted')
    top_hits = pipeline.search_hmm(hmm, gf_sequences)
    hit_list = []
    for hit in top_hits:
        hit_list.append([hit.name.decode('UTF-8'), hmm.name.decode('UTF-8'), hmm.accession.decode('UTF-8'),
                         hit.evalue, hit.score, None, hmm.description.decode('UTF-8')])
    return hit_list


def annot_with_plan7(hmm_list: List[pyhmmer.plan7.HMM], gf_sequences: List[pyhmmer.easel.DigitalSequence],
                     threads: int = 1, disable_bar: bool = False):
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.getLogger().debug("Intiate pool executor")
        futures = [executor.submit(hmm_search_plan7, hmm=hmm, gf_sequences=gf_sequences) for hmm in hmm_list]
        res = []
        logging.getLogger().debug("Run hmm search")
        for future in tqdm(as_completed(futures), unit="hmm", total=len(futures), disable=disable_bar):
            res += future.result()
        return res


def get_hmm_iter_eval(hmm_dir: Generator, eval: float = 0.0001):
    hmms = []
    for hmm_file in hmm_dir:
        hmm = next(pyhmmer.plan7.HMMFile(hmm_file))
        hmm.cutoffs.trusted = (eval, eval)
        hmms.append(hmm)
    return hmms


def get_hmm_iter_cutoffs(hmm_dir: Generator, cutoffs: pd.DataFrame):
    hmms = []
    for hmm_file in hmm_dir:
        hmm = next(pyhmmer.plan7.HMMFile(hmm_file))
        meta_info = cutoffs.loc[hmm.accession.decode('UTF-8')]
        hmm.cutoffs.trusted = (meta_info.eval, meta_info.eval)
        hmm.description = meta_info.description.encode('UTF-8')
        hmms.append(hmm)
    return hmms


def annot_with_hmm(pangenome: Pangenome, hmm_path: Path, method: str = 'hmmsearch', eval: float = None,
                   cutoffs_file: Path = None, tmpdir: Path = Path(tempfile.gettempdir()), threads: int = 1,
                   disable_bar: bool = False):
    """ Launch hmm search against pangenome gene families and HMM

    :param pangenome: Pangenome with gene families
    :param hmm_path: Path to one file with multiple HMM or a directory with one HMM by file
    :param method: Methods used to search HMM
    :param eval:
    :param cutoffs_file:
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param disable_bar: allow to disable progress bar

    :return: dataframe with best result for each pangenome families
    """
    logging.getLogger().info("Begin HMM searching")
    logging.getLogger().debug("Digitalized gene families sequences")
    gf_sequences = digit_seq(pangenome)
    tmp_dir = tempfile.TemporaryDirectory(prefix="hmm_panorama", dir=tmpdir)
    cutoffs_df = None
    # Get list of HMM with Plan7 data model
    if path_is_file(hmm_path):
        split_hmm_file(hmm_path, tmp_dir)
        if eval is not None:
            hmms = get_hmm_iter_eval(hmm_dir=Path(tmp_dir.name).iterdir(), eval=eval)
        elif cutoffs_file is not None:
            cutoffs_df = pd.read_csv(cutoffs_file, delimiter="\t", names=cutoffs_col_names).set_index('accession')
            cutoffs_df['description'] = cutoffs_df["description"].fillna('unknown')
            hmms = get_hmm_iter_cutoffs(hmm_dir=Path(tmp_dir.name).iterdir(), cutoffs=cutoffs_df)
        else:
            raise Exception("You didn't provide e-value or cutoffs file")
    elif path_is_dir(hmm_path):
        if eval is not None:
            hmms = get_hmm_iter_eval(hmm_dir=hmm_path.iterdir(), eval=eval)
        elif cutoffs_file is not None:
            cutoffs_df = pd.read_csv(cutoffs_file, delimiter="\t", names=cutoffs_col_names).set_index('accession')
            cutoffs_df['description'] = cutoffs_df["description"].fillna('unknown')
            hmms = get_hmm_iter_cutoffs(hmm_dir=hmm_path.iterdir(), cutoffs=cutoffs_df)
        else:
            raise Exception("You didn't provide e-value or cutoffs file")
    else:
        raise Exception("Unexpected error")
    # Choice of one method
    if method == 'hmmsearch':
        res = annot_with_hmmsearch(hmms, gf_sequences, cutoffs_df, threads, disable_bar)
    elif method == 'plan7':
        res = annot_with_plan7(hmms, gf_sequences, threads, disable_bar)
    else:
        raise Exception("Methods to search HMM is plan7 or hmmsearch. Please choose one of them.")
    return pd.DataFrame(res, columns=res_col_names)


if __name__ == "__main__":
    # default libraries
    import argparse
    # installed libraries
    from ppanggolin.formats.readBinaries import check_pangenome_info
    # local libraries
    from panorama.utils import check_log, set_verbosity_level


    def test_hmm_search(args: argparse.Namespace):
        pangenome = Pangenome(name="pangenome")
        pangenome.add_file(pangenome_file=args.pangenome)
        check_pangenome_info(pangenome, need_families=True)
        res_df = annot_with_hmm(pangenome=pangenome, hmm_path=args.hmm, method=args.method, tmpdir=args.tmpdir,
                                eval=args.e_value, cutoffs_file=args.cutoffs, threads=args.threads,
                                disable_bar=args.disable_prog_bar)
        print(res_df)
        res_df.to_csv(path_or_buf="hmm_res.tsv", sep="\t")


    main_parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
        formatter_class=argparse.RawTextHelpFormatter)

    req = main_parser.add_argument_group(title="Required for test")
    req.add_argument("-p", "--pangenome", required=True, type=Path,
                     help="Pangenome.h5 file to test annotation with hmmer")
    req.add_argument("--hmm", required=True, type=Path,
                     help="HMM file to test annotation with hmmer")
    opt = main_parser.add_argument_group(title="Optional argument")
    opt.add_argument("--method", required=False, type=str, choices=['hmmsearch', 'plan7'], default='hmmsearch',
                     help="Two methods are available and defined in pyhmmer. Plan7 or hmmsearch")
    exclusive = opt.add_mutually_exclusive_group()
    exclusive.add_argument('--e_value', required=False, type=float, default=None)
    exclusive.add_argument('--cutoffs', required=False, type=Path, default=None)
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
    test_hmm_search(main_args)
