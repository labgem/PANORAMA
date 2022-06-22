#!/usr/bin/env python3
# coding:utf-8

# default libraries
import logging
from pathlib import Path
from typing import List
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile

# installed libraries
import pyhmmer
import pandas as pd

# local libraries
from panorama.pangenomes import Pangenome
from panorama.utils import path_is_dir, path_is_file


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


def hmm_search(hmm_path, gf_sequences: List[pyhmmer.easel.DigitalSequence]) -> List[str, str, float, str]:
    """ Compute a hmmsearch between gene families and one HMM

    :param hmm_path: Path to HMM file
    :param gf_sequences: Digitalized sequences

    :return: result of hmmsearch against gene family sequences
    """
    with pyhmmer.plan7.HMMFile(hmm_path.absolute()) as hmms:
        hmm = next(hmms)
        pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
        top_hits = pipeline.search_hmm(hmm, gf_sequences)
        hit_list = []
        for hit in top_hits:
            hit_list.append([hit.name.decode('UTF-8'), hmm.name.decode('UTF-8'), hit.evalue,
                             hmm.description.decode('UTF-8')])
        return hit_list


def export_to_df(hmm_res: List[List[str, str, float, str]], eval_threshold: float = 0.0001):
    """
    Export result of all HMM comparison to gene family in one dataframe.
    Filter result to keep for each gene family best evalue score greatter than threshold

    :param hmm_res: list of HMM results
    :param eval_threshold: e-value threshold

    :return: Filtered dataframe
    """
    hmm_df = pd.DataFrame(hmm_res, columns=['Gene family name', 'Function name', 'e-value', 'Function description'])
    hmm_filter = hmm_df[hmm_df['e-value'] <= eval_threshold]
    return hmm_filter[hmm_filter['e-value'] == hmm_filter.groupby('Gene family name')['e-value'].transform(min)]


def launch_hmm_search(pangenome: Pangenome, hmm_path: Path, eval_threshold: float = 0.0001,
                      tmpdir: Path = Path(tempfile.gettempdir()), threads: int = 1, disable_bar: bool = False):
    """ Launch hmm search against pangenome gene families and HMM

    :param pangenome: Pangenome with gene families
    :param hmm_path: Path to one file with multiple HMM or a directory with one HMM by file
    :param eval_threshold: e-value threshold to filter results
    :param tmpdir: Path to temporary directory
    :param threads: Number of available threads
    :param disable_bar: allow to disable progress bar

    :return: dataframe with best result for each pangenome families
    """
    logging.getLogger().info("Begin HMM searching")
    logging.getLogger().debug("Digitalized gene families sequences")
    gf_sequences = digit_seq(pangenome)
    tmp_dir = tempfile.TemporaryDirectory(prefix="hmm_panorama", dir=tmpdir)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        logging.getLogger().debug("Intiate pool executor")
        if path_is_file(hmm_path):
            futures = [executor.submit(hmm_search, hmm_path=hmm,
                                       gf_sequences=gf_sequences) for hmm in split_hmm_file(hmm_path, tmp_dir)]
        elif path_is_dir(hmm_path):
            futures = [executor.submit(hmm_search, hmm_path=hmm,
                                       gf_sequences=gf_sequences) for hmm in hmm_path.iterdir()]
        else:
            raise Exception("Unexpected error")
        res = []
        logging.getLogger().debug("Run hmm search")
        for future in tqdm(as_completed(futures), unit="hmm", total=len(futures), disable=disable_bar):
            res += future.result()
    logging.getLogger().debug("Export HMM search results")
    return export_to_df(res, eval_threshold)


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
        launch_hmm_search(pangenome=pangenome, hmm_path=args.hmm, eval_threshold=args.evalue, tmpdir=args.tmpdir,
                          threads=args.threads, disable_bar=args.disable_prog_bar)


    main_parser = argparse.ArgumentParser(
        description="Comparative Pangenomic analyses toolsbox",
        formatter_class=argparse.RawTextHelpFormatter)

    req = main_parser.add_argument_group(title="Required for test")
    req.add_argument("-p", "--pangenome", required=True, type=Path,
                     help="Pangenome.h5 file to test annotation with hmmer")
    req.add_argument("--hmm", required=True, type=Path,
                     help="HMM file to test annotation with hmmer")
    opt = main_parser.add_argument_group(title="Optional argument")
    opt.add_argument("--e_value", required=False, type=float, )
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
