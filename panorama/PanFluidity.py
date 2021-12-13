#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from multiprocessing import get_context
from time import time

# installed libraries
import ppanggolin.formats
from gmpy2 import popcount
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo
from itertools import combinations
from tqdm import tqdm
import tables


def genomes_fluidity(pangenome: Pangenome, disable_bar: bool = False) -> float:
    """ Compute the genomes fluidity from the pangenome

    :param pangenome: pangenome which will be used to compute the genomes fluidity
    :param disable_bar: Disable the progress bar

    :return: Genomes fluidity value from the pangenome
    """
    # check statuses and load info
    logging.getLogger().info("Check informations in pangenome")
    checkPangenomeInfo(pangenome, needAnnotations=True, needFamilies=True, disable_bar=disable_bar)
    logging.getLogger().debug("Compute binaries sequences corresponding to presence / absence of families in organisms")
    pangenome.compute_org_bitarrays()  # Compute binaries corresponding to presence / absence of families in organisms
    t0 = time()
    g_sum = 0
    t4 = time()
    logging.getLogger().debug("Get number of families in each organisms")
    org2_nb_fam = nb_fam_per_org(pangenome, disable_bar)
    print(time() - t4)
    logging.getLogger().info("Compute rate of unique family for each genome combination")
    for c_organisms in tqdm(list(combinations(pangenome.organisms, 2)), unit="combination", disable=disable_bar):
        tot_fam = org2_nb_fam.get(c_organisms[0].name) + org2_nb_fam.get(c_organisms[1].name)
        common_fam = popcount(c_organisms[0].bitarray & c_organisms[1].bitarray) - 1
        g_sum += (tot_fam - 2 * common_fam) / tot_fam
    print("Temps de calcul : ", time() - t0)
    return (2 / (pangenome.number_of_organisms() * (pangenome.number_of_organisms() - 1))) * g_sum


def nb_fam_per_org(pangenome: Pangenome, disable_bar: bool = False) -> dict:
    """
    Create a dictionary with for each organism the number of gene families

    :param pangenome: Pangenome which contain the organisms and gene families
    :param disable_bar: Disable the progress bar
    :return: Dictionary with organisms as key and number of families as value
    """
    org2_nb_fam = dict()
    for org in tqdm(pangenome.organisms, unit='organism', disable=disable_bar):
        org2_nb_fam[org.name] = popcount(org.bitarray)
    return org2_nb_fam


def write_fluidity(pangenome: Pangenome, g_fluidity: float):
    """Write or update the genomes fluidity value in pangenome

    :param pangenome: pangenome where it should be written the genomes fluidity
    :param g_fluidity: genomes_fluidity computed for the pangenome
    """

    h5f = tables.open_file(pangenome.file, "a")
    if "/info" not in h5f:
        logging.getLogger().info("Your pangenome is without information. Information will be write to be relevant with "
                                 "genomes fluidity")
        ppanggolin.formats.writeInfo(pangenome, h5f)

    info_group = h5f.root.info
    info_group._v_attrs.fluidity = g_fluidity
    h5f.close()
    logging.getLogger().info("Done writing the genome fluidity in pangenome")


# TODO Function to normalize genome fluidity
# def genomes_fluidity_norm(p_pangenome: Pangenome, disable_bar: bool = False) -> float:

# TODO Function to compute mash distance between genome

# TODO Function to export results

def mp_genome_fluidity(pangenome_file: str) -> str:
    """ Compute in multiprocessing the genome fluidity of multiple pangenome

    :param pangenome_file: one of the pangenome in the list of pangenome

    :return: Name of the computed pangenome
    """
    # TODO Allow to show progress bar of each process
    pangenome = Pangenome()
    pangenome.addFile(pangenome_file)
    g = genomes_fluidity(pangenome=pangenome, disable_bar=True)
    write_fluidity(pangenome, g)
    return pangenome_file


def launch(args: argparse.Namespace):
    """Allow to launch the panfluidity script from panorama command

    :param args: list of arguments
    """
    if len(args.pangenomes) == 0:
        raise Exception("Not one pangenome was found. Please check path to pangenome files.")
    if args.cpu == 1:
        for pangenome_file in args.pangenomes:
            pangenome = Pangenome()
            pangenome.addFile(pangenome_file)
            g = genomes_fluidity(pangenome=pangenome, disable_bar=args.disable_prog_bar)
            write_fluidity(pangenome, g)
    else:
        logging.getLogger().info("Begin of the genomes fluidity computation")
        with get_context('fork').Pool(args.cpu) as p:
            for pangenome_file in tqdm(p.imap_unordered(mp_genome_fluidity, args.pangenomes),
                                       unit='pangenome', total=len(args.pangenomes), disable=args.disable_prog_bar):
                logging.getLogger().debug(f"{pangenome_file} Done")
        logging.getLogger().info("All the genomes fluidity were computed")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Parser arguments specific to panfluidity command
    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("panfluidity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required")
    required.add_argument('-p', '--pangenomes', required=True, type=str, nargs='+',
                          help="A list of pangenome .h5 files")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")
    optional = parser.add_argument_group(title="Optional arguments",
                                         description="All of the following arguments are optional and"
                                                     " with a default value")
    optional.add_argument("-c", "--cpu", required=False, default=1, type=int, help="Number of available cpus")

    return parser
