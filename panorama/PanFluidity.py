#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
# import statistics
from time import time

# installed libraries
import ppanggolin.formats
from gmpy2 import popcount
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo
from itertools import combinations
from tqdm import tqdm
import tables


def genomes_fluidity(p_pangenome: Pangenome, disable_bar: bool = False) -> float:
    """ Compute the genomes fluidity from the pangenome

    :param p_pangenome: pangenome which will be used to compute the genomes fluidity
    :param disable_bar: Disable the progress bar

    :return: Genomes fluidity value from the pangenome
    """
    # check statuses and load info
    logging.getLogger().info("Check informations in pangenome")
    checkPangenomeInfo(p_pangenome, needAnnotations=True, needFamilies=True, disable_bar=disable_bar)
    logging.getLogger().debug("Compute binaries sequences corresponding to presence / absence of families in organisms")
    p_pangenome.compute_org_bitarrays()  # Compute binaries corresponding to presence / absence of families in organisms
    t0 = time()
    g_sum = 0
    # t1, t2, t3 = ([], [], [])
    t4 = time()
    logging.getLogger().debug("Get number of families in each organisms")
    org2_nb_fam = nb_fam_per_org(p_pangenome, disable_bar)
    print(time() - t4)
    logging.getLogger().info("Compute rate of unique family for each genome combination")
    for c_organisms in tqdm(list(combinations(p_pangenome.organisms, 2)), unit="combination", disable=disable_bar):
        # t1.append(time())
        tot_fam = org2_nb_fam.get(c_organisms[0].name) + org2_nb_fam.get(c_organisms[1].name)
        # t2.append(time())
        common_fam = popcount(c_organisms[0].bitarray & c_organisms[1].bitarray) - 1
        # t3.append(time())
        g_sum += (tot_fam - 2 * common_fam) / tot_fam
    # dif1 = [b - a for a, b in zip(t1, t2)]
    # dif2 = [b - a for a, b in zip(t2, t3)]

    # print(statistics.mean(dif1), statistics.mean(dif2))
    print("Temps de calcul : ", time() - t0)
    return (2 / (p_pangenome.number_of_organisms() * (p_pangenome.number_of_organisms() - 1))) * g_sum


def nb_fam_per_org(p_pangenome: Pangenome, disable_bar: bool = False) -> dict:
    """
    Create a dictionary with for each organism the number of gene families

    :param p_pangenome: Pangenome which contain the organisms and gene families
    :param disable_bar: Disable the progress bar
    :return: Dictionary with organisms as key and number of families as value
    """
    org2_nb_fam = dict()
    for org in tqdm(p_pangenome.organisms, unit='organism', disable=disable_bar):
        org2_nb_fam[org.name] = popcount(org.bitarray)
    return org2_nb_fam


def write_fluidity(pangenome: Pangenome, g_fluidity: float):
    """Write or update the genomes fluidity value in pangenome
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

# TODO allow to give a list of pangenome

def launch(args: argparse.Namespace):
    """Allow to launch the panfluidity script from panorama command

    :param args: list of arguments
    """
    pangenome = Pangenome()
    pangenome.addFile(args.pangenome)
    g = genomes_fluidity(p_pangenome=pangenome, disable_bar=args.disable_prog_bar)
    write_fluidity(pangenome, g)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Parser arguments specific to panfluidity command
    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("panfluidity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")

    return parser
