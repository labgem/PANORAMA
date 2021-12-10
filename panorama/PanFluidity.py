#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
# import statistics
import sys
from time import time

# installed libraries
import pkg_resources
from gmpy2 import popcount
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo
from itertools import combinations
from tqdm import tqdm


# local libraries


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
    print(time()-t4)
    logging.getLogger().info("Compute rate of unique family for each genome combination")
    for c_organisms in tqdm(list(combinations(p_pangenome.organisms, 2)), unit="combination", disable=disable_bar):
        # t1.append(time())
        tot_fam = org2_nb_fam.get(c_organisms[0].name) + org2_nb_fam.get(c_organisms[1].name)
        # t2.append(time())
        common_fam = popcount(c_organisms[0].bitarray & c_organisms[1].bitarray) - 1
        # t3.append(time())
        g_sum += (tot_fam - 2*common_fam) / tot_fam
    # dif1 = [b - a for a, b in zip(t1, t2)]
    # dif2 = [b - a for a, b in zip(t2, t3)]

    # print(statistics.mean(dif1), statistics.mean(dif2))
    print("Temps de calcul : ", time()-t0)
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


# TODO Function to normalize genome fluidity
# def genomes_fluidity_norm(p_pangenome: Pangenome, disable_bar: bool = False) -> float:

# TODO Function to compute mash distance between genome

# TODO Function to export results


def launch(p_args: argparse.Namespace):
    pangenome = Pangenome()
    pangenome.addFile(p_args.pangenome)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute the pan-genome fluidity  from one pangenome",
        formatter_class=argparse.RawTextHelpFormatter
    )

    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required")
    required.add_argument('-p', '--pangenome', required=True, type=str, help="The pangenome .h5 file")
    required.add_argument('-o', '--output', required=True, type=str,
                          help="Output directory where the file(s) will be written")

    optional = parser.add_argument_group(title="Optional arguments",
                                         description="The following arguments are optional")
    optional.add_argument("--disable_bar", required=False, action="store_true",
                          help="disables the progress bars")
    optional.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                          help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")

    args = parser.parse_args()

    level = logging.INFO  # info, warnings and errors, default verbose == 1
    if hasattr(args, "verbose"):
        if args.verbose == 2:
            level = logging.DEBUG  # info, debug, warnings and errors
        elif args.verbose == 0:
            level = logging.WARNING  # only warnings and errors

    logging.basicConfig(level=level,
                        format='%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logging.getLogger().info("Command: " + " ".join([arg for arg in sys.argv]))
    logging.getLogger().info("PPanGGOLiN version: " + pkg_resources.get_distribution("ppanggolin").version)
    # launch(args)

    # printInfo(pangenome=args.pangenome, content=True)

    p = Pangenome()
    p.addFile(args.pangenome)
    g = genomes_fluidity(p_pangenome=p, disable_bar=args.disable_bar)
    print(g)
