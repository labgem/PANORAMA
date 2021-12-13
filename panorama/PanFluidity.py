#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from multiprocessing import get_context
import os

# installed libraries
import ppanggolin.formats
from gmpy2 import popcount
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import checkPangenomeInfo
from itertools import combinations
from tqdm import tqdm
import tables
import pandas as pd


def mkOutdir(output, force):
    if not os.path.exists(output):
        os.makedirs(output)
    elif not force:
        raise FileExistsError(f"{output} already exists. Use -f if you want to overwrite the files in the directory")


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
    g_sum = 0
    logging.getLogger().debug("Get number of families in each organisms")
    org2_nb_fam = nb_fam_per_org(pangenome, disable_bar)
    logging.getLogger().info("Compute rate of unique family for each genome combination")
    for c_organisms in tqdm(list(combinations(pangenome.organisms, 2)), unit="combination", disable=disable_bar):
        tot_fam = org2_nb_fam.get(c_organisms[0].name) + org2_nb_fam.get(c_organisms[1].name)
        common_fam = popcount(c_organisms[0].bitarray & c_organisms[1].bitarray) - 1
        g_sum += (tot_fam - 2 * common_fam) / tot_fam
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

    with tables.open_file(pangenome.file, "a") as h5f:
        if "/info" not in h5f:
            logging.getLogger().info("Your pangenome is without information. Information will be write to be "
                                     "relevant with genomes fluidity")
            ppanggolin.formats.writeInfo(pangenome, h5f)

        info_group = h5f.root.info
        info_group._v_attrs.fluidity = g_fluidity
    logging.getLogger().info("Done writing the genome fluidity in pangenome")


# TODO Function to normalize genome fluidity
# def genomes_fluidity_norm(p_pangenome: Pangenome, disable_bar: bool = False) -> float:

# TODO Function to compute mash distance between genome

# TODO Function to export results (graphique avec coloration par genre ...)
def export_tsv(pangenomes: list, output: str, taxonomies: str = None):
    export_dict = {}
    for pangenome in pangenomes:
        with tables.open_file(pangenome, "r") as h5f:
            infoGroup = h5f.root.info
            if "numberOfOrganisms" in infoGroup._v_attrs._f_list():
                export_dict[pangenome] = [infoGroup._v_attrs['fluidity'],
                                          infoGroup._v_attrs['numberOfOrganisms']]
            else:
                raise Exception("The number of organism is not write in pangenome. Please provide a pangenome with the "
                                "number of organisms write")
    export_df = pd.DataFrame.from_dict(export_dict, orient='index', columns=['Genomes fluidity',
                                                                             "Nb of Organisms"])
    export_df.index.name = "Pangenome"
    if taxonomies is not None:
        tax_df = pd.read_csv(taxonomies, sep="\t", header=0, index_col=0)
        export_df = export_df.merge(tax_df, left_index=True, right_index=True)
    export_df.to_csv(f"{output}/genomes_fluidity.tsv", sep="\t", header=True, index=True, decimal=",",
                     float_format='%.3f')


def check_file_info(pangenome_file: str) -> bool:
    exists = os.path.isfile(pangenome_file)
    if exists:
        with tables.open_file(pangenome_file, "r") as h5f:
            if 'fluidity' in h5f.root.info._v_attrs._f_list():
                return True
            else:
                return False
    else:
        raise FileNotFoundError(f"The {pangenome_file} does not exist")


def skim_pangenome(pangenomes: list) -> list:
    """ Select the pangenome without genomes fluidity computed

    :param pangenomes: list of all the pangenomes

    :return: list of pangenomes without genomes fluidity computed
    """
    logging.getLogger().debug("Get pangenome without genomes fluidity")
    pangenomes_list = []
    for pangenome in pangenomes:
        if not check_file_info(pangenome_file=pangenome):
            pangenomes_list.append(pangenome)
    logging.getLogger().debug("Selection done")
    return pangenomes_list


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
    mkOutdir(output=args.output, force=args.force)
    if len(args.pangenomes) == 0:
        raise Exception("Not one pangenome was found. Please check path to pangenome files.")

    if args.cpu == 1:
        for pangenome_file in skim_pangenome(args.pangenomes):
            logging.getLogger().debug(f"Begin {pangenome_file}")
            pangenome = Pangenome()
            pangenome.addFile(pangenome_file)
            g = genomes_fluidity(pangenome=pangenome, disable_bar=args.disable_prog_bar)
            write_fluidity(pangenome, g)
            logging.getLogger().debug(f"{pangenome_file} Done")
    else:
        logging.getLogger().info("Begin of the genomes fluidity computation")
        skim_list = skim_pangenome(args.pangenomes)
        with get_context('fork').Pool(args.cpu) as p:
            for pangenome_file in tqdm(p.imap_unordered(mp_genome_fluidity, skim_list), unit='pangenome',
                                       total=len(skim_list), disable=args.disable_prog_bar):
                logging.getLogger().debug(f"{pangenome_file} Done")
        logging.getLogger().info("All the genomes fluidity were computed")
    export_tsv(args.pangenomes, args.output, args.taxonomies)


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
    optional.add_argument("-t", "--taxonomies", required=False, default=None, type=str,
                          help="Taxonomies information corresponding to pangenomes. "
                               "Field must be separate by tab with header.")

    return parser
