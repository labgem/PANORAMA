# default libraries
import os
import logging
from pathlib import Path
import csv

# installed libraries
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import write_info
import tables


# File managing system
def mkdir(output: str, force: bool = False) -> Path:
    """Create a directory at the given path

    :param output:
    :param force:

    :raise FileExistError: If the given directory already exist and no force is used
    :raise Exception: Handle all others exception

    :return: Path object to output directory
    """
    try:
        os.makedirs(output)
    except OSError:
        if not force:
            raise FileExistsError(f"{output} already exists."
                                  f"Use --force if you want to overwrite the files in the directory")
        else:
            logging.getLogger().warning(f"{output} already exist and file will be overwrite by the new generated")
            return Path(output)
    except Exception:
        raise Exception("An unexpected error happened. Please report on our GitHub")
    else:
        return Path(output)


def check_tsv_sanity(tsv_path: Path):
    """ Check if the given tsv is readable for the next PANORAMA step

    :param tsv_path:
    """
    pan_to_path = {}
    try:
        file = open(tsv_path.absolute(), 'r')
        tsv = csv.reader(file, delimiter="\t")
    except IOError:
        raise IOError
    except Exception:
        raise Exception("Unexpected error when opening the list of pangenomes")
    else:
        for line in tsv:
            if len(line) < 2:
                raise Exception("Format not readable. You need at least 2 columns (name and path to pangenome)")
            if " " in line[0]:
                raise Exception(f"Your pangenome names contain spaces (The first encountered pangenome name that had "
                                f"this string: '{line[0]}'). To ensure compatibility with all of the dependencies of "
                                f"PPanGGOLiN this is not allowed. Please remove spaces from your genome names.")
            if not Path(line[1]).exists():
                raise IOError(f"The given path {line[1]} not exist")
            pan_to_path[line[0]] = Path(line[1])
        file.close()
        return pan_to_path


def write_pangenome(pangenome: Pangenome, g_fluidity: float = None):
    """Write or update information in pangenome

    :param pangenome: pangenome where it should be written the genomes fluidity
    :param g_fluidity: genomes_fluidity computed for the pangenome
    """

    with tables.open_file(pangenome.file, "a") as h5f:
        if "/info" not in h5f:
            logging.getLogger().info("Your pangenome is without information. Information will be write to be "
                                     "relevant for analyses")
            write_info(pangenome, h5f)

        info_group = h5f.root.info
        if g_fluidity is not None:
            info_group._v_attrs.fluidity = g_fluidity
    logging.getLogger().info("Done writing the genome fluidity in pangenome")


def check_file_info(pangenome_file: str, need_fluidity: bool = False) -> bool:
    exists = os.path.isfile(pangenome_file)
    if exists:
        with tables.open_file(pangenome_file, "r") as h5f:
            info = h5f.root.info._v_attrs._f_list()
            if need_fluidity:
                return True if 'fluidity' in info else False
    else:
        raise FileNotFoundError(f"The {pangenome_file} does not exist")