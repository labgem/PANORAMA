# default libraries
import os
import logging
from ppanggolin.pangenome import Pangenome
from ppanggolin.formats import writeInfo
import tables


# File managing system
def mkdir(output, force):
    if not os.path.exists(output):
        os.makedirs(output)
    elif not force:
        raise FileExistsError(f"{output} already exists. Use -f if you want to overwrite the files in the directory")


def write_pangenome(pangenome: Pangenome, g_fluidity: float = None):
    """Write or update information in pangenome

    :param pangenome: pangenome where it should be written the genomes fluidity
    :param g_fluidity: genomes_fluidity computed for the pangenome
    """

    with tables.open_file(pangenome.file, "a") as h5f:
        if "/info" not in h5f:
            logging.getLogger().info("Your pangenome is without information. Information will be write to be "
                                     "relevant for analyses")
            writeInfo(pangenome, h5f)

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