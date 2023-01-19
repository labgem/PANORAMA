#!/usr/bin/env python3
# coding:utf-8


# default libraries
from __future__ import annotations
import argparse
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from typing import Dict, List, Set

# installed libraries
import networkx as nx
from ppanggolin.context.searchGeneContext import compute_gene_context_graph

# local libraries
from panorama.models import Models, Model, FuncUnit
from panorama.utils import check_tsv_sanity
from panorama.format.read_binaries import check_pangenome_info
from panorama.format.write_binaries import write_pangenome, erase_pangenome
from panorama.system import System
from panorama.region import Module
from panorama.geneFamily import GeneFamily
from panorama.pangenomes import Pangenome


def check_pangenome_detection(pangenome: Pangenome, source: str, force: bool = False, disable_bar: bool = False):
    """ Check and load pangenome information before adding annotation

    :param pangenome: Pangenome object
    :param source: Source used to detect system
    :param disable_bar: Disable progress bar

    :raise KeyError: Provided source is not in the pangenome
    """
    if "systems_sources" in pangenome.status and source in pangenome.status["systems_sources"]:
        if force:
            erase_pangenome(pangenome, systems=True, source=source)
        else:
            raise Exception(f"Systems are already detected based on the annotation source : {source}."
                            f" Use the --force option to erase the already computed systems.")
    try:
        _ = pangenome.status["annotations_sources"]
    except Exception:
        raise Exception("There is no annotation in your pangenome. "
                        "Please see the command annotation before to detect systems")
    else:
        if source in pangenome.status["annotations_sources"]:
            check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_annotations_fam=True,
                                 sources=[source], disable_bar=disable_bar)
        else:
            raise KeyError("Annotation source not in pangenome.")


def read_models(models_path: Path, disable_bar: bool = False) -> Models:
    """Read all json files models in the directory

    :param models_path: path of models directory
    :param disable_bar: Disable progress bar

    :raise KeyError: One or more keys are missing or non-acceptable
    :raise TypeError: One or more value are not with good type
    :raise ValueError: One or more value are not non-acceptable
    :raise Exception: Manage unexpected error
    """
    models = Models()
    models.read(models_path, disable_bar)
    return models


def get_annotation_to_families(pangenome: Pangenome, source: str) -> Dict[str, Set[GeneFamily]]:
    """ Get for eache annotation a set of families with this annotation

    :param pangenome: Pangenome with gene families
    :param source: name of the annotation source

    :return: Dictionnary with for each annotation a set of gene families
    """
    annot2fam = {}
    for gf in pangenome.gene_families:
        annotations = gf.get_source(source)
        if annotations is not None:
            for annotation in annotations:
                if annotation.name in annot2fam:
                    annot2fam[annotation.name].add(gf)
                else:
                    annot2fam[annotation.name] = {gf}
    return annot2fam


def dict_families_context(func_unit: FuncUnit, annot2fam: dict) -> (dict, dict):
    """Recover all families in the function unit

    :param func_unit: function unit object of model
    :param annot2fam: dictionary of annotated families

    :return families: dictionary of families interesting
    """
    families = dict()
    fam2annot = dict()
    for fam_model in func_unit.families:
        if fam_model.name in annot2fam:
            for gf in annot2fam[fam_model.name]:
                families[gf.name] = gf
                if gf.name in fam2annot:
                    fam2annot[gf.name].add(fam_model)
                else:
                    fam2annot[gf.name] = {fam_model}
    return families, fam2annot


def search_fu_with_one_fam(func_unit: FuncUnit, annot2fam: dict, source: str, graph: nx.Graph) -> List[System]:
    """Search defense system with only one family necessary

    :param func_unit: Functional unit
    :param annot2fam: Dictionnary which link annotation with gene family
    :param source: name of the annotation source
    :param graph: connected component graph

    :return: Detected system in pangenome
    """
    detected_systems = []
    for mandatory_fam in func_unit.mandatory:
        if mandatory_fam.name in annot2fam:
            for pan_fam in annot2fam[mandatory_fam.name]:
                if pan_fam not in graph.nodes:
                    detected_systems.append(System(system_id=0, model=func_unit.model, source=source, gene_families={pan_fam}))
    return detected_systems


def verify_param(g: nx.Graph(), fam2annot: dict, model: Model, func_unit: FuncUnit, source: str,
                 detected_systems: list):
    """Check if the models parameters are respected

    :param g: Connected component graph
    :param fam2annot: dictionary of families interesting
    :param model: Defined model
    :param func_unit: Functional unit
    :param source: annotation source
    :param detected_systems: detected system list
    """
    seen = set()

    def extract_cc(node: GeneFamily, graph: nx.Graph, seen: set) -> Set[GeneFamily]:
        """ Get connected component of the gene family

        :param node: node corresponding to gene family
        :param graph: graph of connected component
        :param seen: set of already check gene families

        :return: set of connected component
        """
        nextlevel = {node}
        cc = set()
        while len(nextlevel) > 0:
            thislevel = nextlevel
            nextlevel = set()
            for v in thislevel:
                if v not in seen:
                    cc.add(v)
                    seen.add(v)
                    nextlevel |= set(graph.neighbors(v))
        return cc

    def check_cc(cc: set, fam2annot: dict, func_unit: FuncUnit) -> bool:
        """Check parameters

        :param cc: set of connected component
        :param fam2annot: link between gene family and annotation
        :param func_unit:functional unit

        :return: Boolean true if parameter respected
        """
        count_forbidden, count_mandatory, count_accesory = (0, 0, 0)
        forbidden_list, mandatory_list, accessory_list = (list(map(lambda x: x.name, func_unit.forbidden)),
                                                          list(map(lambda x: x.name, func_unit.mandatory)),
                                                          list(map(lambda x: x.name, func_unit.accessory)))
        for node in cc:
            annotations = fam2annot.get(node.name)
            if annotations is not None:
                for annot in annotations:
                    if annot.type == 'forbidden' and annot.name in forbidden_list:  # if node is forbidden
                        count_forbidden += 1
                        forbidden_list.remove(annot.name)
                        if count_forbidden > func_unit.parameters["max_forbidden"]:
                            return False
                    elif annot.type == 'mandatory' and annot.name in mandatory_list:  # if node is mandatory
                        count_mandatory += 1
                        mandatory_list.remove(annot.name)
                    elif annot.type == 'accessory' and annot.name in accessory_list:  # if node is accessory
                        count_accesory += 1
                        accessory_list.remove(annot.name)
        if count_mandatory >= func_unit.parameters['min_mandatory'] and \
                count_accesory + count_mandatory >= func_unit.parameters['min_total']:
            return True
        else:
            return False

    remove_node = set()
    for node in g.nodes():  # pour chaque noeud du graphe
        if node not in seen:
            cc = extract_cc(node, g, seen)  # extract coonnect component
            if check_cc(cc, fam2annot, func_unit):
                detected_systems.append(System(system_id=0, model=model, gene_families=cc, source=source))
            else:
                remove_node |= cc
    g.remove_nodes_from(remove_node)


def search_system(model: Model, annot2fam: dict, source: str) -> List[System]:
    """Search if model system is in pangenome

    :param model: model to search
    :param annot2fam: Dictionnary with for each annotation a set of gene families
    :param source: name of the annotation source

    :return: Systems detected
    """
    for func_unit in model.func_units:
        detected_systems = []
        families, fam2annot = dict_families_context(func_unit, annot2fam)
        g = compute_gene_context_graph(families, func_unit.parameters["max_separation"] + 1, disable_bar=True)
        if func_unit.parameters['min_total'] == 1:
            detected_systems += search_fu_with_one_fam(func_unit, annot2fam, source, g)
        verify_param(g, fam2annot, model, func_unit, source, detected_systems)
        return detected_systems


def search_systems(models: Models, pangenome: Pangenome, source: str, threads: int = 1, disable_bar: bool = False):
    """
    Search present model in the pangenome
    :param models: Models to search in pangenomes
    :param pangenome: Pangenome with gene families
    :param source: name of the annotation source
    :param threads: number of available threads
    :param disable_bar: Disable progress bar
    """
    annot2fam = get_annotation_to_families(pangenome=pangenome, source=source)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=models.size, unit='model', disable=disable_bar) as progress:
            futures = []
            for model in models:
                future = executor.submit(search_system, model, annot2fam, source)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            detected_systems = []
            for future in futures:
                result = future.result()
                detected_systems += result
    for system in detected_systems:  # TODO pass in mp step
        pangenome.add_system(system)
    pangenome.status["systems"] = "Computed"


def systems_to_module(module: Module, systems: Set[System]):
    for system in systems:
        if system.gene_families.issubset(module.families):
            module.add_system(system)


def systems_to_modules(pangenome: Pangenome, threads: int = 1, disable_bar: bool = False):
    """Associate a model to modules"""
    with ThreadPoolExecutor(max_workers=threads) as executor:
        with tqdm(total=pangenome.number_of_modules(), unit='module', disable=disable_bar) as progress:
            futures = []
            for module in pangenome.modules:
                future = executor.submit(systems_to_module, module, pangenome.systems)
                future.add_done_callback(lambda p: progress.update())
                futures.append(future)

            for future in futures:
                future.result()


def launch(args):
    """
    Launch functions to detect models in pangenomes

    :param args: Argument given
    """
    pan_to_path = check_tsv_sanity(args.pangenomes)
    models = read_models(args.models)
    for pangenome_name, pangenome_info in pan_to_path.items():
        pangenome = Pangenome(name=pangenome_name, taxid=pangenome_info["taxid"])
        pangenome.add_file(pangenome_info["path"])
        check_pangenome_detection(pangenome, source=args.source, force=args.force, disable_bar=args.disable_prog_bar)
        search_systems(models, pangenome, args.source, args.threads, args.disable_prog_bar)
        logging.getLogger().info("System detection Done")
        logging.getLogger().info(f"Write systems in pangenome")
        # systems_to_modules(pangenome=pangenome, threads=args.threads, disable_bar=args.disable_prog_bar)
        write_pangenome(pangenome, pangenome_info["path"], source=args.source, disable_bar=args.disable_prog_bar)
        # write_systems_projection(pangenome=pangenome, output=args.output, threads=args.threads,
        #                          force=args.force, disable_bar=args.disable_prog_bar)
        logging.getLogger().info(f"Systems written")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("detection", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_detection(parser)
    return parser


def parser_detection(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument('-m', '--models', required=True, type=Path, default=None,
                          help="Path to model directory")
    required.add_argument("-s", "--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source where panorama as to select in pangenomes')
    required.add_argument("-o", "--output", required=True, type=Path, nargs='?',
                          help='Output directory')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_detection(main_parser)
    common = main_parser.add_argument_group(title="Common argument")
    common.add_argument("--verbose", required=False, type=int, default=1, choices=[0, 1, 2],
                        help="Indicate verbose level (0 for warning and errors only, 1 for info, 2 for debug)")
    common.add_argument("--log", required=False, type=check_log, default="stdout", help="log output file")
    common.add_argument("-d", "--disable_prog_bar", required=False, action="store_true",
                        help="disables the progress bars")
    common.add_argument('--force', action="store_true",
                        help="Force writing in output directory and in pangenome output file.")
    set_verbosity_level(main_parser.parse_args())
    launch(main_parser.parse_args())
