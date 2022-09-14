#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import json
import tempfile
from pathlib import Path

# installed librairies
from ppanggolin.utils import write_compressed_or_not
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.genome import Organism, Contig
# local librairies
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily


def write_gene_to_families(family: GeneFamily, out_dict: dict):
    for gene in family.genes:
        out_dict["graph"]["nodes"].append({"id": gene.ID,
                                           "attr": {"genomic_type": "CDS",
                                                    "is_fragment": int(gene.is_fragment)}})
        out_dict["graph"]["edges"].append({"from": gene.ID, "to": family.ID, "type": "IN_FAMILY", "attr": {}})
        out_dict["graph"]["node_types"][str(gene.ID)] = ["gene"]


def write_gene_to_contig(contig: Contig, out_dict: dict):
    for gene in contig.genes:
        out_dict["graph"]["edges"].append({"from": gene.ID, "to": contig.name, "type": "IN_CONTIG", "attr": {}})


def write_contig_to_organism(organism: Organism, out_dict: dict):
    for contig in organism.contigs:
        out_dict["graph"]["nodes"].append({"id": f"C_{contig.name}",
                                           "attr": {"is_circular": int(contig.is_circular)}})
        out_dict["graph"]["edges"].append({"from": f"C_{contig.name}", "to": organism.name,
                                           "type": "IN_GENOME", "attr": {}})
        out_dict["graph"]["node_types"][f"C_{contig.name}"] = ["gene"]
        write_gene_to_contig(contig, out_dict)


def write_gene_families(pangenome: Pangenome, out_dict: dict):
    for family in pangenome.gene_families:
        out_dict["graph"]["nodes"].append({"id": f"F_{family.ID}", "attr": {"name": family.name,
                                                                            "partition": family.named_partition,
                                                                            "subpartition": family.partition}})
        out_dict["graph"]["node_types"][f"F_{family.ID}"] = ["FAMILY", family.named_partition]
        out_dict["graph"]["edges"].append({"from": f"F_{family.ID}",
                                           "to": pangenome.name,
                                           "type": ["IN_PANGENOME"],
                                           "attr": {"partition": family.named_partition}})
        write_gene_to_families(family, out_dict)


def write_neighbor_edges(pangenome: Pangenome, out_dict: dict):
    for edge in pangenome.edges:
        part_source, part_target = (pangenome.get_gene_family(edge.source.name).named_partition,
                                    pangenome.get_gene_family(edge.target.name).named_partition)
        out_dict["graph"]["edges"].append(
            {"from": f"F_{edge.source.ID}", "to": f"F_{edge.target.ID}",
             "type": ["NEIGHBOR_OF", f"{part_source}_{part_target}"],
             "attr": {"weight": len(edge.organisms)}})


def write_organisms(pangenome: Pangenome, out_dict: dict):
    for org in pangenome.organisms:
        out_dict["graph"]["nodes"].append({"id": str(org.name), "attr": {}})
        out_dict["graph"]["node_types"][str(org.name)] = ["GENOME"]


def create_dict(pangenome: Pangenome):
    out_dict = {"graph": {"edges": [], "nodes": [], "node_types": {}}}
    out_dict["graph"]["nodes"].append({"id": pangenome.name, "attr": {"parameters": pangenome.parameters}})
    out_dict["graph"]["node_types"][pangenome.name] = ["PANGENOME"]
    write_gene_families(pangenome, out_dict)
    write_neighbor_edges(pangenome, out_dict)
    write_organisms(pangenome, out_dict)

    return out_dict


def write_json(pangenome: Pangenome, output: Path, compress):
    logging.getLogger().info(f"Writing the json file for the {pangenome.name} pangenome graph")
    outname = f"{output.absolute()}/{pangenome.name}.json"
    out_dict = create_dict(pangenome)

    # for mod in pangenome.modules:
    #     out_dict["graph"]["nodes"].append(
    #         {"id": pan_name + '_' + str(mod.ID), "attr": {"nb_fams": str(len(mod.families))}})
    #     out_dict["graph"]["node_types"][pan_name + '_' + str(mod.ID)] = ["Module"]
    #     for family in mod.families:
    #         out_dict["graph"]["edges"].append(
    #             {"from": "F_" + family.name, "to": pan_name + '_' + str(mod.ID), "type": ["IN_MODULE"], "attr": {}})
    #         out_dict["graph"]["edges"].append(
    #             {"from": pan_name + '_' + str(mod.ID), "to": pan_name, "type": ["IN_TAXA"], "attr": {}})
    with write_compressed_or_not(outname, compress) as json_file:
        json.dump(out_dict, json_file)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pangenome = Pangenome(name=args.pangenome.stem)
    pangenome.add_file(args.pangenome)
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True)
    write_json(pangenome, args.out_directory, compress=False)


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("annotation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_export_db(parser)
    return parser


def parser_export_db(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('-p', '--pangenome', required=True, type=Path, nargs='?',
                          help='a pangenome.h5 file')
    required.add_argument('-o', '--out_directory', required=True, type=Path, nargs='?',
                          help='output directory to store json file')


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_export_db(main_parser)
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
