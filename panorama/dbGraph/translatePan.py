#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import json
import tempfile
from pathlib import Path

from ppanggolin.edge import Edge
# installed librairies
from ppanggolin.utils import write_compressed_or_not
from ppanggolin.formats.readBinaries import check_pangenome_info
from ppanggolin.genome import Organism, Contig
from ppanggolin.region import Module, Region, Spot
# local librairies
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily

index = 1
fam_visit = set()


def get_genes(family: GeneFamily):
    genes_list = []
    for gene in family.genes:
        genes_list.append({"name": gene.name,
                           "genomic_type": gene.type,
                           "is_fragment": bool(gene.is_fragment),
                           # Prevent TypeError: Object of type bool_ is not JSON serializable
                           "start": gene.start,
                           "stop": gene.stop,
                           "strand": gene.strand,
                           "product": gene.product})
    return genes_list


def write_gene_to_contig(contig: Contig, out_dict: dict):
    for gene in contig.genes:
        out_dict["graph"]["edges"].append({"from": gene.ID, "to": contig.name,
                                           "type": ["IN_CONTIG"], "attr": {}})


def write_contig_to_organism(organism: Organism, out_dict: dict):
    for contig in organism.contigs:
        out_dict["graph"]["nodes"].append({"id": f"C_{contig.name}",
                                           "attr": {"is_circular": int(contig.is_circular)}})
        out_dict["graph"]["edges"].append({"from": f"C_{contig.name}", "to": organism.name,
                                           "type": ["IN_GENOME"], "attr": {}})
        out_dict["graph"]["node_types"][f"C_{contig.name}"] = ["CONTIG"]
        write_gene_to_contig(contig, out_dict)


def write_families(pangenome: Pangenome, out_dict: dict):
    for family in pangenome.gene_families:
        out_dict["Pangenome"]["Family"].append({"name": family.name,
                                                "partition": family.named_partition,
                                                "subpartition": family.partition,
                                                "Gene": get_genes(family),
                                                "Family": get_neighbor(family)  # Return neighbor with edges weight
                                                })


def get_neighbor(family: GeneFamily):
    global fam_visit
    edge: Edge
    neighbor_list = []
    for edge in family.edges:
        if edge.source.name == family.name:
            if edge.target.name not in fam_visit:
                neighbor_list.append({"weight": len(edge.organisms),
                                      "name": edge.target.name,
                                      "partition": edge.target.named_partition,
                                      "subpartition": edge.target.partition
                                      })
        elif edge.target.name == family.name:
            if edge.source.name not in fam_visit:
                neighbor_list.append({"weight": len(edge.organisms),
                                      "name": edge.source.name,
                                      "partition": edge.source.named_partition,
                                      "subpartition": edge.source.partition
                                      })
        else:
            raise Exception("Source and target name are different from edge's family. "
                            "Please check you import graph data and if the problem persist, post an issue.")
    fam_visit.add(family.name)
    return neighbor_list


def write_organisms(pangenome: Pangenome, out_dict: dict):
    for org in pangenome.organisms:
        out_dict["graph"]["nodes"].append({"id": str(org.name), "attr": {}})
        out_dict["graph"]["node_types"][str(org.name)] = ["GENOME"]

        write_contig_to_organism(organism=org, out_dict=out_dict)


def write_in_rgp(region: Region, out_dict: dict):
    for gene in region.genes:
        out_dict["graph"]["edges"].append({"from": gene.ID, "to": str(region.name),
                                           "type": ["IN_RGP"],
                                           "attr": {"start": True if gene == region.start_gene else False,
                                                    "stop": True if gene == region.stop_gene else False}})


def write_rgp(pangenome: Pangenome, out_dict: dict):
    for rgp in pangenome.regions:
        out_dict["graph"]["nodes"].append({"id": str(rgp.name),
                                           "attr": {"start": rgp.start, "stop": rgp.stop, "score": float(rgp.score),
                                                    "is_whole_contig": rgp.is_whole_contig,
                                                    "is_contig_border": rgp.is_contig_border}})
        out_dict["graph"]["node_types"][str(rgp.name)] = ["RGP"]

        write_in_rgp(region=rgp, out_dict=out_dict)


def write_in_spot(spot: Spot, out_dict: dict):
    for region in spot.regions:
        out_dict["graph"]["edges"].append({"from": f"{str(region.name)}", "to": f"S_{str(spot.ID)}",
                                           "type": ["IN_SPOT"],
                                           "attr": {}})


def write_spot(pangenome: Pangenome, out_dict: dict):
    for spot in pangenome.spots:
        out_dict["graph"]["nodes"].append({"id": f"S_{str(spot.ID)}", "attr": {}})
        out_dict["graph"]["node_types"][f"S_{str(spot.ID)}"] = ["SPOT"]

        write_in_spot(spot=spot, out_dict=out_dict)


def write_modules(pangenome: Pangenome, out_dict: dict):
    module: Module
    for module in pangenome.modules:
        out_dict["Pangenome"]["Module"].append({"module_id": int(module.ID),  # int prevent TypeError: Object of type uint32 is not JSON serializable
                                                "Family": [{"name": family.name,
                                                            "partition": family.named_partition,
                                                            "subpartition": family.partition}
                                                           for family in module.families]
                                                })


def create_dict(pangenome: Pangenome):
    out_dict = {"Pangenome": {"name": pangenome.name, "Family": [], "Edges": [], "Module": []}}
    write_families(pangenome, out_dict)
    # write_neighbor(pangenome, out_dict)
    # write_organisms(pangenome, out_dict)
    # write_rgp(pangenome, out_dict)
    # write_spot(pangenome, out_dict)
    write_modules(pangenome, out_dict)

    return out_dict


def write_json(pangenome: Pangenome, output: Path, compress):
    logging.getLogger().info(f"Writing the json file for the {pangenome.name} pangenome graph")
    outname = f"{output.absolute()}/{pangenome.name}_short2.json"
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
        json.dump(out_dict, json_file, indent=4)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pangenome = Pangenome(name=args.pangenome.stem)
    pangenome.add_file(args.pangenome)
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True, need_partitions=True,
                         need_rgp=True, need_spots=True, need_modules=True)
    write_json(pangenome, args.out_directory, compress=False)
    logging.getLogger().info("Translate pangenome in json Done")


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
