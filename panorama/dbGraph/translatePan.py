#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
from tqdm import tqdm
from pathlib import Path
from typing import Union
import uuid

# installed librairies
from ppanggolin.genome import Organism, Contig
from ppanggolin.edge import Edge
from ppanggolin.region import Module, Region, Spot

# local librairies
from panorama.pangenomes import Pangenome
from panorama.geneFamily import GeneFamily
from panorama.format.read_binaries import check_pangenome_info

fam_visit = set()


def give_gene_tmp_id(pangenome: Pangenome, disable_bar: bool = True):
    for gene in tqdm(pangenome.genes, total=pangenome.number_of_gene(), unit='gene', disable=disable_bar):
        gene.tmp_id = str(uuid.uuid4())


def get_genes(parent: Union[GeneFamily, Region, Contig]):
    genes_list = []
    for gene in parent.genes:
        genes_list.append({"name": gene.name,
                           "genomic_type": gene.type,
                           "is_fragment": bool(gene.is_fragment),
                           # Prevent TypeError: Object of type bool_ is not JSON serializable
                           "start": gene.start,
                           "stop": gene.stop,
                           "strand": gene.strand,
                           "product": gene.product,
                           "local_id": gene.local_identifier,
                           "tmp_id": gene.tmp_id})
    return genes_list


def write_contig(organism: Organism):
    contig_list = []
    for contig in organism.contigs:
        contig_list.append({"name": contig.name,
                            "is_circular": bool(contig.is_circular),
                            "Gene": get_genes(contig)})
    return contig_list


def write_families(pangenome: Pangenome, out_dict: dict):
    family: GeneFamily
    for family in pangenome.gene_families:
        annot = family.get_annot("CARD")
        fam_dic_property = {"name": family.name,
                            "Partition": {"partition": family.named_partition,
                                          "subpartition": family.partition},
                            "annotation": annot[0] if annot is not None else None,
                            "Gene": get_genes(family),
                            "Family": get_neighbor(family)  # Return neighbor with edges weight
                            }
        # fam_dic_property.update(get_neighbor(family))
        # if family.named_partition == "persistent":
        #     out_dict["Pangenome"]["Family"]["Persistent"].append(fam_dic_property)
        # elif family.named_partition == "shell":
        #     out_dict["Pangenome"]["Family"]["Shell"].append(fam_dic_property)
        # elif family.named_partition == "cloud":
        #     out_dict["Pangenome"]["Family"]["Cloud"].append(fam_dic_property)
        # else:
        #     raise Exception("Unrecognized partition name")
        out_dict["Pangenome"]["Family"].append(fam_dic_property)


def get_neighbor(family: GeneFamily):
    global fam_visit
    edge: Edge
    neighbors = []
    for edge in family.edges:
        if edge.source.name == family.name:
            if edge.target.name not in fam_visit:
                annot = edge.target.get_annot("CARD")
                neighbors.append({"weight": len(edge.organisms),
                                  "name": edge.target.name,
                                  "Partition": {"partition": edge.target.named_partition,
                                                "subpartition": edge.target.partition},
                                  "annotation": annot[0] if annot is not None else None})
        elif edge.target.name == family.name:
            if edge.source.name not in fam_visit:
                annot = edge.source.get_annot("CARD")
                neighbors.append({"weight": len(edge.organisms),
                                  "name": edge.source.name,
                                  "Partition": {"partition": edge.source.named_partition,
                                                "subpartition": edge.source.partition},
                                  "annotation": annot[0] if annot is not None else None})
        else:
            raise Exception("Source and target name are different from edge's family. "
                            "Please check you import graph data and if the problem persist, post an issue.")
    fam_visit.add(family.name)

    return neighbors


def write_organisms(pangenome: Pangenome, out_dict: dict):
    for org in pangenome.organisms:
        out_dict["Pangenome"]["Genome"].append({"name": str(org.name),
                                                "Contig": write_contig(org)})


def write_rgp(parent: Union[Pangenome, Spot]):
    rgp_list = []
    for rgp in parent.regions:
        rgp_list.append({"name": str(rgp.name),
                         "start": rgp.start,
                         "stop": rgp.stop,
                         "score": float(rgp.score),
                         "is_whole_contig": rgp.is_whole_contig,
                         "is_contig_border": rgp.is_contig_border,
                         "Gene": get_genes(rgp)})
    return rgp_list


def write_spot(pangenome: Pangenome, out_dict: dict):
    for spot in pangenome.spots:
        out_dict["Pangenome"]["Spot"].append({"name": f"{str(spot.ID)}",
                                              "RGP": write_rgp(parent=spot)})


def write_modules(pangenome: Pangenome, out_dict: dict):
    module: Module
    for module in pangenome.modules:
        module_dict = {"name": int(module.ID),  # int prevent TypeError: Object of type uint32 is not JSON serializable
                       "Family": []}
        for family in module.families:
            annot = family.get_annot("CARD")
            module_dict["Family"].append({"name": family.name,
                                          "Partition": {"partition": family.named_partition,
                                                        "subpartition": family.partition},
                                          "annotation": annot[0] if annot is not None else None})
        out_dict["Pangenome"]["Module"].append(module_dict)


def create_dict(pangenome: Pangenome):
    out_dict = {"Pangenome": {"name": pangenome.name, "taxid": pangenome.taxid,
                              "Family": [],
                              "Partition": [],
                              "Module": [],
                              "RGP": write_rgp(parent=pangenome),
                              "Spot": [], "Genome": []}}
    write_families(pangenome, out_dict)
    write_organisms(pangenome, out_dict)
    write_spot(pangenome, out_dict)
    write_modules(pangenome, out_dict)
    return out_dict


def write_json(pangenome: Pangenome, output: Path, compress):
    logging.getLogger().info(f"Writing the json file for the {pangenome.name} pangenome graph")
    outname = f"{output.absolute()}/{pangenome.name}_short2.json"
    out_dict = create_dict(pangenome)

    # with write_compressed_or_not(outname, compress) as json_file:
    #     json.dump(out_dict, json_file, indent=4)


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pangenome = Pangenome(name=args.pangenome.stem)
    pangenome.add_file(args.pangenome)
    check_pangenome_info(pangenome, need_annotations=True, need_families=True, need_graph=True, need_partitions=True,
                         need_rgp=True, need_spots=True, need_modules=True, need_anntation_fam=True)
    give_gene_tmp_id(pangenome)
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
