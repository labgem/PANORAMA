#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import json
import tempfile
from pathlib import Path

# installed librairies

# local librairies
from panorama.pangenomes import Pangenome


def write_json(pangenome: Pangenome, output, compress, pan_name: str):
    logging.getLogger().info("Writing the json file for the pangenome graph...")
    outname = output + "/pangenomeGraph.json"
    out_dict = {"graph": {"edges": [], "nodes": [], "node_types": {}}}
    out_dict["graph"]["nodes"].append({"id": pangenome.name, "attr": {"parameters"}})
    out_dict["graph"]["node_types"][pan_name] = ["Taxa"]
    for org in pangenome.organisms:
        out_dict["graph"]["nodes"].append({"id": str(org.name), "attr": {}})
        out_dict["graph"]["node_types"][str(org.name)] = ["genome"]
    for edge in pangenome.edges:
        part_source, part_target = (pangenome.get_gene_family(edge.source.name).named_partition,
                                    pangenome.get_gene_family(edge.target.name).named_partition)
        out_dict["graph"]["edges"].append(
            {"from": "F_" + edge.source.name, "to": "F_" + edge.target.name,
             "type": ["NEIGHBOR_OF", f"{part_source}_{part_target}"],
             "attr": {"weight": len(edge.organisms)}})
    for geneFam in list(pangenome.geneFamilies):
        out_dict["graph"]["nodes"].append({"id": "F_" + geneFam.name, "attr": {"nb_genomes": len(edge.organisms),
                                                                               "partition": geneFam.namedPartition,
                                                                               "subpartition": geneFam.partition,
                                                                               "nb_genes": len(geneFam.genes)}})
        out_dict["graph"]["node_types"]["F_" + geneFam.name] = ["GeneFamily", geneFam.namedPartition]
        out_dict["graph"]["edges"].append(
            {"from": "F_" + geneFam.name, "to": pan_name,
             "type": ["IN_TAXA"],
             "attr": {}})
        for gene in geneFam.genes:
            out_dict["graph"]["nodes"].append(
                {"id": gene.ID, "attr": {"genomic_type": "CDS", "is_fragment": int(gene.is_fragment)}})
            out_dict["graph"]["edges"].append({"from": gene.ID, "to": geneFam.ID, "type": "IN_FAMILY", "attr": {}})
            out_dict["graph"]["node_types"][str(gene.ID)] = ["gene"]
            out_dict["graph"]["edges"].append(
                {"from": gene.ID, "to": str(gene.organism.name), "type": "IN_ORG", "attr": {}})
    for mod in pangenome.modules:
        out_dict["graph"]["nodes"].append(
            {"id": pan_name + '_' + str(mod.ID), "attr": {"nb_fams": str(len(mod.families))}})
        out_dict["graph"]["node_types"][pan_name + '_' + str(mod.ID)] = ["Module"]
        for family in mod.families:
            out_dict["graph"]["edges"].append(
                {"from": "F_" + family.name, "to": pan_name + '_' + str(mod.ID), "type": ["IN_MODULE"], "attr": {}})
            out_dict["graph"]["edges"].append(
                {"from": pan_name + '_' + str(mod.ID), "to": pan_name, "type": ["IN_TAXA"], "attr": {}})
    with write_compressed_or_not(outname, compress) as json_file:
        json.dump(out_dict, json_file)


def writeJSONGeneFam(geneFam, json):
    json.write('{' + f'"id": {geneFam.ID}, "attr":' + '{' +
               f'"name": "{geneFam.name}", ' +
               f'"nb_genes": {len(geneFam.genes)}, ' +
               f'"partition": "{geneFam.namedPartition}", "subpartition": "{geneFam.partition}"' + '}')


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    pass


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
    required.add_argument('-s', '--systems', required=True, type=Path,
                          help="Path to systems directory")
    required.add_argument('-p', '--pangenomes', required=True, type=Path, nargs='?',
                          help='A list of pangenome .h5 files in .tsv file')
    required.add_argument("--source", required=True, type=str, nargs="?",
                          help='Name of the annotation source. Default use name of annnotation file or directory.')
    exclusive_mode = required.add_mutually_exclusive_group(required=True)
    exclusive_mode.add_argument('--tsv', type=Path, nargs='?',
                                help='Gene families annotation in TSV file. See our github for more detail about format')
    exclusive_mode.add_argument('--hmm', type=Path, nargs='?',
                                help="File with all HMM or a directory with one HMM by file")
    exclusive_threshold = required.add_mutually_exclusive_group()
    exclusive_threshold.add_argument('--e_value', required=False, type=float, default=None,
                                     help="Set the same e-value for all hmm to assign function")
    exclusive_threshold.add_argument('--cutoffs', required=False, type=Path, default=None,
                                     help="TSV file to set one evalue for each hmm")
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--prediction_size", required=False, type=int, default=1,
                          help="Number of prediction associate with gene families")
    optional.add_argument("--tmpdir", required=False, type=str, nargs='?', default=Path(tempfile.gettempdir()),
                          help="directory for storing temporary files")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of av available threads")


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
