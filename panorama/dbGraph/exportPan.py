#!/usr/bin/env python3
# coding:utf-8

# default libraries
import argparse
import logging
import json
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
# installed librairies

# local librairies
from panorama.panDB import Neo4jDB


def add_pangenome_node(node: dict) -> str:
    query = f"CREATE (n:PANGENOME) SET n.id='{node['id']}'\n"
    attr_pan = {"parameters": dict}
    attr_param = {"annotation": {"remove_Overlapping_CDS": bool, "annotate_RNA": bool, "translation_table": str,
                                 "kingdom": str, "contig_filter": int,
                                 "read_annotations_from_file": bool},
                  "cluster": {"coverage": float, "identity": float, "defragmentation": bool,
                              "translation_table": str, "read_clustering_from_file": bool},
                  "graph": {"removed_high_copy_number_families": bool},
                  "partition": {"beta": float, "free_dispersion": bool, "computed_K": bool,
                                "K": int, "max_node_degree_for_smoothing": int},
                  "RGP": {"persistent_penalty": int, "variable_gain": int, "min_length": int,
                          "min_score": int, "dup_margin": float},
                  "spots": {"set_size": int, "overlapping_match": int, "exact_match": int},
                  "modules": {"size": int, "min_presence": int, "transitive": int,
                              "jaccard": float, "dup_margin": float}}
    for attr_key, attr_val in node["attr"].items():
        if attr_key not in attr_pan:
            raise Exception(f"Unknown pangenome property : {str(attr_key)}")
        if not isinstance(attr_val, attr_pan[attr_key]):
            raise Exception(f"Wrong format value for attribute {str(attr_key)}")
        if attr_key == 'parameters':
            for step, attr_step in attr_val.items():
                if not isinstance(attr_step, dict):
                    raise Exception(f"Wrong format value for property {str(step)} in pangenomes parameters")
                else:
                    for param, value in attr_step.items():
                        if not isinstance(value, attr_param[step][param]):
                            raise Exception(f"Wrong format value for property {str(param)}"
                                            f" in pangenomes parameters {step}")
                        if isinstance(value, str):
                            query += f"SET n.{param} = '{value}'\n"
                        else:
                            query += f"SET n.{param} = {value}\n"
    return query

def add_family_node(node: dict, label: list) -> str:
    query = f"CREATE (n:{':'.join(label)}) SET n.id='{node['id']}'\n"
    attr_fam = {"name": str, "partition": ["persistent", "shell", "cloud"], "subpartition": str}
    for attr_key, attr_val in node["attr"].items():
        if attr_key not in attr_fam:
            raise Exception(f"Unknown family property : {str(attr_key)} for family {node['id']}")
        if attr_key == "partition":
            if attr_val not in attr_fam["partition"]:
                raise Exception(f"family {node['id']} partition is not readable. "
                                f"Allow partitions are {attr_fam['partition']}")
        else:
            if not isinstance(attr_val, attr_fam[attr_key]):
                raise Exception(f"Wrong format value for property {attr_key} in family {node['id']}")
        query += f"SET n.{attr_key} = '{attr_val}'\n"
    return query


def add_in_pangenome_edge(edge: dict, pangenome_data: dict):
    family_label = pangenome_data["graph"]["node_types"][edge["from"]]
    query = f"MATCH (s:{':'.join(family_label)} {{id: '{edge['from']}'}}), " \
            f"(t:PANGENOME {{id: '{edge['to']}'}}) \n" \
            "MERGE (s)-[r:IN_PANGENOME]->(t)\n"
    return query


def add_neighbor_of_edge(edge: dict, pangenome_data: dict):
    source_label = pangenome_data["graph"]["node_types"][edge["from"]]
    target_label = pangenome_data["graph"]["node_types"][edge["to"]]
    query = f"MATCH (s:{':'.join(source_label)} {{id: '{edge['from']}'}}), " \
            f"(t:{':'.join(target_label)} {{id: '{edge['to']}'}}) \n" \
            "MERGE (s)-[r:NEIGHBOR_OF]->(t)\n"
    attr_neighbor = {"weight": int}
    for attr_key, attr_val in edge["attr"].items():
        if attr_key not in attr_neighbor:
            raise Exception("Unknown attribute : " + str(attr_key))
        else:
            if not isinstance(attr_val, attr_neighbor[attr_key]):
                raise Exception(f"Wrong format value for property {attr_key} in neighbor edges")
        if attr_key == "weight":
            query += "SET r.{} = {}\n".format(attr_key, attr_val)
        else:
            query += "SET r.{} = \"{}\"\n".format(attr_key, attr_val)
    return query


def add_gene_node(node: dict) -> str:
    query = f"CREATE (n:GENE) SET n.id='{node['id']}'\n"
    attr_gene = {"genomic_type": str, "is_fragment": int}
    for attr_key, attr_val in node["attr"].items():
        if attr_key not in attr_gene:
            raise Exception(f"Unknown gene property : {str(attr_key)} for gene {node['id']}")
        if not isinstance(attr_val, attr_gene[attr_key]):
            raise Exception(f"Wrong format value for property {attr_key} in gene {node['id']}")
        query += f"SET n.{attr_key} = '{attr_val}'\n"
    return query


def add_in_family_edge(edge: dict, pangenome_data: dict):
    family_label = pangenome_data["graph"]["node_types"][edge["to"]]
    query = f"MATCH (s:GENE {{id: '{edge['from']}'}}), " \
            f"(t:{':'.join(family_label)} {{id: '{edge['to']}'}}) \n" \
            "MERGE (s)-[r:IN_FAMILY]->(t)\n"
    attr_in_fam = {}
    for attr_key, attr_val in edge["attr"].items():
        if attr_key not in attr_in_fam:
            raise Exception("Unknown attribute : " + str(attr_key))
        else:
            if not isinstance(attr_val, attr_in_fam[attr_key]):
                raise Exception(f"Wrong format value for property {attr_key} in neighbor edges")
        query += "SET r.{} = \"{}\"\n".format(attr_key, attr_val)
    return query

def add_in_contig_edge(edge: dict):
    query = f"MATCH (s:GENE {{id: '{edge['from']}'}}), " \
            f"(t:CONTIG {{id: '{edge['to']}'}}) \n" \
            "MERGE (s)-[r:IN_CONTIG]->(t)\n"
    attr_in_contig = {}
    for attr_key, attr_val in edge["attr"].items():
        if attr_key not in attr_in_contig:
            raise Exception("Unknown attribute : " + str(attr_key))
        else:
            if not isinstance(attr_val, attr_in_contig[attr_key]):
                raise Exception(f"Wrong format value for property {attr_key} in neighbor edges")
        query += "SET r.{} = \"{}\"\n".format(attr_key, attr_val)
    return query

def add_contig_node(node: dict) -> str:
    query = f"CREATE (n:CONTIG) SET n.id='{node['id']}'\n"
    attr_contig = {"is_circular": int}
    for attr_key, attr_val in node["attr"].items():
        if attr_key not in attr_contig:
            raise Exception(f"Unknown gene property : {str(attr_key)} for contig {node['id']}")
        if not isinstance(attr_val, attr_contig[attr_key]):
            raise Exception(f"Wrong format value for property {attr_key} in contig {node['id']}")
        query += f"SET n.{attr_key} = '{attr_val}'\n"
    return query


def add_genome_node(node: dict) -> str:
    query = f"CREATE (n:GENOME) SET n.id='{node['id']}'\n"
    attr_genome = {}
    for attr_key, attr_val in node["attr"].items():
        if attr_key not in attr_genome:
            raise Exception(f"Unknown gene property : {str(attr_key)} for contig {node['id']}")
        if not isinstance(attr_val, attr_genome[attr_key]):
            raise Exception(f"Wrong format value for property {attr_key} in contig {node['id']}")
        query += f"SET n.{attr_key} = '{attr_val}'\n"
    return query


def add_in_genome_edge(edge: dict):
    query = f"MATCH (s:CONTIG {{id: '{edge['from']}'}}), " \
            f"(t:GENOME {{id: '{edge['to']}'}}) \n" \
            "MERGE (s)-[r:IN_GENOME]->(t)\n"
    attr_in_genome = {}
    for attr_key, attr_val in edge["attr"].items():
        if attr_key not in attr_in_genome:
            raise Exception("Unknown attribute : " + str(attr_key))
        else:
            if not isinstance(attr_val, attr_in_genome[attr_key]):
                raise Exception(f"Wrong format value for property {attr_key} in neighbor edges")
        query += "SET r.{} = \"{}\"\n".format(attr_key, attr_val)
    return query


def get_queries(pangenome_data: dict):
    """Load the dataset from json."""

    # If the dataset is too large and everything is slow, uncomment the following code,

    # it will remove N (by default 700) random nodes from the data)

    # (do not forget to uncomment 'import numpy' on the top)

    # N = 700 # number of nodes to remove

    # nodes_to_remove = numpy.random.choice([n["id"] for n in json_data["graph"]["nodes"]], N, replace=False)

    # Create nodes
    query_list = []
    for node in tqdm(pangenome_data["graph"]["nodes"], unit="node", desc="Read nodes and create query"):
        node_type = pangenome_data["graph"]["node_types"][node["id"]]
        primary_type = node_type[0]
        # if n["id"] not in nodes_to_remove:
        if primary_type == "PANGENOME":
            query_list.append(add_pangenome_node(node))
        if primary_type == "FAMILY":
            query_list.append(add_family_node(node, node_type))
        if primary_type == "GENE":
            query_list.append(add_gene_node(node))
        if primary_type == "CONTIG":
            query_list.append(add_contig_node(node))
        if primary_type == "GENOME":
            query_list.append(add_genome_node(node))

    for edge in tqdm(pangenome_data["graph"]["edges"], unit="edges", desc="Read edges and create query"):
        primary_type = edge["type"][0]
        if primary_type == "IN_PANGENOME":
            query_list.append(add_in_pangenome_edge(edge, pangenome_data))
        if primary_type == "NEIGHBOR_OF":
            query_list.append(add_neighbor_of_edge(edge, pangenome_data))
        if primary_type == "IN_FAMILY":
            query_list.append(add_in_family_edge(edge, pangenome_data))
        if primary_type == "IN_CONTIG":
            query_list.append(add_in_contig_edge(edge))
        if primary_type == "IN_GENOME":
            query_list.append(add_in_genome_edge(edge))
    return query_list


def load_data(db: Neo4jDB, pangenome_data: dict, threads: int = 1):
    logging.getLogger().debug("Launch Neo4J db")
    queries = get_queries(pangenome_data)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        list(tqdm(executor.map(db.run, queries), total=len(queries), unit='query'))


def launch(args):
    """
    Launch functions to read systems

    :param args: Argument given
    """
    logging.getLogger().info("Read JSON")
    with open(args.json, "r+") as f:
        json_data = json.load(f)

    neo4j_db = Neo4jDB(uri=args.uri,
                       user=args.user,
                       pwd=args.pwd)
    if args.clean:
        logging.getLogger().warning("Clean Neo4J db")
        neo4j_db.clean()
    logging.getLogger().info("Begin load into Neo4J DB")
    load_data(db=neo4j_db, pangenome_data=json_data, threads=args.threads)
    neo4j_db.close()
    logging.getLogger().info("Load finished and connection to Neo4J DB closed")


def subparser(sub_parser) -> argparse.ArgumentParser:
    """
    Subparser to launch PANORAMA in Command line

    :param sub_parser : sub_parser for align command

    :return : parser arguments for align command
    """
    parser = sub_parser.add_parser("dbGraph", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_export(parser)
    return parser


def parser_export(parser):
    """
    Parser for specific argument of annot command

    :param parser: parser for annot argument
    """
    required = parser.add_argument_group(title="Required arguments",
                                         description="All of the following arguments are required :")
    required.add_argument('--json', required=True, type=str, nargs='?',
                          help="JSON file to export in Neo4J")
    required.add_argument('--uri', required=True, type=str, nargs='?',
                          help='URI to connect to the database')
    required.add_argument('--user', required=True, type=str, nargs='?',
                          help='user name to connect to the database')
    required.add_argument('--pwd', required=True, type=str, nargs='?',
                          help='password to connect to the database')
    optional = parser.add_argument_group(title="Optional arguments")
    optional.add_argument("--clean", required=False, action='store_true', default=False,
                          help="Clean the database before to add pangenomes. "
                               "WARNING, all the pangenomes will be removed.")
    optional.add_argument("--threads", required=False, nargs='?', type=int, default=1,
                          help="Number of av available threads")


if __name__ == "__main__":
    from panorama.utils import check_log, set_verbosity_level

    main_parser = argparse.ArgumentParser(description="Comparative Pangenomic analyses toolsbox",
                                          formatter_class=argparse.RawTextHelpFormatter)
    parser_export(main_parser)
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
