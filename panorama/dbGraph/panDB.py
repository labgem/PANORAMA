#!/usr/bin/env python3
# coding:utf-8
# default libraries

import logging
from multiprocessing import Lock

import py2neo
# installed librairies
from py2neo import Graph
from dict2graph import Dict2graph


def custom_pre_func(node):
    return node


def custom_post_func(node: py2neo.Node):
    if node is not None and node.__primarylabel__ == "Gene":
        del node["tmp_id"]
    return node


class PangenomeLoader:

    def __init__(self, pangenome_name: str, pangenome_data: dict, lock: Lock, batch_size: int = 1000):
        self.name = pangenome_name
        self.lock = lock
        self.data = pangenome_data
        self.batch_size = batch_size
        self._build_loader()

    def load(self, graph: Graph):
        assert self.lock is not None, "Lock not Initialized"
        try:
            with self.lock:
                logging.debug("parse")
                self.loader.parse(self.data)
                logging.debug("index")
                self.loader.create_indexes(graph)
                logging.debug("merge")
                self.loader.merge(graph)
        except Exception as error:
            raise Exception(f"Load to Neo4j failed because : {error}")

    def _build_loader(self):
        d2g = Dict2graph()
        # c.config_dict_label_override = config.JSON2GRAPH_LABELOVERRIDE
        # c.config_func_custom_relation_name_generator = DataTransformer.nameRelation
        d2g.config_dict_primarykey_generated_hashed_attrs_by_label = {
            "Pangenome": 'AllAttributes',  # Random id
            "Family": ["name"],
            "Partition": 'AllAttributes',
            # "Persistent": ["name", "subpartition"],
            # "Shell": ["name", "subpartition"],
            # "Cloud": ["name", "subpartition"],
            "Gene": 'AllAttributes',
            "Module": "InnerContent",
            "Spot": "AllContent",
            "RGP": "InnerContent",
            "Genome": "InnerContent",
            "Contig": "AllContent"
        }
        # d2g.config_dict_concat_list_attr = {"Author": {"middle": " "}}
        # c.config_str_collection_hub_label = "{LIST_MEMBER_LABEL}Collection"
        # d2g.config_list_collection_hub_extra_labels = []

        # d2g.config_graphio_batch_size = config.COMMIT_INTERVAL
        # c.config_dict_primarykey_attr_by_label = config.JSON2GRAPH_ID_ATTR
        d2g.config_str_primarykey_generated_attr_name = "hash_id"
        d2g.config_list_blocklist_collection_hubs = [
            # "PersistentCollection",
            # "ShellCollection",
            # "CloudCollection",
            "PangenomeCollection",
            "FamilyCollection",
            "PartitionCollection",
            "GeneCollection",
            "NeighborCollection",
            "ModuleCollection",
            "SpotCollection",
            "RGPCollection",
            "GenomeCollection",
            "ContigCollection",
        ]
        d2g.config_dict_node_prop_to_rel_prop = {"Family": {"weight": ["NEIGHBOR"]}} #,
                                                 # "Shell": {"weight": ["NEIGHBOR"]},
                                                 # "Cloud": {"weight": ["NEIGHBOR"]}}  # ,  "partition": ["IN_MODULE"]}}
        d2g.config_dict_primarykey_attr_by_label = {"Family": ["name"],
                                                    "Gene": ["name"],
                                                    "Partition": ["partition"]}
        d2g.config_dict_reltype_override = {"PANGENOME_HAS_FAMILY": "IS_IN_PANGENOME",
                                            "FAMILY_HAS_GENE": "IS_IN_FAMILY",
                                            "FAMILY_HAS_PARTITION": "HAS_PARTITION",
                                            "FAMILY_HAS_FAMILY": "NEIGHBOR",
                                            "MODULE_HAS_FAMILY": "IS_IN_MODULE",
                                            "SPOT_HAS_RGP": "IS_IN_SPOT",
                                            "RGP_HAS_GENE": "IS_IN_RGP",
                                            "GENOME_HAS_CONTIG": "IS_IN_GENOME",
                                            "CONTIG_HAS_GENE": "IS_IN_CONTIG"}
        d2g.config_list_blocklist_reltypes = ["PANGENOME_HAS_MODULE",
                                              "PANGENOME_HAS_RGP",
                                              "PANGENOME_HAS_SPOT",
                                              "PANGENOME_HAS_GENOME"]
        d2g.config_bool_capitalize_labels = False
        d2g.config_func_node_post_modifier = custom_post_func
        # d2g.config_func_node_pre_modifier = custom_pre_func
        d2g.config_graphio_batch_size = self.batch_size
        self.loader = d2g


if __name__ == "__main__":
    from py2neo import GraphService

    neo4j_graph = Graph(uri="bolt://localhost:7687", user="neo4j",
                        password="PANORAMA2022", name="test")

    print(neo4j_graph.name)
