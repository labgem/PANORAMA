#!/usr/bin/env python3
# coding:utf-8
# default libraries

import logging
from multiprocessing import Lock

# installed librairies
from py2neo import Graph
from dict2graph import Dict2graph

# local librairies
from panorama.pangenomes import Pangenome


def custom_pre_func(node):
    return node


def custom_post_func(node):
    if node is not None and node.__primarylabel__ == "Family":
        node.add_label(node["partition"])
    return node

class PangenomeLoader:

    def __init__(self, pangenome_name: str, pangenome_data: dict, lock: Lock):
        self.name = pangenome_name
        self.lock = lock
        self.data = pangenome_data
        self._build_loader()

    def load(self, graph: Graph):
        assert self.lock is not None, "Lock not Initialized"
        try:
            with self.lock:
                self.loader.parse(self.data)
                self.loader.create_indexes(graph)
                self.loader.merge(graph)
        except Exception as error:
            raise Exception(f"Load to Neo4j failed because : {error}")

    def _build_loader(self):
        d2g = Dict2graph()
        # c.config_dict_label_override = config.JSON2GRAPH_LABELOVERRIDE
        # c.config_func_custom_relation_name_generator = DataTransformer.nameRelation
        d2g.config_dict_primarykey_generated_hashed_attrs_by_label = {
            "Pangenome": None,  # Random id
            "Family": ["name", "partition", "subpartition"],
            "Gene": 'OuterContent',
            "Module": "AllContent"
        }
        # d2g.config_dict_concat_list_attr = {"Author": {"middle": " "}}
        # c.config_str_collection_hub_label = "{LIST_MEMBER_LABEL}Collection"
        # d2g.config_list_collection_hub_extra_labels = []

        # d2g.config_graphio_batch_size = config.COMMIT_INTERVAL
        # c.config_dict_primarykey_attr_by_label = config.JSON2GRAPH_ID_ATTR
        d2g.config_str_primarykey_generated_attr_name = "id"
        d2g.config_list_blocklist_collection_hubs = [
            "FamilyCollection",
            "GeneCollection",
            "NeighborCollection",
            "ModuleCollection"
        ]
        d2g.config_dict_node_prop_to_rel_prop = {"Family": {"weight": ["NEIGHBOR"], "partition": ["IN_MODULE"]}}
        d2g.config_dict_primarykey_attr_by_label = {"Family": ["name", "partition", "subpartition"],
                                                    "Module": ["module_id"]}
        d2g.config_dict_reltype_override = {"PANGENOME_HAS_FAMILY": "IN_PANGENOME",
                                            "FAMILY_HAS_GENE": "IN_FAMILY",
                                            "FAMILY_HAS_FAMILY": "NEIGHBOR",
                                            "MODULE_HAS_FAMILY": "IN_MODULE"}
        d2g.config_list_blocklist_reltypes = ["PANGENOME_HAS_MODULE"]
        d2g.config_bool_capitalize_labels = False
        # d2g.config_dict_label_override = {
        #     "location": "Location",
        #     "cite_spans": "Citation",
        #     "affiliation": "Affiliation",
        # }
        d2g.config_func_node_post_modifier = custom_post_func
        # d2g.config_func_node_pre_modifier = custom_pre_func

        # c.config_dict_property_name_override = config.JSON2GRAPH_PROPOVERRIDE
        self.loader = d2g


if __name__ == "__main__":
    from py2neo import GraphService

    neo4j_graph = Graph(uri="bolt://localhost:7687", user="neo4j",
                        password="PANORAMA2022", name="test")

    print(neo4j_graph.name)