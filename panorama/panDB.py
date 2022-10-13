#!/usr/bin/env python3
# coding:utf-8
# default libraries

import logging
from multiprocessing import Lock
from typing import List

# installed librairies
from py2neo import Graph
from graphio import NodeSet, RelationshipSet
from dict2graph import Dict2graph

# local librairies
from panorama.pangenomes import Pangenome


def custom_pre_func(node):
    return node


def custom_post_func(node):
    if node is not None and node.__primarylabel__ == "Family":
        node.add_label(node["partition"])
    return node


# class panDB:
#     node_sets: List[NodeSet] = []
#     relationship_sets: List[RelationshipSet] = []
#
#     def __init__(self):
#         self.node_sets: List[NodeSet] = []
#         self.relationship_sets: List[RelationshipSet] = []
#
#         _current_nodes = None
#         _current_rels = None
#
#     def get_nodesets(self):
#         for nodeset in self.node_sets:
#             yield nodeset
#
#     def get_relationshipsets(self):
#         for relationshipset in self.relationship_sets:
#             yield relationshipset
#
#     def create_indexes(self, graph: Graph):
#         for rels in self.relationship_sets.values():
#             rels.create_index(graph)
#         for nodes in self.node_sets.values():
#             nodes.create_index(graph)
#
#     def create(self, graph: Graph):
#         for nodes in self.nodeSets.values():
#             nodes.create(graph)
#         for rels in self.relationshipSets.values():
#             rels.create(graph)

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
            "Pangenome": 'AllAttributes',  # Random id
            "Family": 'AllContent',
            "Persistent": ["name", "subpartition"],
            "Shell": ["name", "subpartition"],
            "Cloud": ["name", "subpartition"],
            "Gene": 'AllAttributes',
            "Module": "AllContent",
            "Spot": "AllContent",
            "RGP": "InnerContent",
            "Genome": "AllContent",
            "Contig": "AllContent"
        }
        # d2g.config_dict_concat_list_attr = {"Author": {"middle": " "}}
        # c.config_str_collection_hub_label = "{LIST_MEMBER_LABEL}Collection"
        # d2g.config_list_collection_hub_extra_labels = []

        # d2g.config_graphio_batch_size = config.COMMIT_INTERVAL
        # c.config_dict_primarykey_attr_by_label = config.JSON2GRAPH_ID_ATTR
        d2g.config_str_primarykey_generated_attr_name = "hash_id"
        d2g.config_list_blocklist_collection_hubs = [
            "PersistentCollection",
            "ShellCollection",
            "CloudCollection",
            "GeneCollection",
            "NeighborCollection",
            "ModuleCollection",
            "SpotCollection",
            "RGPCollection",
            "GenomeCollection",
            "ContigCollection",
        ]
        d2g.config_dict_node_prop_to_rel_prop = {"Persistent": {"weight": ["NEIGHBOR"]},
                                                 "Shell": {"weight": ["NEIGHBOR"]},
                                                 "Cloud": {"weight": ["NEIGHBOR"]}}  # ,  "partition": ["IN_MODULE"]}}
        d2g.config_dict_primarykey_attr_by_label = {"Persistent": ["name", "subpartition"],
                                                    "Shell": ["name", "subpartition"],
                                                    "Cloud": ["name", "subpartition"]}
        d2g.config_dict_reltype_override = {"PANGENOME_HAS_FAMILY": "IS_IN_PANGENOME",
                                            "PERSISTENT_HAS_GENE": "IS_IN_FAMILY",
                                            "SHELL_HAS_GENE": "IS_IN_FAMILY",
                                            "CLOUD_HAS_GENE": "IS_IN_FAMILY",
                                            "PERSISTENT_HAS_PERSISTENT": "NEIGHBOR",
                                            "PERSISTENT_HAS_SHELL": "NEIGHBOR",
                                            "PERSISTENT_HAS_CLOUD": "NEIGHBOR",
                                            "SHELL_HAS_PERSISTENT": "NEIGHBOR",
                                            "SHELL_HAS_SHELL": "NEIGHBOR",
                                            "SHELL_HAS_CLOUD": "NEIGHBOR",
                                            "CLOUD_HAS_CLOUD": "NEIGHBOR",
                                            "CLOUD_HAS_PERSISTENT": "NEIGHBOR",
                                            "CLOUD_HAS_SHELL": "NEIGHBOR",
                                            "MODULE_HAS_PERSISTENT": "IS_IN_MODULE",
                                            "MODULE_HAS_SHELL": "IS_IN_MODULE",
                                            "MODULE_HAS_CLOUD": "IS_IN_MODULE",
                                            "SPOT_HAS_RGP": "IS_IN_SPOT",
                                            "RGP_HAS_GENE": "IS_IN_RGP",
                                            "GENOME_HAS_CONTIG": "IS_IN_GENOME",
                                            "CONTIG_HAS_GENE": "IS_IN_CONTIG"}
        d2g.config_list_blocklist_reltypes = ["PANGENOME_HAS_MODULE",
                                              "PANGENOME_HAS_RGP",
                                              "PANGENOME_HAS_SPOT",
                                              "PANGENOME_HAS_GENOME"]
        d2g.config_bool_capitalize_labels = False
        # d2g.config_dict_label_override = {
        #     "location": "Location",
        #     "cite_spans": "Citation",
        #     "affiliation": "Affiliation",
        # }
        # d2g.config_func_node_post_modifier = custom_post_func
        # d2g.config_func_node_pre_modifier = custom_pre_func

        # c.config_dict_property_name_override = config.JSON2GRAPH_PROPOVERRIDE
        self.loader = d2g


class SimilarityLoader:
    def __init__(self, data: dict, lock: Lock):
        self.lock = lock
        self.data = data
        self._build_loader()

    def load(self, graph: Graph):
        assert self.lock is not None, "Lock not Initialized"
        try:
            with self.lock:
                self.loader.parse(self.data)
                self.loader.create_indexes_rels(graph)
                self.loader.merge_rels(graph)

        except Exception as error:
            raise Exception(f"Load to Neo4j failed because : {error}")

    def _build_loader(self):
        d2g = Dict2graph()
        d2g.config_dict_primarykey_generated_hashed_attrs_by_label = {
            "Family": ["name"]
        }
        d2g.config_list_blocklist_collection_hubs = [
            "FamilyCollection"
        ]
        d2g.config_dict_create_merge_depending_scheme = {
            "merge": ["Family", "FAMILY_HAS_FAMILY"],
        }
        d2g.config_dict_node_prop_to_rel_prop = {"Family": {"identity": ["SIMILAR"], "covery": ["SIMILAR"]}}
        d2g.config_dict_primarykey_attr_by_label = {"Family": ["name"]}
        d2g.config_list_blocklist_reltypes = ["SIMILAR_HAS_FAMILY"]
        d2g.config_dict_reltype_override = {"FAMILY_HAS_FAMILY": "SIMILAR"}
        d2g.config_bool_capitalize_labels = False
        d2g.config_graphio_batch_size = 10
        # d2g.config_func_node_post_modifier = custom_post_func
        self.loader = d2g


if __name__ == "__main__":
    from py2neo import GraphService

    neo4j_graph = Graph(uri="bolt://localhost:7687", user="neo4j",
                        password="PANORAMA2022", name="test")

    print(neo4j_graph.name)
