#!/usr/bin/env python3
# coding: utf8

# default libraries

# installed libraries
from ppanggolin.geneFamily import GeneFamily as Fam
# local libraries


class GeneFamily(Fam):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.
    :param family_id: The internal identifier to give to the gene family
    :param name: The name of the gene family (to be printed in output files)
    """

    def __init__(self, family_id: int, name: str):
        super().__init__(family_id, name)
        self.annotation = {}

    def __repr__(self):
        return f"GF {self.ID}: {self.name}"

    def get_annot(self, source: str):
        if self.annotation.get(source) is not None:
            return self.annotation[source]
        else:
            return None

    def add_annotation(self, source: str, annotation: list):
        """ Add annotation to gene family

        :param source: Name of database source
        :param annotation: Identifier of the annotation
        :param force:
        """
        if self.get_annot(source) is not None:
            self.annotation[source].add(annotation)
        else:
            self.annotation[source] = {annotation}