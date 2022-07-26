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

    def add_annotation(self, source: str, annotation: dict, force: bool = False):
        """ Add annotation to gene family

        :param source: Name of database source
        :param annotation: Identifier of the annotation
        :param force:
        """
        if self.annotation.get(source) is not None and not force:
            raise Exception(f"This source is already used in {self.name}. "
                            f"If you really want to rewrite annot use '--force' option.")
        self.annotation[source] = annotation
