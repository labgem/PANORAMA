#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Union

# install libraries
import ppanggolin.geneFamily
from ppanggolin.pangenome import Pangenome as Pan

# local libraries
from panorama.geneFamily import GeneFamily


class Pangenome(Pan):
    """
    This is a class representing pangenome based on PPanGGOLLiN class. It is used as a basic unit for all the analysis
    to access to the different elements of your pangenome, such as organisms, contigs, genes or gene families.
    This class provide some more methods needed to analyse pangenome.

    :param name: Name of the pangenome
    """

    def __init__(self, name):
        """Constructor method.
        """
        super().__init__()
        self.name = name

    def get_annnotation(self, annotation: str, source: str) -> Union[GeneFamily, list]:
        """ Get the annotation in pangenome gene family by name or source of annotation

        :param annotation: Name of the annotation
        :param source: Name of the source

        :return: Gene famillies with the annotation or source
        """
        pass

    def _create_gene_family(self, name: str) -> GeneFamily:
        """Creates a gene family object with the given `name`

        :param name: the name to give to the gene family. Must not exist already.
        :return: the created GeneFamily object
        """

        new_fam = GeneFamily(family_id=self.max_fam_id, name=name)
        self.max_fam_id += 1
        self._famGetter[new_fam.name] = new_fam
        return new_fam

    def get_gene_family(self, name: str) -> GeneFamily:
        return super(Pangenome, self).get_gene_family(name)


class Pangenomes:
    """
    This class represente a group of pangenome object.
    """

    def __init__(self):
        """Constructor method
        """
        self.pangenomes_set = {}

    def add_pangenome(self, pangenome: Pangenome):
        """ Add a pangenome object

        :param pangenome: Pangenome object
        :return:
        """
        self.pangenomes_set[pangenome.name] = pangenome
