#!/usr/bin/env python3
# coding: utf8

# default libraries

# local libraries
from ppanggolin.pangenome import Pangenome as Pan


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
