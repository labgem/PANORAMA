#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Set

# installed libraries

# local libraries
from panorama.detection.models import Model
from panorama.geneFamily import GeneFamily


class System:
    """
    This represents a biological system detected in pangenome.


    :param system_id: source of the annotation
    :param gene_families: source accesion identifier
    """

    def __init__(self, system_id: int, model: Model, gene_families: Set[GeneFamily] = None):
        """Constructor Method
        """
        self.ID = system_id
        self.model = model
        self.gene_families = gene_families

    def __repr__(self):
        return f"System ID: {self.ID}, Name: {self.name}"

    def __len__(self):
        return len(self.gene_families)

    @property
    def name(self):
        return self.model.name

    @property
    def canonical(self):
        return self.model.canonical
