#!/usr/bin/env python3
# coding: utf8

# default libraries
from __future__ import annotations
from typing import Set, Union

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

    def __init__(self, system_id: Union[str, int], model: Model, source: str, gene_families: Set[GeneFamily] = None):
        """Constructor Method
        """
        self.ID = system_id if isinstance(system_id, str) else str(system_id)
        self.model = model
        self.source = source
        self.gene_families = gene_families
        self.canonical = set()

    def __repr__(self):
        return f"System ID: {self.ID}, Name: {self.name}"

    def __len__(self):
        return len(self.gene_families)

    @property
    def name(self):
        return self.model.name

    @property
    def canonical_models(self):
        return self.model.canonical

    def add_canonical(self, system: System):
        system.ID = f"{self.ID}.{chr(97 + len(self.canonical))}"
        self.canonical.add(system)


