#!/usr/bin/env python3
# coding: utf8

# default libraries
import logging
from typing import Dict, Generator, List, Union

# installed libraries
from ppanggolin.region import Module as Mod

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.system import System

class Module(Mod):
    """
    This class represent a hotspot.
    :param module_id: identifier of the module
    :param families: Set of families which define the module
    """
    def __init__(self, module_id: int, families: set = None):
        """
        'core' are gene families that define the module.
        'associated_families' are gene families that you believe are associated to the module in some way,
        but do not define it.
        """
        super().__init__(module_id=module_id, families=families)
        self._systemsGetter = {}
        self._sys2fam = {}

    @property
    def gene_families(self) -> GeneFamily:
        return super().families

    @property
    def systems(self):
        for system in self._systemsGetter.values():
            yield system

    def get_system(self, ID: int):
        return self._systemsGetter[ID]

    def add_system(self, system: System):
        try:
            self.get_system(system.ID)
        except KeyError:
            self._systemsGetter[system.ID] = system
        else:
            if system.name != self._systemsGetter[system.ID].name:
                raise Exception("Two system with same ID but with different name are trying to be added to module."
                                "This error is unexpected. Please report on our GitHub")
            else:
                if system.gene_families != self._systemsGetter[system.ID].gene_families:
                    raise Exception("Two system with same ID and name but with different gene families are trying to be"
                                    " added to module. This error is unexpected. Please report on our GitHub")
