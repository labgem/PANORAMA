#!/usr/bin/env python3
# coding: utf8


# default libraries
from __future__ import annotations
import logging

# installed libraries
from ppanggolin.geneFamily import GeneFamily as Fam
from ppanggolin.edge import Edge
from pyhmmer.plan7 import HMM

# local libraries



class GeneFamily(Fam):
    """
    This represents a single gene family. It will be a node in the pangenome graph, and be aware of its genes and edges.

    :param family_id: The internal identifier to give to the gene family
    :param name: The name of the gene family (to be printed in output files)
    """

    def __init__(self, family_id: int, name: str):
        super().__init__(family_id, name)
        self._hmm = None
        self.profile = None
        self.optimized_profile = None
        self._systems_getter = {}

    # def __init__(self, family_id: int = None, name: str = None, family: Fam = None):
        # if family is not None:  # Allow to cast
        #     if not isinstance(family, Fam):
        #         raise TypeError("PPanGGOLiN Gene Family type is expected")
        #     if family_id is not None or name is not None:
        #         logging.getLogger("Panorama").warning("The provided gene family identifier and name will not be used")
        #     super().__init__(family.ID, family.name)
        #     self._edges = family.edges
        #     self._genePerOrg = family._genePerOrg
        #     self._genes_getter = family._genes_getter
        #     self.removed = family.removed  # for the repeated family not added in the main graph
        #     self.sequence = family.sequence
        #     self.partition = family.partition
        #     self._spots = family.spots
        #     self._modules = family.modules
        #     self.bitarray = family.bitarray
        #     for meta in family.metadata:
        #         self.add_metadata(meta.source, meta)
        #
        # elif family_id is not None and name is not None:
        #     super().__init__(family_id, name)
        # else:
        #     raise AssertionError("You must provide a Gene Family or a name and a family ID")
        # self._hmm = None
        # self.profile = None
        # self.optimized_profile = None

    def __repr__(self):
        return f"GF {self.ID}: {self.name}"

    def __hash__(self):
        return hash((self.name, self.ID))

    def __eq__(self, other: GeneFamily) -> bool:
        """
        Test whether two gene families have the same genes

        :param other: Another system to test equality

        :return: Equal or not

        :raises TypeError: Try to compare a systems with another type object
        """
        if not isinstance(other, GeneFamily):
            raise TypeError(f"Another gene family is expected to be compared to the first one. You give a {type(other)}")
        return set(self.genes) == set(other.genes)

    def __ne__(self, other: GeneFamily) -> bool:
        return not self.__eq__(other)

    def _getattr_from_ppanggolin(self, family: Fam):
        """Get attribute from ppanggolin gene family to set in PANORAMA Gene Family object
        :param family: ppanggolin gene family object
        """
        self._edges = family.edges
        self._genePerOrg = family._genePerOrg
        self._genes_getter = family._genes_getter
        self.removed = family.removed  # for the repeated family not added in the main graph
        self.sequence = family.sequence
        self.partition = family.partition
        self._spots = family.spots
        self._module = family.module
        self.bitarray = family.bitarray
        for meta in family.metadata:
            self.add_metadata(meta.source, meta)

    @property
    def HMM(self) -> HMM:
        """Return gf HMM"""
        return self._hmm

    @HMM.setter
    def HMM(self, hmm: HMM):
        if not isinstance(hmm, HMM):
            raise TypeError(f"Expected type is {HMM.__class__.name}, found type was {type(hmm)}")
        self._hmm = hmm

    @staticmethod
    def recast(family: Fam) -> GeneFamily:
        """Recast a family from PPanGGOLiN into PANORAMA Gene Family

        """
        assert isinstance(family, Fam), "family must be a Gene Family object from PPanGGOLiN"
        panorama_fam = GeneFamily(family_id=family.ID, name=family.name)
        panorama_fam._getattr_from_ppanggolin(family)
        return panorama_fam

    def add_system(self, system):
        """Add a system to the family
        :param system: System to add
        :type system: System
        """
        if system.ID in self._systems_getter and self.get_system(system.ID) != system:
            print(self.get_system(system.ID), system)
            raise KeyError("A different system with the same name already exist in the gene family")
        self._systems_getter[system.ID] = system

    def get_system(self, identifier: int):
        """Get a system by its identifier
        :param identifier: name of the system

        :return: the system searched in the family
        :rtype: System

        :raises KeyError: System with the given name does not exist in the module
        """
        try:
            return self._systems_getter[identifier]
        except KeyError:
            raise KeyError(f"There isn't system with the ID {identifier} in the gene family")
