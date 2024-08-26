#!/usr/bin/env python3
# coding: utf8


# default libraries
from __future__ import annotations
import logging

# installed libraries
from ppanggolin.geneFamily import GeneFamily as Fam
from ppanggolin.genome import Organism
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
        self._units_getter = {}
        self._systems_getter = {}
        self._akin = None

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
            raise TypeError(f"Another gene family is expected to be compared to the first one. "
                            f"You give a {type(other)}")
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

        Todo look at this function is still needed
        """
        assert isinstance(family, Fam), "family must be a Gene Family object from PPanGGOLiN"
        panorama_fam = GeneFamily(family_id=family.ID, name=family.name)
        panorama_fam._getattr_from_ppanggolin(family)
        return panorama_fam

    def add_system_unit(self, unit):
        """Add a system to the family
        :param unit: System to add
        :type unit: System
        """
        if unit in self._systems_getter and self.get_system(unit.ID) != unit:
            logging.getLogger("PANORAMA").error(f"System unit {unit.ID}: {unit.name} can't be added to family "
                                                f"because same ID is known for {self.get_system(unit.ID).name} ")
            raise KeyError(f"A different system with the same name already exist in the gene family {self}")
        self._units_getter[unit.ID] = unit

    def get_system_unit(self, identifier: int):
        """Get a system by its identifier
        :param identifier: name of the system

        :return: the system searched in the family
        :rtype: System

        :raises KeyError: System with the given name does not exist in the module
        """
        try:
            return self._units_getter[identifier]
        except KeyError:
            raise KeyError(f"There isn't system with the ID {identifier} in the gene family")

    def add_system(self, system):
        """Add a system to the family
        :param system: System to add
        :type system: System
        """
        # if system.ID in self._systems_getter and self.get_system(system.ID) != system:
        #     logging.getLogger("PANORAMA").error(f"System {system.ID}: {system.name} can't be added to family "
        #                                         f"because same ID is known for {self.get_system(system.ID).name} ")
        #     raise KeyError(f"A different system with the same name already exist in the gene family {self}")
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

    def is_multigenic(self) -> bool:
        """Check whether the gene family is multigenic"""
        if len(self) == len(set(self.organisms)):
            return True
        else:
            return False

    def is_multigenic_in_org(self, organism: Organism) -> bool:
        """Check if the gene family is multigenic for this organism"""
        if len(set(self.get_genes_per_org(organism))) > 1:
            return True
        else:
            return False

    @property
    def akin(self) -> Akin:
        """Get the akin families with other pangenomes

        Returns:
            Similar families
        """
        if self._akin is None:
            logging.getLogger('PANORAMA').debug(f"Not any akin families has been assigned to {self.name}")
        return self._akin

    @akin.setter
    def akin(self, akin: Akin):
        """Set the akin families with other pangenomes

        Args:
            akin: Similar families

        Raises:
            KeyError: if akin is not instance Akin
        """
        if not isinstance(akin, Akin):
            raise TypeError(f"{akin} is not an instance of Akin.")
        if self._akin is not None and self._akin != akin:
            logging.getLogger("PANORAMA").debug(f"Akin families is already set for {self.name} and "
                                                "a different one is given. Could be an error")
        self._akin = akin


class Akin:
    """
    This class represents a group of gene families that are similar between multiple pangenomes
    """

    def __init__(self, identifier: int, reference: GeneFamily, *gene_families: GeneFamily) -> None:
        self.ID = identifier
        self._families = {}
        self.reference = reference.name
        self.add(reference)
        for gene_family in gene_families:
            self.add(gene_family)

    def __setitem__(self, name: str, family: GeneFamily):
        try:
            _ = self._families[name]
        except KeyError:
            self._families[name] = family
        else:
            raise KeyError(f"Gene family: {name} already exists in the cluster")

    def __getitem__(self, name: str):
        try:
            return self._families[name]
        except KeyError:
            raise KeyError(f"There is no gene family: {name} in the cluster")

    def add(self, family: GeneFamily):
        """
        Add a gene family to the set of akin families

        Args:
            family: the akin gene families
        """
        assert isinstance(family, GeneFamily), "A GeneFamily object is expected to be added to cluster"
        self._families[family.name] = family
        family.akin = self

    def get(self, name: str) -> GeneFamily:
        """
        Get a gene family from the set of akin families

        Args:
            name: the gene family name to get
        """
        return self._families[name]
