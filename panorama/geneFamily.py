#!/usr/bin/env python3
# coding: utf8

"""
This module provides classes to represent gene families and study them
"""

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
    Represents a single gene family. It is a node in the pangenome graph and is aware of its genes and edges.

    Attributes:
        name (str): The name of the gene family to be printed in output files.
        _hmm (HMM, optional): The HMM associated with the gene family.
        profile (optional): The profile associated with the gene family.
        optimized_profile (optional): The optimized profile for the gene family.
        _units_getter (dict): A dictionary to retrieve system units.
        _systems_getter (dict): A dictionary to retrieve systems.
        _akin (Akin, optional): The akin families associated with other pangenomes.
    """

    def __init__(self, family_id: int, name: str):
        """
        Initializes a GeneFamily instance.

        Args:
            family_id (int): The internal identifier of the gene family.
            name (str): The name of the gene family.
        """
        super().__init__(family_id, name)
        self._hmm = None
        self.profile = None
        self.optimized_profile = None
        self._units_getter = {}
        self._systems_getter = {}
        self._akin = None

    def __repr__(self) -> str:
        """
        Returns a string representation of the GeneFamily instance.

        Returns:
            str: The string representation of the GeneFamily instance.
        """
        return f"GF {self.ID}: {self.name}"

    def __hash__(self) -> int:
        """
        Returns the hash of the GeneFamily instance.

        Returns:
            int: The hash value of the GeneFamily instance.
        """
        return hash((self.name, self.ID))

    def __eq__(self, other: GeneFamily) -> bool:
        """
        Checks if two GeneFamily instances are equal based on their genes.

        Args:
            other (GeneFamily): Another GeneFamily instance to compare.

        Returns:
            bool: True if the GeneFamily instances are equal, False otherwise.

        Raises:
            TypeError: If the other object is not a GeneFamily instance.
        """
        if not isinstance(other, GeneFamily):
            raise TypeError(
                f"Expected another GeneFamily instance for comparison, but received {type(other)}"
            )
        return set(self.genes) == set(other.genes)

    def __lt__(self, other: GeneFamily) -> bool:
        if len(self) == len(other):
            return self.ID < other.ID
        return len(self) < len(other)

    def __le__(self, other):
        return len(self) <= len(other)

    def __gt__(self, other: GeneFamily) -> bool:
        if len(self) == len(other):
            return self.ID > other.ID
        return len(self) > len(other)

    def __ge__(self, other: GeneFamily) -> bool:
        return len(self) >= len(other)

    def __ne__(self, other: GeneFamily) -> bool:
        """
        Checks if two GeneFamily instances are not equal.

        Args:
            other (GeneFamily): Another GeneFamily instance to compare.

        Returns:
            bool: True if the GeneFamily instances are not equal, False otherwise.
        """
        return not self.__eq__(other)

    def _getattr_from_ppanggolin(self, family: Fam):
        """
        Copies attributes from a PPanGGOLiN GeneFamily instance to a PANORAMA GeneFamily instance.

        Args:
            family (Fam): A PPanGGOLiN GeneFamily instance.
        """
        self._edges = family.edges
        self._genePerOrg = family._genePerOrg
        self._genes_getter = family._genes_getter
        self.removed = family.removed
        self.sequence = family.sequence
        self.partition = family.partition
        self._spots = family.spots
        self._module = family.module
        self.bitarray = family.bitarray
        for meta in family.metadata:
            self.add_metadata(meta, meta.ID)

    @property
    def HMM(self) -> HMM:
        """
        Gets the HMM associated with the GeneFamily.

        Returns:
            HMM: The HMM associated with the GeneFamily.
        """
        return self._hmm

    @HMM.setter
    def HMM(self, hmm: HMM):
        """
        Sets the HMM for the GeneFamily.

        Args:
            hmm (HMM): The HMM to associate with the GeneFamily.

        Raises:
            TypeError: If the provided hmm is not an HMM instance.
        """
        if not isinstance(hmm, HMM):
            raise TypeError(
                f"Expected an HMM instance, but received {type(hmm)}"
            )
        self._hmm = hmm

    @property
    def module(self):
        """Return module belonging to the family

        Returns:
            panorama.region.Module: module belonging to the family
        """
        return self._module

    @module.setter
    def module(self, module):
        """Set the modules belonging to the family

        Args:
            module (panorama.region.Module): module belonging to the family
        """
        from panorama.region import Module
        if not isinstance(module, Module):
            raise TypeError("Module object is expected to object of Module class")
        self._module = module

    @staticmethod
    def recast(family: Fam) -> GeneFamily:
        """
        Recasts a PPanGGOLiN GeneFamily into a PANORAMA GeneFamily.

        Args:
            family (Fam): A PPanGGOLiN GeneFamily instance.

        Returns:
            GeneFamily: The recast PANORAMA GeneFamily instance.
        """
        assert isinstance(family, Fam), "Expected a PPanGGOLiN GeneFamily instance."
        panorama_fam = GeneFamily(family_id=family.ID, name=family.name)
        panorama_fam._getattr_from_ppanggolin(family)
        return panorama_fam

    # def add_system_unit(self, unit):
    #     """
    #     Adds a system unit to the GeneFamily.
    #
    #     Args:
    #         unit: The system unit to add.
    #
    #     Raises:
    #         KeyError: If a different system unit with the same ID already exists.
    #     """
    #     if unit in self._systems_getter and self.get_system_unit(unit.ID) != unit:
    #         logging.getLogger("PANORAMA").error(
    #             f"System unit {unit.ID}: {unit.name} can't be added to family "
    #             f"because the same ID is known for {self.get_system_unit(unit.ID).name} "
    #         )
    #         raise KeyError(
    #             f"A different system with the same name already exists in the gene family {self}"
    #         )
    #     self._units_getter[unit.ID] = unit
    #
    # def get_system_unit(self, identifier: int):
    #     """
    #     Gets a system unit by its identifier.
    #
    #     Args:
    #         identifier (int): The ID of the system unit to retrieve.
    #
    #     Returns:
    #         System: The system unit with the specified identifier.
    #
    #     Raises:
    #         KeyError: If the system unit with the given ID does not exist.
    #     """
    #     try:
    #         return self._units_getter[identifier]
    #     except KeyError:
    #         raise KeyError(
    #             f"No system with ID {identifier} found in the gene family."
    #         )
    #
    # def del_system_unit(self, identifier: int):
    #     try:
    #         del self._units_getter[identifier]
    #     except KeyError:
    #         raise KeyError(
    #             f"No system with ID {identifier} found in the gene family."
    #         )

    def is_multigenic(self) -> bool:
        """
        Checks whether the GeneFamily is multigenic.

        Returns:
            bool: True if the GeneFamily is multigenic, False otherwise.
        """
        return len(self) == len(set(self.organisms))

    def is_multigenic_in_org(self, organism: Organism) -> bool:
        """
        Checks whether the GeneFamily is multigenic in a specific organism.

        Args:
            organism (Organism): The organism to check.

        Returns:
            bool: True if the GeneFamily is multigenic in the organism, False otherwise.
        """
        return len(set(self.get_genes_per_org(organism))) > 1

    @property
    def akin(self) -> Akin:
        """
        Gets the akin families associated with other pangenomes.

        Returns:
            Akin: The akin families.

        Raises:
            KeyError: If no akin families are assigned.
        """
        if self._akin is None:
            logging.getLogger('PANORAMA').debug(
                f"No akin families assigned to {self.name}."
            )
        return self._akin

    @akin.setter
    def akin(self, akin: Akin):
        """
        Sets the akin families associated with other pangenomes.

        Args:
            akin (Akin): The akin families to set.

        Raises:
            TypeError: If the provided akin is not an Akin instance.
        """
        if not isinstance(akin, Akin):
            raise TypeError(f"{akin} is not an instance of Akin.")
        if self._akin is not None and self._akin != akin:
            logging.getLogger("PANORAMA").debug(
                f"Akin families are already set for {self.name}, and a different one was provided. This could be an error."
            )
        self._akin = akin


class Akin:
    """
    Represents a group of gene families that are similar across multiple pangenomes.

    Attributes:
        ID (int): The identifier of the Akin instance.
        reference (str): The reference gene family name.
        _families (dict): A dictionary of gene families in the Akin group.
    """

    def __init__(self, identifier: int, reference: GeneFamily, *gene_families: GeneFamily) -> None:
        """
        Initializes an Akin instance.

        Args:
            identifier (int): The identifier of the Akin instance.
            reference (GeneFamily): The reference GeneFamily instance.
            *gene_families (GeneFamily): Additional GeneFamily instances to add.
        """
        self.ID = identifier
        self._families = {}
        self.reference = reference.name
        self.add(reference)
        for gene_family in gene_families:
            self.add(gene_family)

    def __setitem__(self, name: str, family: GeneFamily):
        """
        Adds a GeneFamily to the Akin group.

        Args:
            name (str): The name of the GeneFamily.
            family (GeneFamily): The GeneFamily instance to add.

        Raises:
            KeyError: If the GeneFamily with the given name already exists in the Akin group.
        """
        if name in self._families:
            raise KeyError(f"Gene family '{name}' already exists in the cluster.")
        self._families[name] = family

    def __getitem__(self, name: str) -> GeneFamily:
        """
        Retrieves a GeneFamily from the Akin group by name.

        Args:
            name (str): The name of the GeneFamily to retrieve.

        Returns:
            GeneFamily: The GeneFamily instance with the specified name.

        Raises:
            KeyError: If the GeneFamily with the given name does not exist in the Akin group.
        """
        if name not in self._families:
            raise KeyError(f"No gene family '{name}' found in the cluster.")
        return self._families[name]

    def __eq__(self, other: Akin) -> bool:
        if not isinstance(other, Akin):
            raise TypeError(f"{other} is not an instance of Akin.")
        return self.ID == other.ID

    def add(self, family: GeneFamily):
        """
        Adds a GeneFamily to the Akin group.

        Args:
            family (GeneFamily): The GeneFamily instance to add.

        Raises:
            AssertionError: If the provided family is not a GeneFamily instance.
        """
        assert isinstance(family, GeneFamily), "Expected a GeneFamily instance."
        self._families[family.name] = family
        family.akin = self

    def get(self, name: str) -> GeneFamily:
        """
        Retrieves a GeneFamily from the Akin group by name.

        Args:
            name (str): The name of the GeneFamily to retrieve.

        Returns:
            GeneFamily: The GeneFamily instance with the specified name.
        """
        return self._families[name]
