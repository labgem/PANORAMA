#!/usr/bin/env python3
# coding: utf8

# default libraries
from __future__ import annotations
from typing import List, Set, Union, Generator

# installed libraries
from ppanggolin.metadata import MetaFeatures

# local libraries
from panorama.models import Model
from panorama.geneFamily import GeneFamily


class System(MetaFeatures):
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
        self._families_getter = {}
        self.canonical = set()
        if gene_families is not None:
            for family in gene_families:
                self.add_family(family)

    def __hash__(self) -> int:
        """Create a hash value for the region
        """
        return hash(self.name)

    def __repr__(self):
        return f"System ID: {self.ID}, Name: {self.name}"

    def __len__(self):
        return len(self._families_getter)

    def __setitem__(self, name: str, family: GeneFamily):
        """Set a gene family in the system

        :param name: Name of the family
        :param family: Gene family belonging to the system

        :raises TypeError: Family is not instance GeneFamily
        :raises KeyError: Another family with the same name already exists in the system
        """
        if name in self._families_getter and self[name] != family:
            raise KeyError("A different gene family with the same name already exist in the system")
        self._families_getter[name] = family

    def __getitem__(self, name) -> GeneFamily:
        """Get the gene family for the given name in the system

        :param name: Name of the gene family

        :return: Gene family with the given name

        :raises KeyError: Family with the given name does not exist in the system
        """
        try:
            return self._families_getter[name]
        except KeyError:
            raise KeyError(f"There isn't gene family with the name {name} in the system")

    def __eq__(self, other: System) -> bool:
        """
        Test whether two systems objects have the same gene families

        :param other: Another system to test equality

        :return: Equal or not

        :raises TypeError: Try to compare a systems with another type object
        """
        if not isinstance(other, System):
            raise TypeError(f"Another systems is expected to be compared to the first one. You give a {type(other)}")
        return set(self.families) == set(other.families)

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """Get the families in the system
        """
        for family in self._families_getter.values():
            yield family

    @property
    def name(self) -> str:
        """Name of the system inhereted by the model

        :return: name of the system
        """
        return self.model.name

    def canonical_models(self) -> List[str]:
        """List of the canonical models

        :return: canonical models
        """
        return self.model.canonical

    def add_canonical(self, system: System):
        """Add a canonical system"""

        system.ID = f"{self.ID}.{chr(97 + len(self.canonical))}"
        self.canonical.add(system)

    def add_family(self, gene_family: GeneFamily):
        assert isinstance(gene_family, GeneFamily), "GeneFamily object is expected"
        self._families_getter[gene_family.name] = gene_family

    def is_subset(self, system: System):
        """Check if another system is subset of this one
        """
        if set(self.families).issubset(set(system.families)):
            return True
        else:
            return False

    def is_superset(self, system: System):
        if set(self.families).issuperset(set(system.families)):
            return True
        else:
            return False

    def intersection(self, system: System) -> Set[GeneFamily]:
        """Intersection of two systems
        """
        return set(self.families).intersection(set(system.families))

    def difference(self, system: System) -> Set[GeneFamily]:
        """Difference of two systems

        :return: All self gene family that are not in the given system
        """
        return set(self.families).difference(set(system.families))

    def symmetric_difference(self, system: System) -> Set[GeneFamily]:
        """Symmetric difference of two systems

        Args:
            system:

        Returns: All gene family that are not commonly present in self and system
        """
        assert isinstance(system, System), "System object is expected"
        return set(self.families).symmetric_difference(set(system.families))
