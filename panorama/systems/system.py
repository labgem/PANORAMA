#!/usr/bin/env python3
# coding: utf8

"""
This module provides System class which represent a biological system
"""

# default libraries
from __future__ import annotations

from typing import Dict, List, Set, Tuple, Union, Generator

import numpy as np
# installed libraries
from ppanggolin.metadata import MetaFeatures
from ppanggolin.genome import Organism
from ppanggolin.region import Region

# local libraries
from panorama.systems.models import Model, FuncUnit
from panorama.geneFamily import GeneFamily
from panorama.region import Module, Spot


class System(MetaFeatures):
    """
    Represents a biological system detected in a pangenome.

    Attributes:
        ID (str): Identifier for the system.
        model (Model): Model associated with the system.
        source (str): Source of the annotation.
        canonical (set): Set of canonical systems.
    """
    _id_counter = 0

    def __init__(self, model: Model, source: str, system_id: Union[str, int] = None, units: Set[SystemUnit] = None):
        """
        Initializes the system with given parameters.

        Args:
            system_id (Union[str, int]): Identifier for the system.
            model (Model): Model associated with the system.
            source (str): Source of the annotation.
            units (Set[SystemUnit]): A set of system unit
        """
        super().__init__()
        System._id_counter += 1
        self.ID = str(system_id) if system_id is not None else str(System._id_counter)
        self.model = model
        self.source = source
        self._unit_getter = {}
        self.canonical = set()
        self._fam2unit = None
        if units is not None:
            for fu in units:
                self.add_unit(fu)

    def __hash__(self) -> int:
        """
        Creates a hash value for the region.

        Returns:
            int: Hash value.
        """
        return hash(frozenset(self._unit_getter.items()))

    def __repr__(self):
        """
        Returns a string representation of the system.

        Returns:
            str: String representation.
        """
        return f"System ID: {self.ID}, Name: {self.name}"

    def __len__(self):
        """
        Returns the number of gene families in the system.

        Returns:
            int: Number of gene families.
        """
        return len(self._unit_getter)

    def __setitem__(self, name: str, unit: SystemUnit):
        """
        Sets a gene family in the system.

        Args:
            name: Name of the system unit.
            unit: System unit belonging to the system.

        Raises:
            TypeError: If the family is not an instance of GeneFamily.
            KeyError: If another family with the same name already exists in the system.
        """
        if name in self._unit_getter and self[name] != SystemUnit:
            raise KeyError("A different gene family with the same name already exist in the system")
        self._unit_getter[name] = unit

    def __getitem__(self, name: str) -> SystemUnit:
        """
        Gets the system unit for the given name in the system.

        Args:
            name: Name of the unit

        Returns:
            SystemUnit: unit with the given name.

        Raises:
            KeyError: If the unit with the given name does not exist in the system.
        """
        try:
            return self._unit_getter[name]
        except KeyError:
            raise KeyError(f"There isn't any unit with the name {name} in the system")

    def __delitem__(self, name: str):
        """
        Remove the gene family for the given name in the system.

        Args:
            name: Name of the unit

        Raises:
            KeyError: If the unit with the given name does not exist in the system.
        """
        try:
            self._unit_getter.pop(name)
        except KeyError:
            raise KeyError(f"There isn't any unit with the name {name} in the system")

    def __eq__(self, other: System) -> bool:
        """
        Tests whether two system objects have the same units.

        Args:
            other: Another system to test equality.

        Returns:
            bool: True if equal, False otherwise.

        Raises:
            TypeError: If trying to compare a system with another type of object.
        """
        if not isinstance(other, System):
            raise TypeError(f"Another system is expected to be compared to the first one. You gave a {type(other)}")
        return set(self._unit_getter.items()) == set(other._unit_getter.items())

    @property
    def name(self) -> str:
        """
        Name of the system inherited by the model.

        Returns:
            str: Name of the system.
        """
        return self.model.name

    @property
    def units(self) -> Generator[SystemUnit, None, None]:
        """
        Retrieves the units in the system.

        Yields:
            SystemUnit: Generator of unit.
        """
        yield from self._unit_getter.values()

    def add_unit(self, unit: SystemUnit):
        """
        Adds a system unit to the system.

        Args:
            unit: system unit to be added.

        Raises:
            AssertionError: If the functional unit is not an instance of SystemUnit.
        """
        assert isinstance(unit, SystemUnit), "FuncUnit object is expected"
        if unit.name in self._unit_getter:
            if self.get_unit(unit.name).is_superset(unit):
                # The new unit is already in the system
                pass
            elif self.get_unit(unit.name).is_subset(unit):
                del self[unit.name]
                self.add_unit(unit)
        else:
            self[unit.name] = unit
        unit.system = self

    def get_unit(self, name: str) -> SystemUnit:
        """
        Get a system unit by his name

        Args:
            name: Name of the unit.

        Returns:
            The system unit with the given name.
        """
        return self[name]

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves the families in the system.

        Yields:
            GeneFamily: Generator of gene families.
        """
        families = set()
        for unit in self.units:
            families |= set(unit.families)
        yield from families

    @property
    def number_of_families(self):
        """
        Get the number of families in the system

        Returns:
            int: Number of families.
        """
        return sum(len(unit) for unit in self.units)

    def is_superset(self, other: System) -> bool:
        """Checks if the current System includes another System.

        Args:
            other (System): The other System to check inclusion for.

        Returns:
            bool: True if all units in 'other' are included in 'self',
                  False otherwise.
        Raises:
            TypeError: If trying to compare a system with another type of object.
        """
        if not isinstance(other, System):
            raise TypeError(f"Another system is expected to be compared to the first one. You gave a {type(other)}")

        for self_unit in self.units:
            is_superset = False
            for other_unit in other.units:
                if self_unit.is_superset(other_unit):
                    is_superset = True
                    break
            if not is_superset:
                return False
        return True

    def is_subset(self, other: System) -> bool:
        """
        Checks if the current System is included in another System.


        Args:
            other (System): The other System to check inclusion against.

        Returns:
            bool: True if 'self' is included in 'other', False otherwise.
        Raises:
            TypeError: If trying to compare a system with another type of object.
        """
        if not isinstance(other, System):
            raise TypeError(f"Another system is expected to be compared to the first one. You gave a {type(other)}")

        for self_unit in self.units:
            is_subset = False
            for other_unit in other.units:
                if self_unit.is_subset(other_unit):
                    is_subset = True
                    break
            if not is_subset:
                return False
        return True

    def intersection(self, other: System) -> Set[SystemUnit]:
        """Computes the intersection of two System objects.

        Args:
            other (System): The other System to intersect with.

        Returns:
            Dict[str, SystemUnit]: A dictionary with all common units as value and name of unit as keys

        Raises:
            TypeError: If trying to compare a system with another type of object.
        """
        if not isinstance(other, System):
            raise TypeError(f"Another system is expected to be compared to the first one. You gave a {type(other)}")

        intersected_unit_getter = set()
        for s_unit in self.units:
            for o_unit in other.units:
                if s_unit.is_superset(o_unit):
                    intersected_unit_getter.add(s_unit)
                elif s_unit.is_subset(o_unit):
                    intersected_unit_getter.add(o_unit)
        return intersected_unit_getter

    def merge(self, other: System):
        if not isinstance(other, System):
            raise TypeError(f"Another system is expected to be merged with the first one. You gave a {type(other)}")

        unit_names = {unit.name for unit in set(self.units).union(other.units)}
        for name in unit_names:
            if name in other._unit_getter:
                other_unit = other.get_unit(name)
                if name in self._unit_getter:
                    self.get_unit(name).merge(other_unit)
                    other_unit.system = self
                else:
                    self.add_unit(other_unit)


    @property
    def models_families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves the gene families described in the model.

        Yields:
            GeneFamily: Generator of gene families in the model.
        """
        model_families = set()
        for unit in self.units:
            model_families |= set(unit.models_families)
        yield from model_families

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves the organisms where the system families belongs.

        Yields:
            Organism: Generator of organisms.
        """
        organisms = set()
        for unit in self.units:
            organisms |= set(unit.organisms)
        yield from organisms

    @property
    def models_organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves the organisms where the system belongs, considering only families in the model.

        Yields:
            Organism: Generator of organisms in the model.
        """
        model_organisms = set()
        for unit in self.units:
            model_organisms |= set(unit.models_organisms)
        yield from model_organisms

    def canonical_models(self) -> List[str]:
        """
        Lists the canonical models.

        Returns:
            List[str]: List of canonical models.
        """
        return self.model.canonical

    def add_canonical(self, system: System):
        """
        Adds a canonical system.

        Args:
            system: Canonical system to be added.
        """
        already_in = False
        for canon in self.canonical:
            if system.is_subset(canon):
                canon.merge(system)
                already_in = True
            elif system.is_superset(canon):
                system.ID = canon.ID
                system.merge(canon)
                self.canonical.remove(canon)
                self.canonical.add(system)
                already_in = True
        if not already_in:
            self.canonical.add(system)

    def _mk_fam2unit(self):
        self._fam2unit = {}
        for unit in self.units:
            for fam in unit.families:
                self._fam2unit[fam] = unit

    def get_metainfo(self, gene_family: GeneFamily) -> Tuple[str, int]:
        """
        Retrieves metadata for a gene family.

        Args:
            gene_family: Gene family for which metadata is retrieved.

        Returns:
            Tuple[str, int]: Tuple containing annotation source and metadata identifier.
        """
        if self._fam2unit is None:
            self._mk_fam2unit()
        return self._fam2unit[gene_family].get_metainfo(gene_family)

    def annotation_sources(self) -> Set[str]:
        """
        Returns the set of annotation sources.

        Returns:
            Set[str]: Set of annotation sources.
        """
        annotation_sources = set()
        for unit in self.units:
            annotation_sources |= unit.annotation_sources()
        return annotation_sources

    @property
    def modules(self) -> Generator[Module, None, None]:
        """
        Retrieves the modules associated with the system.

        Yields:
            Module: Generator of modules.
        """
        modules = set()
        for unit in self.units:
            modules |= set(unit.modules)
        yield from modules

    def get_module(self, identifier: int) -> Module:
        """
        Retrieves a module by its identifier.

        Args:
            identifier: Identifier of the module.

        Returns:
            Module: Module with the given identifier.

        Raises:
            KeyError: If the module with the given identifier is not associated with the system.
        """
        spot = None
        for unit in self.units:
            try:
                spot = unit.get_module(identifier)
            except KeyError:
                pass
            else:
                break
        if spot is not None:
            return spot
        else:
            raise KeyError(f"Module with name {identifier} is not associated with system {self.ID}")

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """
        Retrieves the spots associated with the system.

        Yields:
            Spot: Generator of spots.
        """
        spots = set()
        for unit in self.units:
            spots |= set(unit.spots)
        yield from spots

    def get_spot(self, identifier: int) -> Spot:
        """
        Retrieves a spot by its identifier.

        Args:
            identifier: Identifier of the spot.

        Returns:
            Spot: Spot with the given identifier.

        Raises:
            KeyError: If the spot with the given identifier is not associated with the system.
        """
        spot = None
        for unit in self.units:
            try:
                spot = unit.get_spot(identifier)
            except KeyError:
                pass
            else:
                break
        if spot is not None:
            return spot
        else:
            raise KeyError(f"Spot with name {identifier} is not associated with system {self.ID}")

    @property
    def regions(self) -> Generator[Region, None, None]:
        """
        Retrieves the regions associated with the system.

        Yields:
            Region: Generator of regions.
        """
        regions = set()
        for unit in self.units:
            regions |= set(unit.regions)
        yield from regions

    def get_region(self, name: str) -> Region:
        """
        Retrieves a region by its name.

        Args:
            name: Name of the region.

        Returns:
            Region: Region with the given name.

        Raises:
            KeyError: If the region with the given name is not associated with the system.
        """
        region = None
        for unit in self.units:
            try:
                region = unit.get_region(name)
            except KeyError:
                pass
            else:
                break
        if region is not None:
            return region
        else:
            raise KeyError(f"Region with name {name} is not associated with system {self.ID}")


class SystemUnit:
    """
    Represents a functional unit of a system detected in a pangenome.

    Attributes:
        functional_unit (FunUnit): FuncUnit model associated with the system unit detected.
        source (str): Source of the functional unit.
    """
    _id_counter = 0

    def __init__(self, functional_unit: FuncUnit, source: str, gene_families: Set[GeneFamily] = None,
                 families_to_metainfo: Dict[GeneFamily, Tuple[str, int]] = None):
        """
        Initialize a functional unit of a system detected in a pangenome.

        Args:
            functional_unit (FunUnit): FuncUnit model associated with the system unit.
            source (str): Source of the functional unit.
            gene_families: Set of gene families in the system.
            families_to_metainfo: Mapping of gene families to their metadata.
        """
        SystemUnit._id_counter += 1
        self.ID = SystemUnit._id_counter
        self.functional_unit = functional_unit
        self.source = source
        self._families_getter = {}
        self._families2metainfo = {}
        self._regions_getter = {}
        self._spots_getter = {}
        self._modules_getter = {}
        self._models_families = None
        self._system = None
        if gene_families:
            for family in gene_families:
                annot_source, meta_id = families_to_metainfo.get(family, ("", 0)) if families_to_metainfo else ("", 0)
                self.add_family(family, annot_source, meta_id)

    def __hash__(self) -> int:
        """
        Return the hash value for the given unit.

        Returns:
            int: Hash value.
        """
        return hash(frozenset(self._families_getter.items()))

    def __repr__(self):
        """
        Returns a string representation of the system.

        Returns:
            str: String representation.
        """
        return f"System unit {self.ID}, name: {self.name}, model: {self.model.name}"

    def __len__(self):
        """
        Returns the number of gene families in the functional unit.

        Returns:
            int: Number of gene families.
        """
        return len(self._families_getter)

    def __setitem__(self, name: str, family: GeneFamily):
        """
        Sets a gene family in the system.

        Args:
            name: Name of the family.
            family: Gene family belonging to the system.

        Raises:
            TypeError: If the family is not an instance of GeneFamily.
            KeyError: If another family with the same name already exists in the system.
        """
        if name in self._families_getter and self[name] != family:
            raise KeyError("A different gene family with the same name already exist in the functional unit")
        self._families_getter[name] = family

    def __getitem__(self, name: str) -> GeneFamily:
        """
        Gets the gene family for the given name in the system.

        Args:
            name: Name of the gene family.

        Returns:
            GeneFamily: Gene family with the given name.

        Raises:
            KeyError: If the family with the given name does not exist in the system.
        """
        try:
            return self._families_getter[name]
        except KeyError:
            raise KeyError(f"There isn't gene family with the name {name} in the system")

    @property
    def system(self) -> Union[System, None]:
        """
        Get the system in which the unit belongs if it has been associated with a system.

        Returns:
            System: System associated with the unit.
        """
        return self._system

    @system.setter
    def system(self, system: System):
        assert isinstance(system, System), "System must be an instance of System."

        self._system = system
        # for family in system.families:
        #     family.add_system(system)

    @property
    def functional_unit(self) -> FuncUnit:
        """
        Get the Functional unit representing the system unit

        Returns:
            FuncUnit: The Functional unit from the model
        """
        return self._fu

    @functional_unit.setter
    def functional_unit(self, func_unit: FuncUnit):
        """
        Set the Functional unit representing the system unit.

        Args:
            func_unit: The Functional unit to be set

        Raises:
            TypeError: If the functional unit is not an instance of FuncUnit.
        """
        if not isinstance(func_unit, FuncUnit):
            raise TypeError("A FuncUnit instance is expected")
        self._fu = func_unit

    @property
    def model(self) -> Model:
        """
        Get the model from where coming the functional unit representing the system unit.

        Returns:

        """
        return self.functional_unit.model

    @property
    def name(self) -> str:
        """
        Name of the system unit inherited by the functional unit.

        Returns:
            str: Name of the system unit.
        """
        return self.functional_unit.name

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves the families in the system.

        Yields:
            GeneFamily: Generator of gene families.
        """
        yield from self._families_getter.values()

    def add_family(self, gene_family: GeneFamily, annotation_source: str = "", metadata_id: int = 0):
        """
        Adds a gene family to the system.

        Args:
            gene_family: Gene family to be added.
            annotation_source: Source of the annotation.
            metadata_id: Metadata identifier.

        Raises:
            AssertionError: If the gene family is not an instance of GeneFamily.
        """
        assert isinstance(gene_family, GeneFamily), "GeneFamily object is expected"
        self._families_getter[gene_family.name] = gene_family
        self._families2metainfo[gene_family] = (annotation_source, metadata_id)
        # gene_family.add_system_unit(self)

    def __eq__(self, other: SystemUnit) -> bool:
        """
        Tests whether two system objects have the same gene families.

        Args:
            other: Another system to test equality.

        Returns:
            bool: True if equal, False otherwise.

        Raises:
            TypeError: If trying to compare a unit with another type of object.
        """
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be compared to the first one. "
                            f"You gave a {type(other)}")
        return set(self.models_families) == set(other.models_families)

    def is_superset(self, other: SystemUnit):
        """Checks if the current unit includes another one.

        Args:
            other (SystemUnit): The other SysUnit to check inclusion for.

        Returns:
            bool: True if all families in 'other' are included in 'self',
                  False otherwise.
        Raises:
            TypeError: If trying to compare a unit with another type of object.
        """
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be compared to the first one. "
                            f"You gave a {type(other)}")

        if set(self.models_families).issuperset(other.models_families):
            return True
        return False

    def is_subset(self, other: SystemUnit):
        """Checks if the current unit is included in another one.

        Args:
            other (SystemUnit): The other unit to check inclusion against.

        Returns:
            bool: True if 'self' is included in 'other', False otherwise.
        Raises:
            TypeError: If trying to compare a unit with another type of object.
        """
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be compared to the first one. "
                            f"You gave a {type(other)}")

        if set(self.models_families).issubset(other.models_families):
            return True
        return False

    def intersection(self, other: SystemUnit) -> Set[GeneFamily]:
        """Computes the intersection of gene families between two units.

        Args:
            other (SystemUnit): The other unit to intersect with.

        Returns:
            Dict[str, GeneFamily]: A Dictionary of common gene families
        Raises:
            TypeError: If trying to compare a unit with another type of object.
        """
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be compared to the first one. "
                            f"You gave a {type(other)}")
        return set(other.models_families).intersection(set(self.models_families))

    def difference(self, other: SystemUnit) -> Set[GeneFamily]:
        """Computes the difference between two units (i. e. all gene families that are in this unit but not the others.)

        Args:
            other (SystemUnit): The other SystemUnit to compute the difference with.

        Returns:
            Dict[str, GeneFamily]: A Dictionary of non-common gene families
        Raises:
            TypeError: If trying to compare a unit with another type of object.
        """
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be compared to the first one. "
                            f"You gave a {type(other)}")

        return set(self.families).difference(set(other.families))

    def symmetric_difference(self, other: SystemUnit) -> Set[GeneFamily]:
        """Computes the difference between two units. (i. e. all gene families that are in exactly one of the units.)

        Args:
            other (SystemUnit): The other SystemUnit to compute the difference with.

        Returns:
            Dict[str, GeneFamily]: A Dictionary of non-common gene families
        Raises:
            TypeError: If trying to compare a unit with another type of object.
        """
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be compared to the first one. "
                            f"You gave a {type(other)}")

        return set(other.families).symmetric_difference(set(self.families))

    def merge(self, other: SystemUnit):
        if not isinstance(other, SystemUnit):
            raise TypeError(f"Another system unit is expected to be merged with. You gave a {type(other)}")

        for family in other.difference(self):
            self.add_family(family)
            # family.del_system_unit(other.ID)

    def _get_models_families(self) -> Set[GeneFamily]:
        families = set()
        for family, metainfo in self._families2metainfo.items():
            if metainfo[1] != 0:
                families.add(family)
        return families

    @property
    def models_families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves the gene families described in the model.

        Yields:
            GeneFamily: Generator of gene families in the model.
        """
        if self._models_families is None:
            self._models_families = self._get_models_families()
        yield from self._models_families

    @property
    def nb_model_families(self) -> int:
        return len(set(self.models_families))

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves the organisms where the system unit families belongs.

        Yields:
            Organism: Generator of organisms.
        """
        organisms = set()
        for family in self.families:
            organisms |= set(family.organisms)
        yield from organisms

    @property
    def nb_organisms(self) -> int:
        return len(set(self.organisms))

    @property
    def models_organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves the organisms where the system belongs, considering only families in the model.

        Yields:
            Organism: Generator of organisms in the model.
        todo: Try to use organisms bitarray
        """
        matrix = np.zeros((self.nb_organisms, self.nb_model_families))
        org2idx = {org: i for i, org in enumerate(self.organisms)}
        for j, family in enumerate(self.models_families):
            for org in family.organisms:
                matrix[org2idx[org], j] = 1
        idx2org = {i: org for org, i in org2idx.items()}

        yield from {idx2org[i] for i in np.nonzero(np.sum(matrix == 1, axis=1) >= self.functional_unit.min_total)[0]}

    def get_metainfo(self, gene_family: GeneFamily) -> Tuple[str, int]:
        """
        Retrieves metadata for a gene family.

        Args:
            gene_family: Gene family for which metadata is retrieved.

        Returns:
            Tuple[str, int]: Tuple containing annotation source and metadata identifier.
        """
        return self._families2metainfo[gene_family]

    def annotation_sources(self) -> Set[str]:
        """
        Returns the set of annotation sources.

        Returns:
            Set[str]: Set of annotation sources.
        """
        return {metainfo[0] for metainfo in self._families2metainfo.values() if metainfo[0] != ''}

    @property
    def modules(self) -> Generator[Module, None, None]:
        """
        Retrieves the modules associated with the unit.

        Yields:
            Module: Generator of modules.
        """
        if not self._modules_getter:
            self._asso_modules()
        yield from self._modules_getter.values()

    def get_module(self, identifier: int) -> Module:
        """
        Retrieves a module by its identifier.

        Args:
            identifier: Identifier of the module.

        Returns:
            Module: Module with the given identifier.

        Raises:
            KeyError: If the module with the given identifier is not associated with the unit.
        """
        try:
            return self._modules_getter[identifier]
        except KeyError:
            raise KeyError(f"Module with identifier {identifier} is not associated with unit {self.ID}")

    def add_module(self, module: Module):
        """
        Adds a module to the unit.

        Args:
            module: Module to be added.

        Raises:
            Exception: If another module with the same identifier is already associated with the unit.
        """
        try:
            mod_in = self.get_module(identifier=module.ID)
        except KeyError:
            self._modules_getter[module.ID] = module
            module.add_unit(self)
        else:
            if module != mod_in:
                raise Exception(
                    f"Another module with identifier {module.ID} is already associated with unit {self.ID}. "
                    f"This is unexpected. Please report an issue on our GitHub")

    def _asso_modules(self):
        """
        Associates modules to the unit based on the families.
        """
        for family in self.families:
            if family.module is not None:
                self.add_module(family.module)

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """
        Retrieves the spots associated with the unit.

        Yields:
            Spot: Generator of spots.
        """
        if not self._spots_getter:
            self._make_spot_getter()
        yield from self._spots_getter.values()

    def _make_spot_getter(self):
        """
        Creates the spot getter.
        """
        for region in self.regions:
            if region.spot is not None:
                self.add_spot(region.spot)

    def get_spot(self, identifier: int) -> Spot:
        """
        Retrieves a spot by its identifier.

        Args:
            identifier: Identifier of the spot.

        Returns:
            Spot: Spot with the given identifier.

        Raises:
            KeyError: If the spot with the given identifier is not associated with the unit.
        """
        try:
            return self._spots_getter[identifier]
        except KeyError:
            raise KeyError(f"Spot with identifier {identifier} is not associated with unit {self.ID}")

    def add_spot(self, spot: Spot):
        """
        Adds a spot to the unit.

        Args:
            spot: Spot to be added.

        Raises:
            Exception: If another spot with the same identifier is already associated with the unit.
        """
        try:
            spot_in = self.get_spot(identifier=spot.ID)
        except KeyError:
            self._spots_getter[spot.ID] = spot
        else:
            if spot != spot_in:
                raise Exception(f"Another spot with identifier {spot.ID} is already associated with unit {self.ID}. "
                                f"This is unexpected. Please report an issue on our GitHub")

    @property
    def regions(self) -> Generator[Region, None, None]:
        """
        Retrieves the regions associated with the unit.

        Yields:
            Region: Generator of regions.
        """
        yield from self._regions_getter.values()

    def get_region(self, name: str) -> Region:
        """
        Retrieves a region by its name.

        Args:
            name: Name of the region.

        Returns:
            Region: Region with the given name.

        Raises:
            KeyError: If the region with the given name is not associated with the unit.
        """
        try:
            return self._regions_getter[name]
        except KeyError:
            raise KeyError(f"Region with identifier {name} is not associated with unit {self.ID}")

    def add_region(self, region: Region):
        """
        Adds a region to the unit.

        Args:
            region: Region to be added.

        Raises:
            Exception: If another region with the same identifier is already associated with the unit.
        """
        try:
            region_in = self.get_region(name=region.name)
        except KeyError:
            self._regions_getter[region.name] = region
        else:
            if region != region_in:
                raise Exception(
                    f"Another region with identifier {region.name} is already associated with unit {self.ID}. "
                    f"This is unexpected. Please report an issue on our GitHub")
