#!/usr/bin/env python3
# coding: utf8

"""
This module defines classes for representing biological systems detected in pangenomes, including System, SystemUnit, and ClusterSystems.
 It provides methods for managing system units, gene families, modules, spots, and regions, as well as utilities for comparing,
 merging, and clustering systems across multiple pangenomes.
"""

# default libraries
from __future__ import annotations

from typing import Dict, List, Set, Tuple, Union, Generator

import numpy as np
from functools import wraps

# installed libraries
from ppanggolin.metadata import MetaFeatures
from ppanggolin.genome import Organism
from ppanggolin.region import Region

# local libraries
from panorama.systems.models import Model, FuncUnit
from panorama.geneFamily import GeneFamily
from panorama.region import Module, Spot


def check_instance_of_system_unit(method):
    """
    Decorator to ensure that a provided argument is an instance of SystemUnit.

    Args:
        method (Callable): The method to be wrapped in type-checking functionality.

    Returns:
        Callable: The wrapped method with type-checking functionality added.

    Raises:
        TypeError: If the `other` argument passed to the wrapped method is not an
            instance of SystemUnit.
    """
    @wraps(method)
    def wrapper(self, other):
        if not isinstance(other, SystemUnit):
            raise TypeError(
                f"Another system unit is expected to be compared to the first one. "
                f"You gave a {type(other)}"
            )
        return method(self, other)

    return wrapper


class SystemUnit(MetaFeatures):
    """
    Represents a functional unit detected in a pangenome system, associating a FuncUnit model with gene families,
    annotation sources, and metadata. Provides methods for managing gene families, modules, spots, and regions,
    as well as set operations and metadata retrieval for gene families within the unit.

    Attributes:
        functional_unit (FunUnit): FuncUnit model associated with the system unit detected.
        source (str): Source of the functional unit.
    """

    _id_counter = 0

    def __init__(
        self,
        functional_unit: FuncUnit,
        source: str,
        gene_families: Set[GeneFamily] = None,
        families_to_metainfo: Dict[GeneFamily, Tuple[str, int]] = None,
    ):
        """
        Initialize a SystemUnit instance representing a functional unit detected in a pangenome.

        Args:
            functional_unit (FuncUnit): The FuncUnit model associated with this system unit.
            source (str): Source of the functional unit.
            gene_families (Set[GeneFamily], optional): Set of gene families in the system unit.
            families_to_metainfo (Dict[GeneFamily, Tuple[str, int]], optional): Mapping of gene families to their annotation source and metadata ID.
        """
        super().__init__()
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
                annot_source, meta_id = (
                    families_to_metainfo.get(family, ("", 0))
                    if families_to_metainfo
                    else ("", 0)
                )
                self.add_family(family, annot_source, meta_id)

    def __hash__(self) -> int:
        """
        Compute a hash value for the SystemUnit based on its gene families.

        Returns:
            int: The hash value representing the current set of gene families in the unit.
        """
        return hash(frozenset(self._families_getter.items()))

    def __repr__(self):
        """
        Return a string representation of the SystemUnit, including its ID, name, and associated model name.

        Returns:
            str: String representation of the SystemUnit.
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
        Assigns a GeneFamily to the system by name.

        Args:
            name (str): The name of the gene family.
            family (GeneFamily): The GeneFamily instance to assign.

        Raises:
            TypeError: If the provided family is not a GeneFamily instance.
            KeyError: If a different GeneFamily with the same name already exists.
        """
        if not isinstance(family, GeneFamily):
            raise TypeError(
                f"A GeneFamily object is expected. You provided a {type(family)}."
            )
        if name in self._families_getter and self[name] != family:
            raise KeyError(
                "A different gene family with the same name already exist in the functional unit"
            )
        self._families_getter[name] = family

    def __getitem__(self, name: str) -> GeneFamily:
        """
        Retrieve a GeneFamily by its name from the system unit.

        Args:
            name (str): The name of the gene family to retrieve.

        Returns:
            GeneFamily: The gene family associated with the given name.

        Raises:
            KeyError: If no gene family with the specified name exists in the system unit.
        """
        try:
            return self._families_getter[name]
        except KeyError:
            raise KeyError(
                f"There isn't gene family with the name {name} in the system"
            )

    @property
    def system(self) -> Union[System, None]:
        """
        Return the System instance this unit is associated with, or None if not assigned.

        Returns:
            System or None: The associated System object, or None if the unit is not linked to any system.
        """
        return self._system

    @system.setter
    def system(self, system: System):
        """
        Set the System instance this unit is associated with.
        Args:
            system (System): The System instance to associate with this unit.

        Raises:
            AssertionError: The given system is not an object of class System.
        """
        if not isinstance(system, System):
            raise TypeError(f"System must be an instance of {System.__name__}.")
        self._system = system

    @property
    def functional_unit(self) -> FuncUnit:
        """
        Return the Functional unit model associated with this system unit.

        Returns:
            FuncUnit: The functional unit model.

        Raises:
            AttributeError: If the functional unit is not set.
        """
        if not hasattr(self, "_fu"):
            raise AttributeError("functional_unit is not set.")
        return self._fu

    @functional_unit.setter
    def functional_unit(self, func_unit: FuncUnit):
        """
        Set the Functional unit model associated with this system unit.

        Args:
            func_unit (FuncUnit): The Functional unit to be set

        Raises:
            TypeError: If the functional unit is not an instance of FuncUnit.
        """
        if not isinstance(func_unit, FuncUnit):
            raise TypeError(
                f"A FuncUnit instance is expected. You provided a {type(func_unit)}."
            )
        self._fu = func_unit

    @property
    def model(self) -> Model:
        """
        Return the Model instance associated with the functional unit of this system unit.

        Returns:
            Model: The model in which the functional unit is defined.
        """
        return self.functional_unit.model

    @property
    def name(self) -> str:
        """
        Returns the name of the system unit, as defined by its associated functional unit.

        Returns:
            str: Name of the system unit.
        """
        return self.functional_unit.name

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """
        Returns a generator yielding all gene families associated with this system unit.

        Yields:
            GeneFamily: Each gene family in the unit.
        """
        yield from self._families_getter.values()

    def add_family(
        self, gene_family: GeneFamily, annotation_source: str = "", metadata_id: int = 0
    ):
        """
        Adds a GeneFamily to the system and associates it with annotation source and metadata ID.

        Args:
            gene_family (GeneFamily): The gene family to add.
            annotation_source (str, optional): The source of the annotation. Defaults to "".
            metadata_id (int, optional): The metadata identifier. Defaults to 0.

        Raises:
            AssertionError: If gene_family is not an instance of GeneFamily.
        """
        assert isinstance(gene_family, GeneFamily), "GeneFamily object is expected"
        self._models_families = None  # New family added need to reset model families
        self._families_getter[gene_family.name] = gene_family
        self._families2metainfo[gene_family] = (annotation_source, metadata_id)
        # gene_family.add_system_unit(self)

    @check_instance_of_system_unit
    def __eq__(self, other: SystemUnit) -> bool:
        """
        Determine if two SystemUnit instances are equal by comparing their sets of model gene families.

        Args:
            other (SystemUnit): The other SystemUnit to compare with.

        Returns:
            bool: True if both units have the same model gene families, False otherwise.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        return set(self.models_families) == set(other.models_families)

    @check_instance_of_system_unit
    def is_superset(self, other: SystemUnit):
        """
        Check if this SystemUnit is a superset of another, i.e., contains all model gene families of the other unit.

        Args:
            other (SystemUnit): The unit to compare against.

        Returns:
            bool: True if all model gene families in 'other' are present in this unit, False otherwise.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        return set(self.models_families).issuperset(other.models_families)

    @check_instance_of_system_unit
    def is_subset(self, other: SystemUnit):
        """
        Check if this SystemUnit is a subset of another, i.e., contains all model gene families of the other unit.

        Args:
            other (SystemUnit): The unit to compare against.

        Returns:
            bool: True if this unit is a subset of the other, False otherwise.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        return set(self.models_families).issubset(other.models_families)

    @check_instance_of_system_unit
    def intersection(self, other: SystemUnit) -> Set[GeneFamily]:
        """Return the set of model gene families common to both this SystemUnit and another.

        Args:
            other (SystemUnit): The other SystemUnit to compare with.

        Returns:
            Set[GeneFamily]: Set of gene families present in both units.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        return set(other.models_families).intersection(set(self.models_families))

    @check_instance_of_system_unit
    def difference(self, other: SystemUnit) -> Set[GeneFamily]:
        """
        Return the set of gene families present in this SystemUnit but not in another.

        Args:
            other (SystemUnit): The SystemUnit to compare against.

        Returns:
            Set[GeneFamily]: Gene families unique to this unit.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        return set(self.families).difference(set(other.families))

    @check_instance_of_system_unit
    def symmetric_difference(self, other: SystemUnit) -> Set[GeneFamily]:
        """
        Return the set of gene families that are present in exactly one of this SystemUnit or another.

        Args:
            other (SystemUnit): The SystemUnit to compare with.

        Returns:
            Set[GeneFamily]: Gene families unique to either this unit or the other.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        return set(other.families).symmetric_difference(set(self.families))

    @check_instance_of_system_unit
    def merge(self, other: SystemUnit):
        """
        Merge another SystemUnit into this one by adding gene families present in the other unit but not in this unit.

        Args:
            other (SystemUnit): The SystemUnit to merge with.

        Raises:
            TypeError: If 'other' is not a SystemUnit instance.
        """
        for family in other.difference(self):
            self.add_family(family)
            # family.del_system_unit(other.ID)

    def _get_models_families(self) -> Set[GeneFamily]:
        """
        Return the set of gene families in this SystemUnit that are associated with a nonzero metadata ID.

        Returns:
            Set[GeneFamily]: Set of gene families with a nonzero metadata identifier.
        """
        families = set()
        for family, metainfo in self._families2metainfo.items():
            if metainfo[1] != 0:
                families.add(family)
        return families

    @property
    def models_families(self) -> Generator[GeneFamily, None, None]:
        """
        Return a generator yielding all gene families in this SystemUnit that are associated with a nonzero metadata ID.

        Yields:
            GeneFamily: Each gene family described in the model.
        """
        if self._models_families is None:
            self._models_families = self._get_models_families()
        yield from self._models_families

    @property
    def nb_model_families(self) -> int:
        """
        Return the number of unique gene families in this SystemUnit that are associated with a nonzero metadata ID.

        Returns:
            int: Number of distinct model-associated gene families.
        """
        return len(set(self.models_families))

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """
        Return a generator yielding all unique Organism instances associated with the gene families in this SystemUnit.

        Yields:
            Organism: Each unique organism linked to any gene family in the unit.
        """
        organisms = set()
        for family in self.families:
            organisms |= set(family.organisms)
        yield from organisms

    @property
    def nb_organisms(self) -> int:
        """
        Return the number of unique Organism instances associated with the gene families in this SystemUnit.

        Returns:
            int: Number of distinct organisms linked to the unit.
        """
        return len(set(self.organisms))

    @property
    def models_organisms(self) -> Generator[Organism, None, None]:
        """
        Return a generator yielding all unique Organism instances present in at least `min_total` model gene families within this SystemUnit.

        Yields:
            Organism: Each organism meeting the minimum model family presence threshold.

        Note:
            Considers only gene families associated with a nonzero metadata ID (model families).
            Attempts to use a matrix approach for efficient computation.
            TODO: Try to use organisms bitarray for optimization.
        """
        matrix = np.zeros((self.nb_organisms, self.nb_model_families))
        org2idx = {org: i for i, org in enumerate(self.organisms)}
        for j, family in enumerate(self.models_families):
            for org in family.organisms:
                matrix[org2idx[org], j] = 1
        idx2org = {i: org for org, i in org2idx.items()}

        yield from {
            idx2org[i]
            for i in np.nonzero(
                np.sum(matrix == 1, axis=1) >= self.functional_unit.min_total
            )[0]
        }

    def get_metainfo(self, gene_family: GeneFamily) -> Tuple[str, int]:
        """
        Return the annotation source and metadata ID associated with the given gene family.

        Args:
            gene_family (GeneFamily): The gene family for which to retrieve metadata.

        Returns:
            Tuple[str, int]: A tuple containing the annotation source and metadata identifier.
        """
        return self._families2metainfo[gene_family]

    def annotation_sources(self) -> Set[str]:
        """
        Retrieve a set of unique annotation source names from the system.

        Returns:
            Set[str]: A set containing all non-empty annotation source names.
        """
        return {
            metainfo[0]
            for metainfo in self._families2metainfo.values()
            if metainfo[0] != ""
        }

    @property
    def modules(self) -> Generator[Module, None, None]:
        """
        Generator that yields all modules associated with the unit.

        Yields:
            Module: Each module associated with the unit.
        """
        if not self._modules_getter:
            self._asso_modules()
        yield from self._modules_getter.values()

    def get_module(self, identifier: int) -> Module:
        """
        Retrieve the module associated with the given identifier.

        Args:
            identifier (int): The identifier of the module to retrieve.

        Returns:
            Module: The module corresponding to the given identifier.

        Raises:
            KeyError: If no module with the specified identifier is associated with this unit.
        """
        try:
            return self._modules_getter[identifier]
        except KeyError:
            raise KeyError(
                f"Module with identifier {identifier} is not associated with unit {self.ID}"
            )

    def add_module(self, module: Module):
        """
        Associate a module with this unit.

        Args:
            module (Module): The module to add.

        Raises:
            Exception: If a different module with the same identifier is already associated with this unit.
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
                    f"This is unexpected. Please report an issue on our GitHub"
                )

    def _asso_modules(self):
        """
        Associates modules to the unit based on the gene families present.
        """
        for family in self.families:
            if family.module is not None:
                self.add_module(family.module)

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """
        Retrieves the spots associated with the unit.

        Yields:
            Spot: Each spot associated with the unit.
        """
        if len(self._spots_getter) == 0:
            self._make_spot_getter()
        yield from self._spots_getter.values()

    def _make_spot_getter(self):
        """
        Creates and populates the spot getter with spots associated to the unit, either from regions or gene families.
        """
        if len(self._regions_getter) > 0:
            for region in self.regions:
                if region.spot is not None:
                    self.add_spot(region.spot)
        else:
            spots = set()
            for gf in self.families:
                spots |= set(gf.spots)

            for spot in spots:
                self.add_spot(spot)

    def get_spot(self, identifier: int) -> Spot:
        """
        Retrieves a spot by its identifier.

        Args:
            identifier (int): Identifier of the spot.

        Returns:
            Spot: The spot with the given identifier.

        Raises:
            KeyError: If the spot is not associated with the unit.
        """
        try:
            return self._spots_getter[identifier]
        except KeyError:
            raise KeyError(
                f"Spot with identifier {identifier} is not associated with unit {self.ID}"
            )

    def add_spot(self, spot: Spot):
        """
        Adds a spot to the unit.

        Args:
            spot (Spot): The spot to add.

        Raises:
            Exception: If a different spot with the same identifier is already associated with the unit.
        """
        try:
            spot_in = self.get_spot(identifier=spot.ID)
        except KeyError:
            self._spots_getter[spot.ID] = spot
        else:
            if spot != spot_in:
                raise Exception(
                    f"Another spot with identifier {spot.ID} is already associated with unit {self.ID}. "
                    f"This is unexpected. Please report an issue on our GitHub"
                )

    @property
    def regions(self) -> Generator[Region, None, None]:
        """
        Retrieves the regions associated with the unit.

        Yields:
            Region: Each region associated with the unit.
        """
        yield from self._regions_getter.values()

    def get_region(self, name: str) -> Region:
        """
        Retrieves a region by its name.

        Args:
            name (str): Name of the region.

        Returns:
            Region: The region with the given name.

        Raises:
            KeyError: If the region is not associated with the unit.
        """
        try:
            return self._regions_getter[name]
        except KeyError:
            raise KeyError(
                f"Region with identifier {name} is not associated with unit {self.ID}"
            )

    def add_region(self, region: Region):
        """
        Adds a region to the unit.

        Args:
            region (Region): The region to add.

        Raises:
            Exception: If a different region with the same identifier is already associated with the unit.
        """
        try:
            region_in = self.get_region(name=region.name)
        except KeyError:
            self._regions_getter[region.name] = region
        else:
            if region != region_in:
                raise Exception(
                    f"Another region with identifier {region.name} is already associated with unit {self.ID}. "
                    f"This is unexpected. Please report an issue on our GitHub"
                )


class System(MetaFeatures):
    """
    Represents a biological system detected in a pangenome, composed of one or more functional units.

    A System groups gene families that co-occur according to a defined model. It manages its internal
    system units (instances of SystemUnit), maintains references to organisms, gene families, regions,
    modules, and annotation sources, and supports set-based operations such as union, intersection,
    and inclusion with other systems.

    Attributes:
        ID (str): Unique identifier for the system. Automatically generated if not provided.
        model (Model): The model from which the system structure is derived.
        source (str): Source of the annotation or prediction (e.g., experimental, inferred).
        canonical (Set[System]): Set of canonical systems used to group variants of the same system.
        pangenome (Optional[Pangenome]): The pangenome object the system belongs to.
        cluster_id (Optional[int]): Identifier used to group homologous systems across pangenomes.
    """

    _id_counter = 0

    def __init__(
        self,
        model: Model,
        source: str,
        system_id: Union[str, int] = None,
        units: Set[SystemUnit] = None,
    ):
        """
        Initializes a System object with a model and optional functional units.

        Args:
            model (Model): The model defining the structure and required components of the system.
            source (str): Source of the system annotation (e.g., predicted, curated).
            system_id (Union[str, int], optional): Unique identifier for the system. If not provided,
                an incremental ID is assigned automatically.
            units (Set[SystemUnit], optional): Optional set of system units (functional components)
                to initialize the system with.

        Raises:
            TypeError: If any unit in `units` is not an instance of SystemUnit.
        """
        super().__init__()
        System._id_counter += 1
        self.ID = str(system_id) if system_id is not None else str(System._id_counter)
        self.model = model
        self.source = source
        self._unit_getter = {}
        self.canonical = set()
        self._fam2unit = None
        self.pangenome = None
        self.cluster_id = None
        if units is not None:
            for fu in units:
                self.add_unit(fu)

    def __hash__(self) -> int:
        """
        Computes a hash based on the set of system units.

        Returns:
            int: A hash representing the system.
        """
        return hash(frozenset(self._unit_getter.items()))

    def __repr__(self) -> str:
        """
        Returns a human-readable representation of the system.

        Returns:
            str: A summary string containing the system ID and model name.
        """
        return f"System ID: {self.ID}, Name: {self.name}"

    def __len__(self) -> int:
        """
        Returns the number of units in the system.

        Returns:
            int: Number of system units.
        """
        return len(self._unit_getter)

    def __setitem__(self, name: str, unit: SystemUnit):
        """
        Adds a system unit to the system under a given name.

        Args:
            name (str): Name under which to register the unit.
            unit (SystemUnit): The unit to register.

        Raises:
            TypeError: If `unit` is not an instance of SystemUnit.
            KeyError: If another unit with the same name already exists and differs.
        """
        if not isinstance(unit, SystemUnit):
            raise TypeError(
                f"A SystemUnit object is expected. You provided a {type(unit)}."
            )
        if name in self._unit_getter and self[name] != unit:
            raise KeyError("A different system unit with the same name already exists.")
        self._unit_getter[name] = unit

    def __getitem__(self, name: str) -> SystemUnit:
        """
        Retrieves a system unit by its name.

        Args:
            name (str): The name of the unit.

        Returns:
            SystemUnit: The system unit corresponding to the name.

        Raises:
            KeyError: If no unit with the given name exists.
        """
        try:
            return self._unit_getter[name]
        except KeyError:
            raise KeyError(f"There isn't any unit with the name {name} in the system")

    def __delitem__(self, name: str):
        """
        Removes a system unit by its name.

        Args:
            name (str): The name of the unit to remove.

        Raises:
            KeyError: If no unit with the given name exists.
        """
        try:
            self._unit_getter.pop(name)
        except KeyError:
            raise KeyError(f"There isn't any unit with the name {name} in the system")

    def __eq__(self, other: System) -> bool:
        """
        Compares this system to another for structural equality.

        Args:
            other (System): Another system to compare with.

        Returns:
            bool: True if systems are structurally identical.

        Raises:
            TypeError: If `other` is not a System.
        """
        if not isinstance(other, System):
            raise TypeError(
                f"Another system is expected to be compared to the first one. You gave a {type(other)}"
            )
        return set(self._unit_getter.items()) == set(other._unit_getter.items())

    @property
    def name(self) -> str:
        """
        Return the name of the system unit, as defined by its associated functional unit.

        Returns:
            str: The name of the system unit.
        """
        return self.model.name

    @property
    def units(self) -> Generator[SystemUnit, None, None]:
        """
        Return a generator yielding all SystemUnit instances contained in the system.

        Yields:
            SystemUnit: Each unit in the system.
        """
        yield from self._unit_getter.values()

    def add_unit(self, unit: SystemUnit):
        """
        Adds a system unit to the system, replacing it if it is a superset.

        Args:
            unit (SystemUnit): The unit to add.

        Raises:
            AssertionError: If the provided unit is not a SystemUnit.
        """
        assert isinstance(unit, SystemUnit), "SystemUnit object is expected"
        if unit.name in self._unit_getter:
            existing = self.get_unit(unit.name)
            if existing.is_superset(unit):
                return
            elif existing.is_subset(unit):
                del self[unit.name]
        self[unit.name] = unit
        unit.system = self

    def get_unit(self, name: str) -> SystemUnit:
        """
        Fetches a unit by name.

        Args:
            name (str): Name of the system unit.

        Returns:
            SystemUnit: The corresponding unit.
        """
        return self[name]

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves all gene families associated with the system.

        Returns:
            Generator[GeneFamily, None, None]: All gene families from all units.
        """
        families = set()
        for unit in self.units:
            families |= set(unit.families)
        yield from families

    @property
    def number_of_families(self):
        """
        Computes the total number of gene families in the system.

        Returns:
            int: Count of gene families across all units.
        """
        return sum(len(unit) for unit in self.units)

    @property
    def number_of_model_gene_families(self):
        """
        Computes the total number of model gene families in the system.

        Returns:
            int: Count of model gene families across all units.
        """
        return sum(unit.nb_model_families for unit in self.units)

    def is_superset(self, other: System) -> bool:
        """
        Checks if this system contains all units of another.

        Args:
            other (System): System to compare against.

        Returns:
            bool: True if self is a superset of the other.

        Raises:
            TypeError: If `other` is not a System.
        """
        if not isinstance(other, System):
            raise TypeError(
                f"Another system is expected to be compared to the first one. You gave a {type(other)}"
            )

        return all(
            any(self_unit.is_superset(other_unit) for self_unit in self.units)
            for other_unit in other.units
        )

    def is_subset(self, other: System) -> bool:
        """
        Checks if this system is fully contained in another.

        Args:
            other (System): System to compare against.

        Returns:
            bool: True if self is a subset of the other.

        Raises:
            TypeError: If `other` is not a System.
        """
        if not isinstance(other, System):
            raise TypeError(
                f"Another system is expected to be compared to the first one. You gave a {type(other)}"
            )

        return all(
            any(other_unit.is_superset(self_unit) for other_unit in other.units)
            for self_unit in self.units
        )

    def intersection(self, other: System) -> Set[SystemUnit]:
        """
        Computes the common units between this system and another.

        Args:
            other (System): Another system to intersect with.

        Returns:
            Set[SystemUnit]: Units shared between both systems.

        Raises:
            TypeError: If `other` is not a System.
        """
        if not isinstance(other, System):
            raise TypeError(
                f"Another system is expected to be compared to the first one. You gave a {type(other)}"
            )

        return {
            s_unit if s_unit.is_superset(o_unit) else o_unit
            for s_unit in self.units
            for o_unit in other.units
            if s_unit.is_superset(o_unit) or s_unit.is_subset(o_unit)
        }

    def merge(self, other: System):
        """
        Merges another system into this one by unifying their units.

        Args:
            other (System): The system to merge into this one.

        Raises:
            TypeError: If `other` is not a System.
        """
        if not isinstance(other, System):
            raise TypeError(
                f"Another system is expected to be merged with the first one. You gave a {type(other)}"
            )

        unit_names = {unit.name for unit in self.units}.union(
            {unit.name for unit in other.units}
        )
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
        Retrieves all gene families defined by the model.

        Returns:
            Generator[GeneFamily, None, None]: Model gene families.
        """
        model_families = set()
        for unit in self.units:
            model_families |= set(unit.models_families)
        yield from model_families

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves all organisms where gene families from the system were found.

        Returns:
            Generator[Organism, None, None]: Organisms represented in the system.
        """
        organisms = set()
        for unit in self.units:
            organisms |= set(unit.organisms)
        yield from organisms

    @property
    def models_organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves organisms matching model gene family requirements.

        Returns:
            Generator[Organism, None, None]: Organisms satisfying the system's model constraints.
        """
        model_organisms = set()
        for unit in self.units:
            model_organisms |= set(unit.models_organisms)
        yield from model_organisms

    def canonical_models(self) -> List[str]:
        """
        Returns the canonical model names associated with this system.

        Returns:
            List[str]: Canonical model names.
        """
        return self.model.canonical

    def add_canonical(self, system: System):
        """
        Adds a canonical system to this instance. If a similar canonical system already exists,
        merges or replaces it.

        Args:
            system (System): Canonical system to incorporate.
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
        """Internal method to create a mapping from gene family to system unit."""
        self._fam2unit = {}
        for unit in self.units:
            for fam in unit.families:
                self._fam2unit[fam] = unit

    def get_metainfo(self, gene_family: GeneFamily) -> Tuple[str, int]:
        """
        Retrieves metadata associated with a gene family (e.g., annotation source).

        Args:
            gene_family (GeneFamily): The gene family to query.

        Returns:
            Tuple[str, int]: Annotation source and metadata ID.
        """
        if self._fam2unit is None:
            self._mk_fam2unit()
        return self._fam2unit[gene_family].get_metainfo(gene_family)

    def annotation_sources(self) -> Set[str]:
        """
        Collects all annotation sources used in the system.

        Returns:
            Set[str]: Unique annotation source identifiers.
        """
        annotation_sources = set()
        for unit in self.units:
            annotation_sources |= unit.annotation_sources()
        return annotation_sources

    @property
    def modules(self) -> Generator[Module, None, None]:
        """
        Retrieves all modules associated with the system.

        Returns:
            Generator[Module, None, None]: Modules involved in the system.
        """
        modules = set()
        for unit in self.units:
            modules |= set(unit.modules)
        yield from modules

    def get_module(self, identifier: int) -> Module:
        """
        Retrieves a module by its unique identifier.

        Args:
            identifier (int): Module ID.

        Returns:
            Module: The corresponding module.

        Raises:
            KeyError: If not found.
        """
        for unit in self.units:
            try:
                return unit.get_module(identifier)
            except KeyError:
                continue
        raise KeyError(
            f"Module with ID {identifier} is not associated with system {self.ID}"
        )

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """
        Retrieves all genomic spots linked to the system.

        Returns:
            Generator[Spot, None, None]: All spots.
        """
        spots = set()
        for unit in self.units:
            spots |= set(unit.spots)
        yield from spots

    def get_spot(self, identifier: int) -> Spot:
        """
        Retrieves a spot by its unique identifier.

        Args:
            identifier (int): Spot ID.

        Returns:
            Spot: Corresponding spot.

        Raises:
            KeyError: If not found.
        """
        for unit in self.units:
            try:
                return unit.get_spot(identifier)
            except KeyError:
                continue
        raise KeyError(
            f"Spot with ID {identifier} is not associated with system {self.ID}"
        )

    @property
    def regions(self) -> Generator[Region, None, None]:
        """
        Retrieves all genomic regions associated with the system.

        Returns:
            Generator[Region, None, None]: Regions where system units were found.
        """
        regions = set()
        for unit in self.units:
            regions |= set(unit.regions)
        yield from regions

    def get_region(self, name: str) -> Region:
        """
        Retrieves a region by its name.

        Args:
            name (str): Region name.

        Returns:
            Region: Matching region.

        Raises:
            KeyError: If the region is not associated.
        """
        for unit in self.units:
            try:
                return unit.get_region(name)
            except KeyError:
                continue
        raise KeyError(
            f"Region with name '{name}' is not associated with system {self.ID}"
        )


class ClusterSystems:
    """
    Represents a cluster of systems that are considered homologous or functionally equivalent
    across different pangenomes.

    Attributes:
        ID (int): Unique identifier for the cluster.
        _systems_getter (Dict[Tuple[str, str], System]): Dictionary mapping (pangenome name, system ID)
            to System instances.
    """

    def __init__(self, identifier: int, *systems: System):
        """
        Initializes a ClusterSystems object and adds provided systems.

        Args:
            identifier (int): The identifier for the system cluster.
            *systems (System): Variable number of systems to add initially.
        """
        self.ID = identifier
        self._systems_getter = {}
        for system in systems:
            self.add(system)

    def __setitem__(self, key: Tuple[str, str], value: System):
        """
        Inserts a system into the cluster.

        Args:
            key (Tuple[str, str]): Tuple containing (pangenome name, system ID).
            value (System): The system instance to insert.

        Raises:
            KeyError: If the key already exists in the cluster.
        """
        if key in self._systems_getter:
            raise KeyError(
                f"System {key} already exists in conserved systems with ID {self.ID}"
            )
        self._systems_getter[key] = value

    def __getitem__(self, key: Tuple[str, str]) -> System:
        """
        Retrieves a system by its (pangenome name, system ID) key.

        Args:
            key (Tuple[str, str]): The lookup key.

        Returns:
            System: The system corresponding to the key.

        Raises:
            KeyError: If the key does not exist.
        """
        try:
            return self._systems_getter[key]
        except KeyError:
            raise KeyError(f"System {key} is not in conserved system with ID {self.ID}")

    def add(self, system: System) -> None:
        """
        Adds a system to the cluster.

        Args:
            system (System): The system to add.

        Raises:
            AssertionError: If the input is not a System instance.
            KeyError: If a system with the same key already exists in the cluster.
        """
        assert isinstance(
            system, System
        ), f"System object is expected, got {type(system)}"
        self[(system.pangenome.name, system.ID)] = system
        system.cluster_id = self.ID

    def get(self, system_id: str, pangenome_name: str) -> System:
        """
        Retrieves a system from the cluster using its system ID and pangenome name.

        Args:
            system_id (str): ID of the system.
            pangenome_name (str): Name of the pangenome.

        Returns:
            System: The system object.

        Raises:
            AssertionError: If system_id is not a string.
            KeyError: If the system is not found in the cluster.
        """
        assert isinstance(
            system_id, str
        ), f"System id should be a string, got {type(system_id)}"
        return self[(pangenome_name, system_id)]

    @property
    def systems(self) -> Generator[System, None, None]:
        """
        Generator yielding all systems in the cluster.

        Returns:
            Generator[System, None, None]: Systems in the cluster.
        """
        yield from self._systems_getter.values()

    def pangenomes(self) -> List[str]:
        """
        Returns the list of all pangenome names represented in the cluster.

        Returns:
            List[str]: List of pangenome names.
        """
        return [key[0] for key in self._systems_getter.keys()]
