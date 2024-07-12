#!/usr/bin/env python3
# coding: utf8

"""
This module provides System class which represent a biological system
"""

# default libraries
from __future__ import annotations
from typing import Dict, List, Set, Tuple, Union, Generator

# installed libraries
import networkx as nx
from ppanggolin.metadata import MetaFeatures
from ppanggolin.genome import Organism
from ppanggolin.region import Region

# local libraries
from panorama.systems.models import Model
from panorama.geneFamily import GeneFamily
from panorama.region import Module, Spot


class System(MetaFeatures):
    """
    Represents a biological system detected in a pangenome.

    Attributes:
        ID (str): Identifier for the system.
        model (Model): Model associated with the system.
        system_source (str): Source of the annotation.
        canonical (set): Set of canonical systems.
    """

    def __init__(self, system_id: Union[str, int], model: Model, source: str,
                 gene_families: Set[GeneFamily] = None, families_to_metainfo: Dict[GeneFamily, Tuple[str, int]] = None):
        """
        Initializes the system with given parameters.

        Args:
            system_id: Identifier for the system.
            model: Model associated with the system.
            source: Source of the annotation.
            gene_families: Set of gene families in the system.
            families_to_metainfo: Mapping of gene families to their metadata.
        """
        super().__init__()
        self.ID = str(system_id)
        self.model = model
        self.system_source = source
        self._families_getter = {}
        self._families2metainfo = {}
        self._models_families = None
        self._graph = None
        self._regions_getter = {}
        self._spots_getter = {}
        self._modules_getter = {}
        self.canonical = set()
        if gene_families:
            for family in gene_families:
                annot_source, meta_id = families_to_metainfo.get(family, ("", 0)) if families_to_metainfo else ("", 0)
                self.add_family(family, annot_source, meta_id)

    def __hash__(self) -> int:
        """
        Creates a hash value for the region.

        Returns:
            int: Hash value.
        """
        return hash(self.name)

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
            raise KeyError("A different gene family with the same name already exist in the system")
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

    def __eq__(self, other: System) -> bool:
        """
        Tests whether two system objects have the same gene families.

        Args:
            other: Another system to test equality.

        Returns:
            bool: True if equal, False otherwise.

        Raises:
            TypeError: If trying to compare a system with another type of object.
        """
        if not isinstance(other, System):
            raise TypeError(f"Another system is expected to be compared to the first one. You gave a {type(other)}")
        return set(self.families) == set(other.families)

    @property
    def name(self) -> str:
        """
        Name of the system inherited by the model.

        Returns:
            str: Name of the system.
        """
        return self.model.name

    @property
    def graph(self) -> nx.Graph:
        """
        Retrieves the graph associated with the system.

        Returns:
            nx.Graph: Graph object.

        Raises:
            AssertionError: If the system has not been associated with a graph.
        """
        if self._graph is None:
            raise AssertionError("System has not been associated with a graph")
        return self._graph

    @graph.setter
    def graph(self, graph: nx.Graph):
        """
        Sets the graph for the system.

        Args:
            graph: Graph to be associated with the system.

        Raises:
            AssertionError: If the graph is not an instance of nx.Graph.
        """
        assert isinstance(graph, nx.Graph), f"Graph must be an instance of nx.Graph, not {type(graph)}"
        self._graph = graph

    @property
    def source(self) -> str:
        """
        Alias to return the source name for the system.

        Returns:
            str: Source name.
        """
        return self.system_source

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves the families in the system.

        Yields:
            GeneFamily: Generator of gene families.
        """
        yield from self._families_getter.values()

    def _get_models_families(self):
        for family, metainfo in self._families2metainfo.items():
            if metainfo[1] != 0:
                yield family

    @property
    def models_families(self) -> Generator[GeneFamily, None, None]:
        """
        Retrieves the gene families described in the model.

        Yields:
            GeneFamily: Generator of gene families in the model.
        """
        if self._models_families is None:
            self._models_families = set(self._get_models_families())
        yield from self._models_families

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves the organisms where the system belongs.

        Yields:
            Organism: Generator of organisms.
        """
        organisms = set()
        for family in self.families:
            organisms |= set(family.organisms)
        yield from organisms

    @property
    def models_organisms(self) -> Generator[Organism, None, None]:
        """
        Retrieves the organisms where the system belongs, considering only families in the model.

        Yields:
            Organism: Generator of organisms in the model.
        """
        organisms = set()
        for family in self.models_families:
            organisms |= set(family.organisms)
        yield from organisms

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
                already_in = True
            elif system.is_superset(canon):
                system.ID = canon.ID
                self.canonical.remove(canon)
                self.canonical.add(system)
                already_in = True
        if not already_in:
            self.canonical.add(system)

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

    def is_subset(self, system: System) -> bool:
        """
        Checks if another system is a subset of this one.

        Args:
            system: System to be checked.

        Returns:
            bool: True if subset, False otherwise.
        """
        return set(self.families).issubset(set(system.families))

    def is_superset(self, system: System) -> bool:
        """
        Checks if another system is a superset of this one.

        Args:
            system: System to be checked.

        Returns:
            bool: True if superset, False otherwise.
        """
        return set(self.families).issuperset(set(system.families))

    def intersection(self, system: System) -> Set[GeneFamily]:
        """
        Finds the intersection of two systems.

        Args:
            system: System to intersect with.

        Returns:
            Set[GeneFamily]: Set of common gene families.
        """
        return set(self.families).intersection(set(system.families))

    def difference(self, system: System) -> Set[GeneFamily]:
        """
        Finds the difference between two systems.

        Args:
            system: System to find the difference with.

        Returns:
            Set[GeneFamily]: Set of gene families in the current system but not in the given system.
        """
        return set(self.families).difference(set(system.families))

    def symmetric_difference(self, system: System) -> Set[GeneFamily]:
        """
        Finds the symmetric difference between two systems.

        Args:
            system: System to find the symmetric difference with.

        Returns:
            Set[GeneFamily]: Set of gene families that are not commonly present in both systems.

        Raises:
            AssertionError: If the given system is not an instance of System.
        """
        assert isinstance(system, System), "System object is expected"
        return set(self.families).symmetric_difference(set(system.families))

    @property
    def modules(self) -> Generator[Module, None, None]:
        """
        Retrieves the modules associated with the system.

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
            KeyError: If the module with the given identifier is not associated with the system.
        """
        try:
            return self._modules_getter[identifier]
        except KeyError:
            raise KeyError(f"Module with identifier {identifier} is not associated with system {self.ID}")

    def add_module(self, module: Module):
        """
        Adds a module to the system.

        Args:
            module: Module to be added.

        Raises:
            Exception: If another module with the same identifier is already associated with the system.
        """
        try:
            mod_in = self.get_module(identifier=module.ID)
        except KeyError:
            self._modules_getter[module.ID] = module
            module.add_system(self)
        else:
            if module != mod_in:
                raise Exception(f"Another module with identifier {module.ID} is already associated with system {self.ID}. "
                                f"This is unexpected. Please report an issue on our GitHub")

    def _asso_modules(self):
        """
        Associates modules to the system based on the families.
        """
        for family in self.families:
            if family.module is not None:
                self.add_module(family.module)

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """
        Retrieves the spots associated with the system.

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
            KeyError: If the spot with the given identifier is not associated with the system.
        """
        try:
            return self._spots_getter[identifier]
        except KeyError:
            raise KeyError(f"Spot with identifier {identifier} is not associated with system {self.ID}")

    def add_spot(self, spot: Spot):
        """
        Adds a spot to the system.

        Args:
            spot: Spot to be added.

        Raises:
            Exception: If another spot with the same identifier is already associated with the system.
        """
        try:
            spot_in = self.get_spot(identifier=spot.ID)
        except KeyError:
            self._spots_getter[spot.ID] = spot
        else:
            if spot != spot_in:
                raise Exception(f"Another spot with identifier {spot.ID} is already associated with system {self.ID}. "
                                f"This is unexpected. Please report an issue on our GitHub")

    @property
    def regions(self) -> Generator[Region, None, None]:
        """
        Retrieves the regions associated with the system.

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
            KeyError: If the region with the given name is not associated with the system.
        """
        try:
            return self._regions_getter[name]
        except KeyError:
            raise KeyError(f"Region with identifier {name} is not associated with system {self.ID}")

    def add_region(self, region: Region):
        """
        Adds a region to the system.

        Args:
            region: Region to be added.

        Raises:
            Exception: If another region with the same identifier is already associated with the system.
        """
        try:
            region_in = self.get_region(name=region.name)
        except KeyError:
            self._regions_getter[region.name] = region
        else:
            if region != region_in:
                raise Exception(f"Another region with identifier {region.name} is already associated with system {self.ID}. "
                                f"This is unexpected. Please report an issue on our GitHub")
