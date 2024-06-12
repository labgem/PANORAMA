#!/usr/bin/env python3
# coding: utf8

# default libraries
from __future__ import annotations
from typing import Dict, List, Set, Union, Generator

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
    This represents a biological system detected in pangenome.


    :param system_id: source of the annotation
    :param gene_families: source accesion identifier
    """

    def __init__(self, system_id: Union[str, int], model: Model, source: str, gene_families: Set[GeneFamily] = None,
                 families_source: List[str] = None, families_to_metadata_id: Dict[GeneFamily, int] = None):
        """Constructor Method
        """
        self.ID = system_id if isinstance(system_id, str) else str(system_id)
        self.model = model
        self.system_source = source
        self.families_sources = families_source if families_source is not None else [source]
        self._families_getter = {}
        self._families2metadata_id = {}
        self._models_families = None
        self._graph = None
        self._regions_getter = {}
        self._spots_getter = {}
        self._modules_getter = {}
        self.canonical = set()
        if gene_families is not None:
            for family in gene_families:
                meta_id = 0
                if families_to_metadata_id is not None:
                    meta_id = families_to_metadata_id[family] if family in families_to_metadata_id else 0
                self.add_family(family, meta_id)

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
    def name(self) -> str:
        """Name of the system inhereted by the model

        :return: name of the system
        """
        return self.model.name

    @property
    def graph(self) -> nx.Graph:
        if self._graph is None:
            raise AssertionError("System has not be associated with a graph")
        return self._graph

    @graph.setter
    def graph(self, graph: nx.Graph):
        assert isinstance(graph, nx.Graph), f"Graph must be an instance of nx.Graph, not {type(graph)}"
        self._graph = graph

    @property
    def source(self) -> str:
        """Alias to return the source name for the system
        """
        return self.system_source

    @property
    def families(self) -> Generator[GeneFamily, None, None]:
        """Get the families in the system
        """
        for family in self._families_getter.values():
            yield family

    def _get_models_families(self):
        self._models_families = {}
        families_name = set()
        for fam in self.model.families:
            families_name.add(fam.name)
            families_name |= fam.exchangeable

        for family in self.families:
            for families_source in self.families_sources:
                metadata = family.get_metadata_by_source(families_source)
                if metadata is not None:
                    if any(value.protein_name in families_name for value in metadata):
                        yield family

    @property
    def models_families(self) -> Generator[GeneFamily, None, None]:
        """Get the gene families that are describing in the model
        """
        if self._models_families is None:
            self._models_families = set(self._get_models_families())
        yield from self._models_families

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """Get the organisms where the system belongs
        """
        organisms = set()
        for family in self.families:
            organisms |= set(family.organisms)
        yield from organisms

    @property
    def models_organisms(self) -> Generator[Organism, None, None]:
        """Get the organisms where the system belongs but take into account only families that are in the model
        """
        organisms = set()
        for family in self.models_families:
            organisms |= set(family.organisms)
        yield from organisms

    def canonical_models(self) -> List[str]:
        """List of the canonical models

        :return: canonical models
        """
        return self.model.canonical

    def add_canonical(self, system: System):
        """Add a canonical system"""
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
            system.ID = f"{self.ID}.{chr(97 + len(self.canonical))}"
            self.canonical.add(system)

    def add_family(self, gene_family: GeneFamily, metadata_id: int = 0):
        assert isinstance(gene_family, GeneFamily), "GeneFamily object is expected"
        self._families_getter[gene_family.name] = gene_family
        self._families2metadata_id[gene_family] = metadata_id

    def get_metadata_id(self, gene_family: GeneFamily):
        return self._families2metadata_id[gene_family]

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

    @property
    def modules(self) -> Generator[Module, None, None]:
        """Get the modules associated to the system
        """
        if self._modules_getter == {}:
            self._asso_modules()
        yield from self._modules_getter.values()

    def get_module(self, identifier: int) -> Module:
        try:
            module = self._modules_getter[identifier]
        except KeyError:
            raise KeyError(f"Module with identifier {identifier} is not associated to system {self.ID}")
        else:
            return module

    def add_module(self, module: Module):
        try:
            mod_in = self.get_module(identifier=module.ID)
        except KeyError:
            self._modules_getter[module.ID] = module
            module.add_system(self)
        else:
            if module != mod_in:
                raise Exception(f"Another module with identifier {module.ID} is already associated to system {self.ID}."
                                f"This is unexpected. Please report an issue on our GitHub")

    def _asso_modules(self):
        for family in self.families:
            if family.module is not None:
                self.add_module(family.module)

    @property
    def spots(self) -> Generator[Spot, None, None]:
        if len(self._spots_getter) == 0:
            self._make_spot_getter()
        yield from self._spots_getter.values()

    def _make_spot_getter(self):
        spots = set()
        for region in self.regions:
            if region.spot is not None:
                self.add_spot(region.spot)

    def get_spot(self, identifier: int):
        try:
            spot = self._spots_getter[identifier]
        except KeyError:
            raise KeyError(f"Module with identifier {identifier} is not associated to system {self.ID}")
        else:
            return spot

    def add_spot(self, spot: Spot):
        try:
            spot_in = self.get_spot(identifier=spot.ID)
        except KeyError:
            self._spots_getter[spot.ID] = spot
        else:
            if spot != spot_in:
                raise Exception(f"Another spot with identifier {spot.ID} is already associated to system {self.ID}."
                                f"This is unexpected. Please report an issue on our GitHub")

    @property
    def regions(self) -> Generator[Region, None, None]:
        yield from self._regions_getter.values()

    def get_region(self, name: str) -> Region:
        try:
            region = self._regions_getter[name]
        except KeyError:
            raise KeyError(f"RGP with identifier {name} is not associated to system {self.ID}")
        else:
            return region

    def add_region(self, region: Region):
        try:
            region_in = self.get_region(name=region.name)
        except KeyError:
            self._regions_getter[region.name] = region
        else:
            if region != region_in:
                raise Exception(f"Another region with identifier {region.name} is already associated to system {self.ID}."
                                f"This is unexpected. Please report an issue on our GitHub")
