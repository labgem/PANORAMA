#!/usr/bin/env python3
# coding: utf8

# default libraries
from typing import Generator, List, Tuple

# installed libraries
from ppanggolin.genome import Organism
from ppanggolin.region import Spot as Hotspot, Module as Mod

# local libraries
from panorama.geneFamily import GeneFamily
# from panorama.systems.system import System


class Spot(Hotspot):
    """"""

    def __init__(self, spot_id: int):
        """Constructor method
        """
        super().__init__(spot_id)
        self.pangenome = None
        self.conserved_id = None
        self._organisms_getter = None

    @property
    def conserved(self) -> True:
        """Return True if the spot is conserved between pangenomes, False otherwise."""
        return True if self.conserved_id is not None else False

    def _get_organisms(self):
        self._organisms_getter = {}
        for region in self.regions:
            self._organisms_getter[region.organism.name] = region.organism

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        if self._organisms_getter is None:
            self._get_organisms()
        yield from self._organisms_getter.values()


class ConservedSpots:
    def __init__(self, identifier: int, *spots: Spot):
        """Constructor method
        """
        self.ID = identifier
        self._spots_getter = {}
        for spot in spots:
            self.add(spot)

    def __setitem__(self, key: Tuple[str, int], value: Spot):
        try:
            self._spots_getter[key]
        except KeyError:
            self._spots_getter[key] = value
        else:
            raise KeyError(f"Spot {key} already exist in conserved spots with ID {self.ID}")

    def __getitem__(self, key: Tuple[str, int]) -> Spot:
        try:
            return self._spots_getter[key]
        except KeyError:
            raise KeyError(f"Spot {key} is not in conserved spot with ID {self.ID}")

    def add(self, spot: Spot) -> None:
        """
        Add a spot in the conserved set of spots

        Args:
            spot: spot to add in the object
        """
        assert isinstance(spot, Spot), f"Spot object is expected, given type is {type(spot)}"
        self[(spot.pangenome.name, spot.ID)] = spot
        spot.conserved_id = self.ID

    def get(self, spot_id: int, pangenome_name: str) -> Spot:
        """
        Get a spot in the conserved set of spots
        Args:
            spot_id: spot identifier
            pangenome_name: name of the pangenome from which the spot belongs

        Returns:
            The spot with the given id and pangenome
        """
        assert isinstance(spot_id, int), f"Spot id should be an integer, given type is {type(spot_id)}"
        return self[(pangenome_name, spot_id)]

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """Generator of the spots in the conserved object"""
        for spot in self._spots_getter.values():
            yield spot

    def pangenomes(self) -> List[str]:
        """
        Get the list of pangenomes where the conserved spot belongs

        Returns:
            List of pangenome
        """
        return [k[0] for k in self._spots_getter.keys()]


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
    def organisms(self):
        organisms = set()
        for family in self.families:
            organisms |= family.organisms
        return organisms

    @property
    def gene_families(self) -> GeneFamily:
        return super().families

    @property
    def systems(self):
        for system in self._systemsGetter.values():
            yield system

    def get_system(self, identifier: int):
        try:
            system = self._systemsGetter[identifier]
        except KeyError:
            raise KeyError(f"System {identifier} is not associated to module {self.ID}")
        else:
            return system

    def add_system(self, system):
        try:
            self.get_system(system.ID)
        except KeyError:
            self._systemsGetter[system.ID] = system
        else:
            if system.name != self._systemsGetter[system.ID].name:
                raise Exception("Two system with same ID but with different name are trying to be added to module."
                                "This error is unexpected. Please report on our GitHub")
            else:
                if system.families != self._systemsGetter[system.ID].families:
                    raise Exception("Two system with same ID and name but with different gene families are trying to be"
                                    " added to module. This error is unexpected. Please report on our GitHub")
