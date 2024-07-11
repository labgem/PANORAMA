"""
This module contains classes to represent regions, spots, conserved spots, and modules in a pangenome.
"""

# default libraries
from typing import Generator, List, Tuple

# installed libraries
from ppanggolin.genome import Organism
from ppanggolin.region import Region as RGP, Spot as Hotspot, Module as Mod

# local libraries
from panorama.geneFamily import GeneFamily


class Region(RGP):
    """
    Represents a region in a pangenome.

    Args:
        name (str): The name of the region.
    """

    def __init__(self, name: str):
        """
        Constructor method.

        Args:
            name (str): The name of the region.
        """
        super().__init__(name)


class Spot(Hotspot):
    """
    Represents a spot in a pangenome.
    """

    def __init__(self, spot_id: int):
        """
        Constructor method.

        Args:
            spot_id (int): The identifier of the spot.
        """
        super().__init__(spot_id)
        self.pangenome = None
        self.conserved_id = None
        self._organisms_getter = None

    @property
    def conserved(self) -> bool:
        """
        Check if the spot is conserved between pangenomes.

        Returns:
            bool: True if the spot is conserved, False otherwise.
        """
        return True if self.conserved_id is not None else False

    def _get_organisms(self):
        self._organisms_getter = {}
        for region in self.regions:
            self._organisms_getter[region.organism.name] = region.organism

    @property
    def organisms(self) -> Generator[Organism, None, None]:
        """
        Generator of the organisms that contain the spot.

        Yields:
            Organism: The next organism that contains the spot.
        """
        if self._organisms_getter is None:
            self._get_organisms()
        yield from self._organisms_getter.values()

    @property
    def number_of_organisms(self):
        """
        Get the number of organisms that contain the spot.

        Returns:
            int: The number of organisms.
        """
        if self._organisms_getter is None:
            self._get_organisms()
        return len(self._organisms_getter)


class ConservedSpots:
    """
    Represents a set of conserved spots across multiple pangenomes.
    """

    def __init__(self, identifier: int, *spots: Spot):
        """
        Constructor method.

        Args:
            identifier (int): The identifier of the conserved spots set.
            *spots (Spot): The spots to add to the conserved set.
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
            raise KeyError(f"Spot {key} already exists in conserved spots with ID {self.ID}")

    def __getitem__(self, key: Tuple[str, int]) -> Spot:
        try:
            return self._spots_getter[key]
        except KeyError:
            raise KeyError(f"Spot {key} is not in conserved spot with ID {self.ID}")

    def add(self, spot: Spot) -> None:
        """
        Add a spot to the conserved set of spots.

        Args:
            spot (Spot): The spot to add to the object.

        Raises:
            AssertionError: If the spot is not an instance of the Spot class.
            KeyError: If the spot already exists in the conserved spots with the same ID.
        """
        assert isinstance(spot, Spot), f"Spot object is expected, given type is {type(spot)}"
        self[(spot.pangenome.name, spot.ID)] = spot
        spot.conserved_id = self.ID

    def get(self, spot_id: int, pangenome_name: str) -> Spot:
        """
        Get a spot from the conserved set of spots.

        Args:
            spot_id (int): The identifier of the spot.
            pangenome_name (str): The name of the pangenome from which the spot belongs.

        Returns:
            Spot: The spot with the given id and pangenome.

        Raises:
            AssertionError: If the spot id is not an integer.
            KeyError: If the spot is not in the conserved set of spots.
        """
        assert isinstance(spot_id, int), f"Spot id should be an integer, given type is {type(spot_id)}"
        return self[(pangenome_name, spot_id)]

    @property
    def spots(self) -> Generator[Spot, None, None]:
        """
        Generator of the spots in the conserved object.

        Yields:
            Spot: The next spot in the conserved object.
        """
        for spot in self._spots_getter.values():
            yield spot

    def pangenomes(self) -> List[str]:
        """
        Get the list of pangenomes where the conserved spot belongs.

        Returns:
            List[str]: The list of pangenome names.
        """
        return [k[0] for k in self._spots_getter.keys()]


class Module(Mod):
    """
    Represents a module in a pangenome.

    Args:
        module_id (int): The identifier of the module.
        families (set, optional): The set of families that define the module.
    """

    def __init__(self, module_id: int, families: set = None):
        """
        Constructor method.

        Args:
            module_id (int): The identifier of the module.
            families (set, optional): The set of families that define the module. Defaults to None.
        """
        super().__init__(module_id=module_id, families=families)
        self._systemsGetter = {}
        self._sys2fam = {}

    @property
    def organisms(self):
        """
        Get the set of organisms that contain the module.

        Returns:
            set: The set of organisms.
        """
        organisms = set()
        for family in self.families:
            organisms |= set(family.organisms)
        return organisms

    @property
    def number_of_organisms(self):
        """
        Get the number of organisms that contain the module.

        Returns:
            int: The number of organisms.
        """
        return len(set(self.organisms))

    @property
    def gene_families(self) -> GeneFamily:
        """
        Get the set of gene families that define the module.

        Returns:
            GeneFamily: The set of gene families.
        """
        return super().families

    @property
    def systems(self):
        """
        Generator of the systems associated with the module.

        Yields:
            System: The next system associated with the module.
        """
        for system in self._systemsGetter.values():
            yield system

    def get_system(self, identifier: int):
        """
        Get a system associated with the module.

        Args:
            identifier (int): The identifier of the system.

        Returns:
            System: The system with the given identifier.

        Raises:
            KeyError: If the system is not associated with the module.
        """
        try:
            system = self._systemsGetter[identifier]
        except KeyError:
            raise KeyError(f"System {identifier} is not associated to module {self.ID}")
        else:
            return system

    def add_system(self, system):
        """
        Add a system to the module.

        Args:
            system (System): The system to add to the module.

        Raises:
            Exception: If a system with the same ID but different name or gene families is already associated with the module.
        """
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
