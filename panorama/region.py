"""
This module contains classes to represent regions, spots, conserved spots, and modules in a pangenome.
"""

# default libraries
from typing import Generator, List, Set, Tuple

# installed libraries
from ppanggolin.genome import Organism
from ppanggolin.region import (
    GeneContext as GeneCont,
)
from ppanggolin.region import (
    Module as Mod,
)
from ppanggolin.region import (
    Region as RGP,
)
from ppanggolin.region import (
    Spot as Hotspot,
)

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
        except KeyError as error:
            raise KeyError(f"Spot {key} is not in conserved spot with ID {self.ID}") from error

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
        self._unit_getter = {}
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
    def gene_families(self) -> Generator[GeneFamily, None, None]:
        """
        Get the set of gene families that define the module.

        Returns:
            GeneFamily: The set of gene families.
        """
        return super().families

    @property
    def units(self):
        """
        Generator of the systems associated with the module.

        Yields:
            Generator[SystemUnit]: The next system associated with the module.
        """
        yield from self._unit_getter.values()

    def get_unit(self, identifier: int):
        """
        Get a unit associated with the module.

        Args:
            identifier (int): The identifier of the unit.

        Returns:
            System: The unit with the given identifier.

        Raises:
            KeyError: If the unit is not associated with the module.
        """
        try:
            unit = self._unit_getter[identifier]
        except KeyError as error:
            raise KeyError(f"System {identifier} is not associated to module {self.ID}") from error
        else:
            return unit

    def add_unit(self, unit):
        """
        Add a system to the module.

        Args:
            unit (System): The system to add to the module.

        Raises:
            Exception: If a system with the same ID but different name or
            gene families is already associated with the module.
        """
        try:
            unit_in = self.get_unit(unit.ID)
        except KeyError:
            self._unit_getter[unit.ID] = unit
        else:
            if unit.name != unit_in.name:
                raise Exception(
                    "Two system with same ID but with different name are trying to be added to module."
                    "This error is unexpected. Please report on our GitHub"
                )
            else:
                if unit.families != unit_in.families:
                    raise Exception(
                        "Two system with same ID and name but with different gene families are trying to be"
                        " added to module. This error is unexpected. Please report on our GitHub"
                    )

    @property
    def systems(self):
        """
        Generator of the systems associated with the module.

        Yields:
            System: The next system associated with the module.
        """
        systems = set()
        for unit in self.units:
            systems.add(unit.system)
        yield from systems


class GeneContext(GeneCont):
    """
    A class used to represent a gene context

    """

    def __init__(
        self,
        pangenome,
        gc_id: int,
        families: Set[GeneFamily],
        families_of_interest: Set[GeneFamily],
    ):
        """
        :param gc_id : identifier of the Gene context
        :param families: Gene families included in the GeneContext
        """

        super().__init__(
            gc_id=f"{pangenome.name}_{gc_id}",
            families=families,
            families_of_interest=families_of_interest,
        )

        self.pangenome = pangenome.name

    def summarize(self) -> dict:
        """
        Summarize gene context information in a dict

        :return: dict with gene context info.
        """

        return {
            "GeneContext ID": self.ID,
            "pangenome": self.pangenome,
            "Gene Family count": len(self.families),
            "Partitions": ";".join({f.named_partition for f in self.families}),
        }

    def __gt__(self, other):
        return self.ID > other.ID

    def __lt__(self, other):
        return self.ID < other.ID
