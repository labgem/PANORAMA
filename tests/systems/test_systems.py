#!/usr/bin/env python3
# coding: utf8

"""
Unit tests for the system module containing SystemUnit, System, and ClusterSystems classes.

This test suite covers all major functionality including:
- SystemUnit creation, manipulation, and set operations
- System creation, unit management, and comparisons
- ClusterSystems functionality for grouping homologous systems
- Integration with gene families, organisms, modules, spots, and regions
"""
from random import randint, choice
from typing import Set

import pytest

from ppanggolin.metadata import Metadata
from ppanggolin.genome import Organism, Gene

# Import the classes to test
from panorama.geneFamily import GeneFamily
from panorama.region import Region, Spot, Module
from panorama.systems.system import SystemUnit, System, ClusterSystems
from panorama.systems.models import Family, FuncUnit, Model


class TestFixture:
    """Base class for fixture and Mock definitions."""

    @pytest.fixture
    def model(self):
        """
        Fixture for creating a test instance of the Model class.

        Returns:
            Model: A pre-configured instance of the Model class.
        """
        model = Model(
            name="test_model",
            min_mandatory=1,
            min_total=1,
            canonical=["canonical_1", "canonical_2"],
        )
        return model

    @pytest.fixture
    def functional_unit(self, model):
        """Create a mock FuncUnit for testing."""
        fu = FuncUnit(name="test_unit", presence="mandatory", min_total=2)
        fu.model = model
        return fu


class TestSystemUnit(TestFixture):
    """Test class for SystemUnit functionality."""

    _organims_counter = 0
    _spots_counter = 0
    _modules_counter = 0

    def create_gene_family(self, name, identifier):
        """
        Creates a gene family with associated organisms, a module, and spots.

        Args:
            name: The name of the gene family to be created.
            identifier: The unique identifier for the gene family.

        Returns:
            MockGeneFamily: An object representing the created gene family.
        """
        gf = GeneFamily(name=name, family_id=identifier)
        for _ in range(2):
            self._organims_counter += 1
            gf._genePerOrg[Organism(name=f"organism_{self._organims_counter}")] = set()
        self._modules_counter += 1
        gf.module = Module(module_id=self._modules_counter)
        for _ in range(2):
            self._spots_counter += 1
            gf.add_spot(Spot(spot_id=self._spots_counter))
        return gf

    @pytest.fixture
    def gene_families(self):
        """Create a mock GeneFamily for testing."""
        return {self.create_gene_family("gf1", 1), self.create_gene_family("gf2", 2)}

    def test_init_basic(self, functional_unit):
        """Test basic SystemUnit initialization."""
        unit = SystemUnit(functional_unit, "test_source")

        assert unit.functional_unit == functional_unit
        assert unit.source == "test_source"
        assert unit._families_getter == {}
        assert unit._families2metainfo == {}
        assert unit._regions_getter == {}
        assert unit._spots_getter == {}
        assert unit._modules_getter == {}
        assert unit._model_families is None
        assert unit._system is None
        assert unit.ID == 1  # Auto-generated ID
        assert len(unit) == 0  # No gene families initially

    @pytest.fixture
    def system_unit(self, functional_unit):
        """Create a basic SystemUnit for testing."""
        return SystemUnit(functional_unit, "test_source")

    def test_repr(self, system_unit):
        """Test string representation of SystemUnit."""
        assert (
            repr(system_unit)
            == f"System unit {system_unit.ID}, name: {system_unit.name}, model: {system_unit.model.name}"
        )

    def test_system_property(self, system_unit):
        """Test system property getter and setter."""
        assert system_unit.system is None

        sys = System(system_unit.model, system_unit.source)
        system_unit.system = sys
        assert system_unit.system == sys

    @pytest.mark.parametrize("invalid", ["invalid", 1, 1.0, {1, 2}, [1, 2], (1, 2)])
    def test_system_setter_invalid_type(self, system_unit, invalid):
        """Test that setting invalid system types raises TypeError."""
        with pytest.raises(
            TypeError, match=f"System must be an instance of {System.__name__}."
        ):
            system_unit.system = invalid

    def test_functional_unit_property(self, system_unit, functional_unit):
        """Test functional unit property."""
        assert system_unit.functional_unit == functional_unit

    @pytest.mark.parametrize("invalid", ["invalid", 1, 1.0, {1, 2}, [1, 2], (1, 2)])
    def test_functional_unit_setter_invalid_type(self, system_unit, invalid):
        """Test that setting invalid functional unit types raises TypeError."""
        with pytest.raises(TypeError):
            system_unit.functional_unit = invalid

    def test_model_property(self, system_unit, functional_unit):
        """Test model property access through functional unit."""
        assert system_unit.model == functional_unit.model

    def test_name_property(self, functional_unit, system_unit):
        """Test name property access through functional unit."""
        assert system_unit.name == functional_unit.name

    def test_setitem_getitem(self, system_unit, gene_families):
        """Test setting and getting gene families by name."""
        gf1, gf2 = gene_families
        system_unit[gf1.name] = gf1
        system_unit[gf2.name] = gf2
        assert system_unit._families_getter[gf1.name] == gf1
        assert system_unit._families_getter[gf2.name] == gf2

    @pytest.mark.parametrize("invalid", ["invalid", 1, 1.0, {1, 2}, [1, 2], (1, 2)])
    def test_setitem_invalid_type(self, system_unit, invalid):
        """Test that setting non-GeneFamily raises TypeError."""
        with pytest.raises(
            TypeError,
            match=f"A GeneFamily object is expected. You provided a {type(invalid)}.",
        ):
            system_unit["invalid"] = invalid

    def test_setitem_duplicate_name_different_family(self, system_unit, gene_families):
        """Test that setting different family with existing name raises KeyError."""
        gf1, _ = gene_families
        system_unit[gf1.name] = gf1
        gf3 = GeneFamily(name=gf1.name, family_id=3)
        gf3.add(Gene(gene_id="test_gene"))

        with pytest.raises(
            KeyError,
            match="A different gene family with the same name already exist in the functional unit",
        ):
            system_unit[gf3.name] = gf3

    def test_getitem_nonexistent(self, system_unit):
        """Test that getting non-existent family raises KeyError."""
        with pytest.raises(KeyError):
            _ = system_unit["nonexistent"]

    def test_add_family(self, system_unit, gene_families):
        """Test adding gene family to system unit."""
        gf1, gf2 = gene_families
        gf1.add_metadata(Metadata("test_source", annot="test_func1"), 1)
        gf2.add_metadata(Metadata("other_source", annot="test_func2"), 1)
        gf2.add_metadata(Metadata("test_source", annot="test_func3"), 2)
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)

        assert system_unit[gf1.name] == gf1
        assert system_unit[gf2.name] == gf2

    @pytest.mark.parametrize("invalid", ["invalid", 1, 1.0, {1, 2}, [1, 2], (1, 2)])
    def test_add_family_invalid_type(self, system_unit, invalid):
        """Test that adding non-GeneFamily raises AssertionError."""
        with pytest.raises(AssertionError, match="GeneFamily object is expected"):
            system_unit.add_family(invalid, "test_source", 1)

    def test_init_with_gene_families(self, functional_unit, gene_families):
        """Test SystemUnit initialization with gene families."""
        gf1, gf2 = gene_families
        families_to_metainfo = {gf1: ("test_source", 1), gf2: ("test_source", 2)}

        unit = SystemUnit(
            functional_unit, "test_source", gene_families, families_to_metainfo
        )

        assert len(unit) == 2
        assert unit[gf1.name] == gf1
        assert unit[gf2.name] == gf2
        assert unit._families2metainfo == {
            gf1: ("test_source", 1),
            gf2: ("test_source", 2),
        }

    def test_hash_and_equality(self, functional_unit, gene_families):
        """Test SystemUnit hashing and equality comparison."""
        unit1 = SystemUnit(functional_unit, "source1")
        unit2 = SystemUnit(functional_unit, "source2")
        gf1, gf2 = gene_families
        unit1.add_family(gf1, "test_source", 1)
        unit2.add_family(gf1, "test_source", 1)
        # Empty units should be equal and have same hash
        assert hash(unit1) == hash(unit2)

        # Add family to one unit
        unit1.add_family(gf2, "test_source", 2)

        # Now they should be different
        assert hash(unit1) != hash(unit2)

    def test_families_property(self, system_unit, gene_families):
        """Test families generator property."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        families = set(system_unit.families)
        assert len(families) == 2
        assert families == {gf1, gf2}

    def test_get_model_gene_families(self, system_unit, gene_families):
        """Test getting gene families."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)

        assert system_unit._get_model_families() == {gf1, gf2}

    def test_model_families_property(self, system_unit, gene_families):
        """Test model families property."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        model_families = set(system_unit.model_families)
        assert model_families == {gf1, gf2}
        assert system_unit._model_families == model_families

    def test_equality_comparison(self, functional_unit, gene_families):
        """Test equality comparison between SystemUnits."""
        gf1, gf2 = gene_families
        fam2meta_info = {gf1: ("test_source", 1), gf2: ("test_source", 2)}
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)

        assert unit1 == unit2

        unit1.add_family(GeneFamily(3, "new_family"), "test_source", 1)

        assert unit1 != unit2

    def test_is_superset(self, functional_unit, gene_families):
        """Test superset comparison between SystemUnits."""
        gf1, gf2 = gene_families
        fam2meta_info = {gf1: ("test_source", 1), gf2: ("test_source", 2)}
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", {gf1}, fam2meta_info)

        assert unit1.is_superset(unit2)
        assert not unit2.is_superset(unit1)

    def test_is_subset(self, functional_unit, gene_families):
        """Test subset comparison between SystemUnits."""
        gf1, gf2 = gene_families
        fam2meta_info = {gf1: ("test_source", 1), gf2: ("test_source", 2)}
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", {gf1}, fam2meta_info)

        assert not unit1.is_subset(unit2)
        assert unit2.is_subset(unit1)

    def test_intersection(self, functional_unit, gene_families):
        """Test intersection between SystemUnits."""
        gf1, gf2 = gene_families
        fam2meta_info = {gf1: ("test_source", 1), gf2: ("test_source", 2)}
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", {gf1}, fam2meta_info)

        intersection = unit1.intersection(unit2)
        assert len(intersection) == 1
        assert intersection == {gf1}

    def test_difference(self, functional_unit, gene_families):
        """Test difference between SystemUnits."""
        gf1, gf2 = gene_families
        fam2meta_info = {gf1: ("test_source", 1), gf2: ("test_source", 2)}
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", {gf1}, fam2meta_info)

        assert unit1.difference(unit2) == {gf2}
        assert unit2.difference(unit1) == set()

    def test_symmetric_difference(self, functional_unit, gene_families):
        """Test symmetric difference between SystemUnits."""
        gf1, gf2 = gene_families
        fam2meta_info = {gf1: ("test_source", 1), gf2: ("test_source", 2)}
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", {gf1}, fam2meta_info)

        assert (
            unit1.symmetric_difference(unit2)
            == unit2.symmetric_difference(unit1)
            == {gf2}
        )

    def test_merge(self, functional_unit, gene_families):
        """Test merging two SystemUnits."""
        gf1, gf2 = gene_families
        gf3 = GeneFamily(3, "new_family")
        fam2meta_info = {
            gf1: ("test_source", 1),
            gf2: ("test_source", 2),
            gf3: ("test_source", 1),
        }
        unit1 = SystemUnit(functional_unit, "test_source", gene_families, fam2meta_info)
        unit2 = SystemUnit(functional_unit, "test_source", {gf1, gf3}, fam2meta_info)
        initial_len = len(unit1)

        unit1.merge(unit2)

        assert len(unit1) == initial_len + 1
        assert unit1["new_family"] == gf3
        assert set(unit1.families) == {gf1, gf2, gf3}

    def test_organisms_property(self, system_unit, gene_families):
        """Test organisms property."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        organisms = list(system_unit.organisms)

        assert len(organisms) == gf1.number_of_organisms + gf2.number_of_organisms
        assert set(organisms) == set(gf1.organisms).union(set(gf2.organisms))

    def test_nb_organisms_property(self, system_unit, gene_families):
        """Test number of organisms property."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)

        assert (
            system_unit.nb_organisms
            == gf1.number_of_organisms + gf2.number_of_organisms
        )

    def test_model_organisms_property(self, system_unit, gene_families):
        """Test model organisms property."""
        gf1, gf2 = gene_families
        gf3 = GeneFamily(3, "new_family")
        gf3.add_metadata(Metadata("test_source", annot="test_func3"), 1)
        for org in set(gf1.organisms).union(set(gf2.organisms)):
            gf3._genePerOrg[org] = {}
        gf4 = GeneFamily(4, "context_family")
        organism = Organism("test_organism")
        gf4._genePerOrg[organism] = set()
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        system_unit.add_family(gf3, "test_source", 1)
        system_unit.add_family(gf4)

        assert set(system_unit.model_organisms) == set(gf1.organisms).union(
            set(gf2.organisms)
        ).union(set(gf3.organisms))
        assert organism not in system_unit.model_organisms

    def test_annotation_sources(self, system_unit, gene_families):
        """Test annotation sources retrieval."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "other_source", 1)

        sources = system_unit.annotation_sources()
        assert "test_source" in sources
        assert "other_source" in sources

    def test_get_metainfo(self, system_unit, gene_families):
        """Test getting metainfo from system unit."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)

        assert system_unit.get_metainfo(gf1) == ("test_source", 1)
        assert system_unit.get_metainfo(gf2) == ("test_source", 2)

    def test_add_module(self, system_unit):
        """Test adding module to system unit."""

        module = Module(module_id=1)
        system_unit.add_module(module)

        assert system_unit._modules_getter[1] == module

    def test_get_module(self, system_unit):
        """Test getting module from system unit."""
        module = Module(module_id=1)
        system_unit.add_module(module)
        assert system_unit.get_module(1) == module

    def test_get_module_nonexistent(self, system_unit):
        """Test getting non-existent module raises KeyError."""
        with pytest.raises(KeyError):
            system_unit.get_module(999)

    def test_modules_property(self, system_unit, gene_families):
        """Test associating module from gene families to system unit."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)

        assert set(system_unit.modules) == {gf1.module, gf2.module}

    def test_add_spot(self, system_unit):
        """Test adding spot to system unit."""
        spot = Spot(spot_id=1)
        system_unit.add_spot(spot)
        assert system_unit._spots_getter[1] == spot

    def test_get_spot(self, system_unit):
        """Test getting spot from system unit."""
        spot = Spot(spot_id=1)
        system_unit.add_spot(spot)
        assert system_unit.get_spot(1) == spot

    def test_get_spot_nonexistent(self, system_unit):
        """Test getting non-existent spot raises KeyError."""
        with pytest.raises(KeyError):
            system_unit.get_spot(999)

    def test_spots_property_with_gf(self, system_unit, gene_families):
        """Test getting spots associated with gene families from system unit."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        assert set(system_unit.spots) == set(gf1.spots).union(set(gf2.spots))

    def test_add_region(self, system_unit):
        """Test adding region to system unit."""
        rgp = Region("test_region")
        system_unit.add_region(rgp)

        assert system_unit._regions_getter["test_region"] == rgp

    def test_get_region(self, system_unit):
        """Test getting region from system unit."""
        rgp = Region("test_region")
        system_unit.add_region(rgp)
        assert system_unit.get_region("test_region") == rgp

    def test_get_region_nonexistent(self, system_unit):
        """Test getting non-existent region raises KeyError."""
        with pytest.raises(KeyError):
            system_unit.get_region("nonexistent")

    def test_regions_property(self, system_unit):
        """Test getting regions from system unit."""
        rgp1, rgp2 = Region("rgp1"), Region("rgp2")
        system_unit.add_region(rgp1)
        system_unit.add_region(rgp2)
        assert set(system_unit.regions) == {rgp1, rgp2}

    def test_spots_property_with_region(self, system_unit):
        """Test getting spots associated with regions from system unit."""
        rgp1, rgp2, rgp3, rgp4 = (
            Region("rgp1"),
            Region("rgp2"),
            Region("rgp3"),
            Region("rgp4"),
        )
        spot1, spot2 = Spot(spot_id=3), Spot(spot_id=4)
        rgp1.spot = spot1
        rgp2.spot = spot2
        rgp3.spot = spot2
        rgp4.spot = spot2
        system_unit.add_region(rgp1)
        system_unit.add_region(rgp2)
        system_unit.add_region(rgp3)
        system_unit.add_region(rgp4)
        assert set(system_unit.spots) == {spot1, spot2}


class TestSystem(TestFixture):
    """Test class for System functionality."""

    @pytest.fixture
    def unit(self, functional_unit, model) -> SystemUnit:
        """
        Fixture for creating a SystemUnit instance with test data.

        Args:
            functional_unit: The functional unit to be used in the SystemUnit.
            model: The model associated with the test setup.

        Returns:
            SystemUnit: A SystemUnit instance populated with test data.
        """
        gene_families = {GeneFamily(j, f"test_gf_{j}") for j in range(randint(2, 5))}
        fam2metainfo = {gf: ("test_source", randint(1, 5)) for gf in gene_families}
        unit = SystemUnit(
            functional_unit=functional_unit,
            source="test_source",
            gene_families=gene_families,
            families_to_metainfo=fam2metainfo,
        )
        return unit

    @pytest.fixture
    def units(self, model) -> Set[SystemUnit]:
        """
        Creates a fixture that generates a set of SystemUnit objects for a given model.

        Args:
            model: The model with which the functional units and system units are
                associated.

        Returns:
            Set[SystemUnit]: A set containing the generated SystemUnit objects.
        """
        units = set()
        counter_fam = 0
        spot_id = set()
        for i in range(1, 5):
            fu = FuncUnit(name=f"test_unit_{i}", presence="mandatory", min_total=1)
            fu.model = model
            gene_families = {
                GeneFamily(counter_fam + j, f"test_family_{counter_fam+j}")
                for j in range(randint(1, 5))
            }
            for j, gf in enumerate(gene_families):
                gf._genePerOrg = {
                    Organism(name=f"test_organism_{i}"): set()
                    for i in range(randint(4, 10))
                }
                if randint(0, 1):  # Add a module or not randomly
                    gf.module = Module(counter_fam + j)
                if randint(0, 1):
                    for j in range(randint(1, 5)):
                        new_id = randint(1, 100)
                        while new_id in spot_id:
                            new_id = randint(1, 100)
                        spot_id.add(new_id)
                        gf.add_spot(Spot(spot_id=new_id))

            counter_fam += len(gene_families)
            fam2metainfo = {gf: ("test_source", randint(1, 5)) for gf in gene_families}
            unit = SystemUnit(
                functional_unit=fu,
                source="test_source",
                gene_families=gene_families,
                families_to_metainfo=fam2metainfo,
            )
            units.add(unit)
        return units

    @pytest.fixture
    def system(self, model):
        """Create a basic System for testing."""
        return System(model, "test_source")

    def test_init_basic(self, model):
        """Test basic System initialization."""
        system = System(model, "test_source")

        assert system.model == model
        assert system.source == "test_source"
        assert system.ID == "1"
        assert system._unit_getter == {}
        assert system.canonical == set()
        assert system._fam2unit is None
        assert system.pangenome is None
        assert system.cluster_id is None

    def test_init_with_custom_id(self, model):
        """Test System initialization with custom ID."""
        system = System(model, "test_source", system_id="custom_id")
        assert system.ID == "custom_id"

    def test_init_with_units(self, model, units):
        """Test System initialization with units."""
        system = System(model, "test_source", units=units)

        assert system._unit_getter == {unit.name: unit for unit in units}

    def test_hash(self, units):
        """Test System hashing."""
        units_list = list(units)
        system_1 = System(units_list[0].model, "test_source", units=set(units_list[:2]))
        hash1 = hash(system_1)
        system_2 = System(units_list[2].model, "test_source", units=set(units_list[2:]))
        hash2 = hash(system_2)

        assert hash1 != hash2

    def test_repr(self, system):
        """Test System string representation."""
        assert repr(system) == f"System ID: {system.ID}, Name: {system.name}"

    def test_len(self, system, units):
        """Test System length (number of units)."""
        assert len(system) == 0
        for unit in units:
            system._unit_getter[unit.name] = unit
        assert len(system) == 4

    def test_setitem_getitem(self, system, unit):
        """Test setting and getting units by name."""
        system[unit.name] = unit
        assert system[unit.name] == unit

    def test_setitem_invalid_type(self, system):
        """Test that setting non-SystemUnit raises TypeError."""
        with pytest.raises(TypeError):
            system["invalid"] = "not_a_system_unit"

    def test_getitem_nonexistent(self, system):
        """Test that getting non-existent unit raises KeyError."""
        with pytest.raises(KeyError):
            _ = system["nonexistent"]

    def test_delitem(self, system, unit):
        """Test deleting unit from system."""
        system[unit.name] = unit
        del system[unit.name]

        with pytest.raises(KeyError):
            _ = system[unit.name]

    def test_delitem_nonexistent(self, system):
        """Test that deleting non-existent unit raises KeyError."""
        with pytest.raises(KeyError):
            del system["nonexistent"]

    def test_equality(self, model, unit, units):
        """Test System equality comparison."""
        system1 = System(model, "source1", units=units)
        system2 = System(model, "source1", units=units)

        # Empty systems should be equal
        assert system1 == system2

        # Add new unit to one system
        system1[unit.name] = unit
        assert system1 != system2

    def test_equality_invalid_type(self, system):
        """Test that comparing with non-System raises TypeError."""
        with pytest.raises(TypeError):
            _ = system == "not_a_system"

    def test_name_property(self, system):
        """Test name property access through model."""
        assert system.name == system.model.name == "test_model"

    def test_units_property(self, system, units):
        """Test units generator property."""
        system._unit_getter = {unit.name: unit for unit in units}

        units_set = set(system.units)
        assert len(units_set) == 4
        assert units_set == units

    def test_add_unit(self, system, unit):
        """Test adding unit to system."""
        unit._system = None  # Reset system reference
        system.add_unit(unit)

        assert system[unit.name] == unit
        assert unit.system == system

    def test_add_unit_invalid_type(self, system):
        """Test that adding non-SystemUnit raises AssertionError."""
        with pytest.raises(TypeError):
            system.add_unit("not_a_system_unit")

    def test_add_unit_subset_do_nothing(self, system, unit):
        """Test that adding subset unit does not change anything."""
        gene_families = set(list(unit.families)[:-1])
        subset_unit = SystemUnit(
            unit.functional_unit,
            "test_source",
            gene_families=gene_families,
            families_to_metainfo=unit._families2metainfo,
        )

        assert subset_unit.is_subset(unit) is True

        system.add_unit(unit)

        assert system[unit.name] == unit

        system.add_unit(subset_unit)

        assert system[unit.name] == system[subset_unit.name] == unit
        assert system[unit.name] == system[subset_unit.name] != subset_unit

    def test_add_unit_superset_replacement(self, system, unit, model):
        """Test that adding superset unit replaces existing subset."""
        # Create two units where one is superset of the other
        new_gene_family = GeneFamily(100, "new_family")
        fam2metainfo = unit._families2metainfo
        fam2metainfo[new_gene_family] = ("test_source", 1)
        superset_unit = SystemUnit(
            unit.functional_unit,
            "test_source",
            gene_families=set(unit.families).union({new_gene_family}),
            families_to_metainfo=fam2metainfo,
        )

        assert superset_unit.is_superset(unit) is True

        system.add_unit(unit)

        assert system[unit.name] == unit

        system.add_unit(superset_unit)  # Should replace unit1

        assert system[unit.name] == system[superset_unit.name] == superset_unit

    def test_get_unit(self, system, unit):
        """Test getting unit by name."""
        system.add_unit(unit)
        assert system.get_unit(unit.name) == unit

    def test_families_property(self, system, units):
        """Test families property aggregating from all units."""
        for unit in units:
            system.add_unit(unit)
        families = set(system.families)

        assert families == {fam for unit in units for fam in unit.families}

    def test_number_of_families(self, system, units):
        """Test counting total families across units."""
        for unit in units:
            system.add_unit(unit)
        families = set(system.families)

        assert system.number_of_families == len(
            {fam for unit in units for fam in unit.families}
        )

    def test_model_families_property(self, system, units):
        """Test families property aggregating from all units."""
        for unit in units:
            system.add_unit(unit)

        model_families = set(system.model_families)

        assert model_families == {fam for unit in units for fam in unit.model_families}

    def test_number_of_model_gene_families(self, system, units):
        """Test counting model gene families across units."""
        for unit in units:
            system.add_unit(unit)

        model_families = set(system.model_families)

        assert system.number_of_model_gene_families == len(model_families)

    def test_is_superset(self, model, units):
        """Test superset comparison between Systems."""

        def create_superset_units() -> tuple[Set[int], Set[SystemUnit]]:
            """
            Generates a superset of units by augmenting the input units with new gene families.

            Returns:
                tuple[Set[int], Set[SystemUnit]]: A tuple containing:
                    - A set of integer choices where each value represents the random decision for
                      whether to modify the corresponding unit.
                    - A set of `SystemUnit` objects representing the resulting superset of units,
                      including the augmented units where applicable.
            """
            choices = set()
            new_units = set()
            for unit in units:
                choice = randint(0, 1)
                choices.add(choice)
                if choice:
                    fam_id = randint(50, 100)
                    new_gene_family = GeneFamily(fam_id, f"new_family_{fam_id}")
                    fam2metainfo = unit._families2metainfo
                    fam2metainfo[new_gene_family] = ("test_source", 1)
                    superset_unit = SystemUnit(
                        unit.functional_unit,
                        "test_source",
                        gene_families=set(unit.families).union({new_gene_family}),
                        families_to_metainfo=fam2metainfo,
                    )
                    new_units.add(superset_unit)
                else:
                    new_units.add(unit)
            return choices, new_units

        do_create_superset_units = set()
        new_sys_units = set()
        while do_create_superset_units != {0, 1}:
            do_create_superset_units, new_sys_units = create_superset_units()

        system1 = System(model, "test_source")
        system2 = System(model, "test_source")

        for sys_unit in units:
            system1.add_unit(sys_unit)
        for sys_unit in new_sys_units:
            system2.add_unit(sys_unit)

        assert system2.is_superset(system1)

    def test_is_subset(self, model, units):
        """Test subset comparison between Systems."""

        def create_subset_units() -> tuple[Set[int], Set[SystemUnit]]:
            """
            Generates a superset of units by augmenting the input units with new gene families.

            This function randomly determines, for each unit, whether to create a modified version
            of the unit by appending a new gene family to it. The original or newly created unit is
            then added to the resulting superset. The function returns the random choices and the
            augmented set of units.

            Returns:
                tuple[Set[int], Set[SystemUnit]]: A tuple containing:
                    - A set of integer choices where each value represents the random decision for
                      whether to modify the corresponding unit.
                    - A set of `SystemUnit` objects representing the resulting superset of units,
                      including the augmented units where applicable.
            """
            choices = set()
            new_units = set()
            for unit in units:
                choice = randint(0, 1)
                choices.add(choice)
                if choice:
                    unit_families = list(unit.families)
                    unit_families.pop(randint(0, len(unit_families) - 1))
                    new_gene_families = set(unit_families)
                    fam2metainfo = unit._families2metainfo
                    superset_unit = SystemUnit(
                        unit.functional_unit,
                        "test_source",
                        gene_families=new_gene_families,
                        families_to_metainfo=fam2metainfo,
                    )
                    new_units.add(superset_unit)
                else:
                    new_units.add(unit)
            return choices, new_units

        do_create_subset_units = set()
        new_sys_units = set()
        while do_create_subset_units != {0, 1}:
            do_create_subset_units, new_sys_units = create_subset_units()

        system1 = System(model, "test_source")
        system2 = System(model, "test_source")

        for sys_unit in units:
            system1.add_unit(sys_unit)
        for sys_unit in new_sys_units:
            system2.add_unit(sys_unit)

        assert system2.is_subset(system1)

    def test_intersection(self, model, units):
        """Test intersection between Systems."""
        system1 = System(model, "source1")
        system2 = System(model, "source2")

        units_list = list(units)
        for unit in units_list[:-1]:
            system1.add_unit(unit)

        for unit in units_list[1:]:
            system2.add_unit(unit)

        intersection = system1.intersection(system2)
        assert units_list[1:-1] == list(intersection)

    def test_is_superset_subset_with_missing_units(self, model, units):
        """Test superset comparison between Systems, when one or more units are missing."""
        system1 = System(model, "test_source")
        system2 = System(model, "test_source")

        for sys_unit in units:
            system1.add_unit(sys_unit)
        for sys_unit in list(units)[: -randint(2, len(units) - 1)]:
            system2.add_unit(sys_unit)

        assert system2.is_subset(system1)
        assert system1.is_superset(system2)

    def test_merge(self, model, units):
        """Test merging two Systems."""
        system1 = System(model, "source1")
        system2 = System(model, "source2")

        units_list = list(units)
        for unit in units_list[:-1]:
            system1.add_unit(unit)

        for unit in units_list[1:]:
            system2.add_unit(unit)

        system1.merge(system2)

        assert len(system1) == len(set(units))
        assert system1.get_unit(units_list[-1].name) == units_list[-1]

        with pytest.raises(KeyError):
            system2.get_unit(units_list[0].name)

    def test_organisms_property(self, system, units):
        """Test organisms property aggregating from all units."""
        for unit in units:
            system.add_unit(unit)
        organisms = set(system.organisms)
        expected_organisms = {
            org for unit in units for fam in unit.families for org in fam.organisms
        }

        assert organisms == expected_organisms

    def test_model_organisms_property(self, system, units):
        """Test model organisms property aggregating from all units."""
        for unit in units:
            system.add_unit(unit)
        organisms = set(system.model_organisms)
        expected_organisms = {
            org
            for unit in units
            for fam in unit.model_families
            for org in fam.organisms
        }

        assert organisms == expected_organisms

    def test_canonical_models(self, system, model):
        """Test canonical models retrieval."""
        canonical_models = system.canonical_models()
        assert canonical_models == model.canonical == ["canonical_1", "canonical_2"]

    def test_add_canonical(self, system, model):
        """Test adding canonical system."""
        canonical_model = Model(name=model.canonical[0], min_mandatory=1, min_total=1)
        canonical_system = System(canonical_model, "test_source")

        system.add_canonical(canonical_system)
        assert canonical_system in system.canonical

    def test_add_canonical_already_in_and_is_subset(self, system, model, units):
        """Test adding canonical system that is already in the system."""
        canonical_model = Model(name=model.canonical[0], min_mandatory=1, min_total=1)
        canonical_system = System(canonical_model, "test_source", system_id="Canon_1")
        for unit in units:
            canonical_system.add_unit(unit)
        system.add_canonical(canonical_system)

        assert canonical_system in system.canonical
        sub_canonical_system = System(
            canonical_model, "test_source", system_id="Canon_2"
        )
        for unit in list(units)[:-1]:
            sub_canonical_system.add_unit(unit)

        system.add_canonical(sub_canonical_system)

        assert sub_canonical_system not in system.canonical

    def test_add_canonical_already_in_and_is_superset(self, system, model, units, unit):
        """Test adding canonical system that is already in the system."""
        canonical_model = Model(name=model.canonical[0], min_mandatory=1, min_total=1)
        canonical_system = System(canonical_model, "test_source", system_id="Canon_1")
        super_canonical_system = System(
            canonical_model, "test_source", system_id="Canon_2"
        )
        for sys_unit in units:
            canonical_system.add_unit(sys_unit)
            super_canonical_system.add_unit(sys_unit)
        system.add_canonical(canonical_system)

        assert canonical_system in system.canonical

        super_canonical_system.add_unit(unit)
        system.add_canonical(super_canonical_system)

        assert super_canonical_system in system.canonical
        assert canonical_system not in system.canonical

    def test_get_metainfo(self, system, unit):
        """Test getting metadata info for gene family."""

        system.add_unit(unit)
        fam2metainfo = unit._families2metainfo
        for gf in system.families:
            metainfo = system.get_metainfo(gf)
            assert metainfo == fam2metainfo[gf]

    def test_annotation_sources(self, system, model, units):
        """Test getting annotation sources from system."""
        for unit in units:
            system.add_unit(unit)

        new_fu = FuncUnit(name=f"new_test_unit", presence="mandatory", min_total=1)
        new_fu.model = model
        gene_families = {
            GeneFamily(100 + j, f"test_family_{100+j}") for j in range(randint(1, 5))
        }
        fam2metainfo = {gf: ("new_test_source", randint(1, 5)) for gf in gene_families}
        new_unit = SystemUnit(
            functional_unit=new_fu,
            source="new_test_source",
            gene_families=gene_families,
            families_to_metainfo=fam2metainfo,
        )
        system.add_unit(new_unit)

        sources = system.annotation_sources()
        assert "test_source" in sources
        assert "new_test_source" in sources

    def test_modules_property(self, system, units):
        """Test modules property aggregating from all units."""
        for unit in units:
            system.add_unit(unit)
        modules = set(system.modules)
        expected_modules = {
            fam.module
            for unit in units
            for fam in unit.families
            if fam.module is not None
        }
        assert modules == expected_modules

    def test_get_module(self, system, units):
        """Test getting module"""
        for unit in units:
            system.add_unit(unit)
        for fam in system.families:
            if fam.module is not None:
                module = system.get_module(fam.module.ID)
                assert module == fam.module

    def test_get_module_nonexistent(self, system, units):
        """Test getting module not in system raise an error"""
        for unit in units:
            system.add_unit(unit)
        with pytest.raises(KeyError):
            system.get_module(12520)

    def test_spots_property(self, system, units):
        """Test spots property aggregating spots from all units."""

        expected_spots = set()
        for i, unit in enumerate(units):
            system.add_unit(unit)
            spot = Spot(spot_id=i)
            unit.add_spot(spot)
            expected_spots.add(spot)

        assert expected_spots == set(system.spots)

    def test_spots_property_with_family(self, system, units):
        """Test getting spot from system unit."""

        for unit in units:
            system.add_unit(unit)

        expected_spots = {
            spot for unit in units for fam in unit.families for spot in fam.spots
        }

        assert expected_spots == set(system.spots)

    def test_get_spots_with_family(self, system, units):
        """Test getting spots from system unit."""
        for unit in units:
            system.add_unit(unit)
            unit._spots_getter = {}

        for unit in system.units:
            for fam in unit.families:
                for spot in fam.spots:
                    assert system.get_spot(spot.ID) == spot

    def test_get_spot_nonexistent(self, system):
        """Test getting non-existent spot raises KeyError."""
        with pytest.raises(KeyError):
            system.get_spot(999)

    def test_regions_property(self, system, units):
        """Test regions property aggregating regions from all units."""
        expected_rgps = set()
        for i, unit in enumerate(units):
            system.add_unit(unit)
            rgp = Region(name=f"RGP_{i}")
            unit.add_region(rgp)
            expected_rgps.add(rgp)

        assert expected_rgps == set(system.regions)

    def test_get_regions(self, system, units):
        """Test getting regions from system unit."""
        for i, unit in enumerate(units):
            system.add_unit(unit)
            rgp = Region(name=f"RGP_{i}")
            unit.add_region(rgp)

        for unit in system.units:
            for rgp in unit.regions:
                assert system.get_region(rgp.name) == rgp

    def test_get_spots_with_region(self, system, units):
        """Test getting spots from system unit."""
        spots = {Spot(i) for i in range(5)}
        for unit in units:
            system.add_unit(unit)
            for i in range(randint(1, 10)):
                rgp = Region(name=f"RGP_{i}")
                unit.add_region(rgp)
                if randint(0, 1):
                    spot = choice(tuple(spots))
                    rgp.spot = spot
                    spot.add(rgp)
            unit._spots_getter = {}

        for spot in spots:
            assert system.get_spot(spot.ID) == spot
