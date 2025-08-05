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
        model = Model(name="test_model", min_mandatory=1, min_total=1)
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
            gf._genePerOrg[Organism(name=f"organism_{self._organims_counter}")] = {}
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
        assert unit._models_families is None
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

    def test_get_models_gene_families(self, system_unit, gene_families):
        """Test getting gene families."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)

        assert system_unit._get_models_families() == {gf1, gf2}

    def test_model_families_property(self, system_unit, gene_families):
        """Test model families property."""
        gf1, gf2 = gene_families
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        model_families = set(system_unit.models_families)
        assert model_families == {gf1, gf2}
        assert system_unit._models_families == model_families

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

    def test_models_organisms_property(self, system_unit, gene_families):
        """Test models organisms property."""
        gf1, gf2 = gene_families
        gf3 = GeneFamily(3, "new_family")
        gf3.add_metadata(Metadata("test_source", annot="test_func3"), 1)
        for org in set(gf1.organisms).union(set(gf2.organisms)):
            gf3._genePerOrg[org] = {}
        gf4 = GeneFamily(4, "context_family")
        organism = Organism("test_organism")
        gf4._genePerOrg[organism] = {}
        system_unit.add_family(gf1, "test_source", 1)
        system_unit.add_family(gf2, "test_source", 2)
        system_unit.add_family(gf3, "test_source", 1)
        system_unit.add_family(gf4)

        assert set(system_unit.models_organisms) == set(gf1.organisms).union(
            set(gf2.organisms)
        ).union(set(gf3.organisms))
        assert organism not in system_unit.models_organisms

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
