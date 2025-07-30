#!/usr/bin/env python3
# coding:utf-8

"""
Unit tests for the biological systems detection module.

This test suite provides comprehensive coverage for the classes and functions
defined in the biological systems detection module, including validation functions,
Models, FuncUnit, Family classes and their interactions.
"""
import re

import pytest
from pathlib import Path
from unittest.mock import mock_open, patch

# Import the module under test
# Assuming the module is named 'bio_systems' - adjust import as needed
from panorama.systems.models import (
    check_key,
    check_parameters,
    check_dict,
    Models,
    Model,
    FuncUnit,
    Family,
    _BasicFeatures,
    _FuFamFeatures,
    _ModFuFeatures,
)


class TestCheckKeys:
    """Test cases for keys validation functions."""

    required_keys = {"a", "b", "c"}

    def test_check_key_valid_keys(self):
        """Test check_key with all required keys present."""
        data = {"a": 1, "b": 2, "c": 3}
        # Should not raise any exception
        check_key(data, self.required_keys)

    def test_check_key_missing_keys(self):
        """Test check_key with missing required keys."""
        data = {"a": 1}

        with pytest.raises(KeyError) as exc_info:
            check_key(data, self.required_keys)

        assert "the following keys are required" in str(exc_info.value)
        assert "b" in str(exc_info.value)
        assert "c" in str(exc_info.value)

    def test_check_key_no_missing_keys_with_extra_keys_in_data(self):
        data = {"a": 1, "b": 2, "c": 3, "d": 4}
        try:
            check_key(data, self.required_keys)
        except KeyError:
            pytest.fail("check_key raised KeyError unexpectedly!")


class TestCheckParameters:
    """Test cases for parameters' validation utility functions."""

    mandatory_keys = {"min_mandatory", "min_total"}

    def test_check_parameters_valid_params(self):
        """Test check_parameters with valid parameter dictionary."""
        param_dict = {
            "min_mandatory": 2,
            "min_total": 5,
            "transitivity": 1,
            "window": 3,
        }

        # Should not raise any exception
        check_parameters(param_dict, self.mandatory_keys)

    def test_check_parameters_missing_mandatory(self):
        """Test check_parameters with missing mandatory keys."""
        param_dict = {"min_mandatory": 2}

        with pytest.raises(KeyError):
            check_parameters(param_dict, self.mandatory_keys)

    @pytest.mark.parametrize("invalid", [None, "invalid", 1.0, set()])
    def test_check_parameters_invalid_quorum_types(self, invalid):
        """Test check_parameters with invalid parameter types."""
        param_dict = {"min_mandatory": invalid, "min_total": 5}  # Should be int

        with pytest.raises(TypeError) as exc_info:
            check_parameters(param_dict, self.mandatory_keys)

        assert "min_mandatory must be an integer" in str(exc_info.value)

    def test_check_parameters_invalid_quorum_values(self):
        """Test check_parameters with invalid parameter values."""
        param_dict = {"min_mandatory": -2, "min_total": 5}  # Should be >= -1

        with pytest.raises(ValueError) as exc_info:
            check_parameters(param_dict, self.mandatory_keys)

        assert "min_mandatory must be >= -1" in str(exc_info.value)

    @pytest.mark.parametrize("invalid", [None, "invalid", 1.0, set()])
    def test_check_parameters_invalid_duplicate_types(self, invalid):
        """Test check_parameters with invalid parameter types."""
        param_dict = {
            "min_mandatory": invalid,
            "min_total": 5,
            "duplicate": invalid,
        }  # Should be int

        with pytest.raises(TypeError):
            check_parameters(param_dict, self.mandatory_keys)

    def test_check_parameters_invalid_duplicate_values(self):
        """Test check_parameters with invalid parameter values."""
        param_dict = {
            "min_mandatory": 1,
            "min_total": 5,
            "duplicate": -1,
        }  # Should be >= 0

        with pytest.raises(ValueError):
            check_parameters(param_dict, self.mandatory_keys)

    def test_check_parameters_unexpected_key(self):
        """Test check_parameters with unexpected parameter key."""
        param_dict = {"min_mandatory": 2, "min_total": 5, "unexpected_key": 1}

        with pytest.raises(KeyError) as exc_info:
            check_parameters(param_dict, self.mandatory_keys)

        assert "Unexpected parameter key: unexpected_key" in str(exc_info.value)


class TestCheckDict:
    """Test cases for model keys validation utility functions."""

    mandatory_keys = {"name", "presence"}
    param_keys = {"min_mandatory", "min_total"}

    def test_check_dict_valid_structure(self):
        """Test check_dict with a valid dictionary structure."""
        data_dict = {
            "name": "test_name",
            "presence": "mandatory",
            "parameters": {"min_mandatory": 1, "min_total": 2},
            "func_units": ["fu1", "fu2"],
            "families": ["fam1", "fam2"],
            "canonical": ["canonical1", "canonical2"],
            "exchangeable": ["fam1", "fam2", "fam3", "fam4"],
        }

        # Should not raise any exception
        check_dict(data_dict, self.mandatory_keys, self.param_keys)

    @pytest.mark.parametrize("invalid", [None, 1, 1.0, set()])
    def test_check_dict_invalid_name_type(self, invalid):
        """Test check_dict with an invalid name type."""
        data_dict = {"name": invalid, "presence": "mandatory"}

        with pytest.raises(TypeError):
            check_dict(data_dict, self.mandatory_keys)

    @pytest.mark.parametrize("invalid", [None, 1, 1.0, set()])
    def test_check_dict_invalid_presence_type(self, invalid):
        """Test check_dict with an invalid presence type."""
        data_dict = {"name": "test_name", "presence": invalid}

        with pytest.raises(TypeError):
            check_dict(data_dict, self.mandatory_keys)

    def test_check_dict_invalid_presence_value(self):
        """Test check_dict with an invalid presence value."""
        data_dict = {"name": "test_name", "presence": "invalid_presence"}

        with pytest.raises(ValueError):
            check_dict(data_dict, self.mandatory_keys)

    @pytest.mark.parametrize("field_name", ["func_units", "families"])
    @pytest.mark.parametrize("invalid", [None, "invalid", 1, 1.0, set()])
    def test_check_dict_invalid_fu_fam_type(self, field_name, invalid):
        """Test check_dict with an invalid functional unit or family type."""
        data_dict = {"name": "test_name", "presence": "mandatory", field_name: invalid}

        with pytest.raises(TypeError):
            check_dict(data_dict, self.mandatory_keys)

    @pytest.mark.parametrize("field_name", ["func_units", "families"])
    def test_check_dict_invalid_fu_fam_value(self, field_name):
        """Test check_dict with an invalid functional unit or family type."""
        data_dict = {"name": "test_name", "presence": "mandatory", field_name: []}

        with pytest.raises(ValueError):
            check_dict(data_dict, self.mandatory_keys)

    @pytest.mark.parametrize("field_name", ["canonical", "exchangeable"])
    @pytest.mark.parametrize("invalid", [None, "invalid", 1, 1.0, set()])
    def test_check_dict_invalid_canon_exc_type(self, field_name, invalid):
        """Test check_dict with an invalid functional unit or family type."""
        data_dict = {"name": "test_name", "presence": "mandatory", field_name: invalid}

        with pytest.raises(TypeError):
            check_dict(data_dict, self.mandatory_keys)

    @pytest.mark.parametrize("invalid", [[1], ["a", 1], [1, "a"]])
    def test_check_dict_invalid_exchangeable_value(self, invalid):
        """Test check_dict with an invalid exchangeable value."""
        data_dict = {
            "name": "test_name",
            "presence": "mandatory",
            "exchangeable": invalid,
        }

        with pytest.raises(ValueError):
            check_dict(data_dict, self.mandatory_keys)

    @pytest.mark.parametrize("invalid", [None, "invalid", 1, 1.0, set()])
    def test_check_dict_invalid_same_strand_type(self, invalid):
        """Test check_dict with an invalid same_strand type."""
        data_dict = {
            "name": "test_name",
            "presence": "mandatory",
            "same_strand": invalid,
        }
        with pytest.raises(TypeError):
            check_dict(data_dict, self.mandatory_keys)

    def test_check_dict_unexpected_key(self):
        """Test check_dict with unexpected parameter key."""
        data_dict = {"name": "test_name", "presence": "mandatory", "unexpected_key": 1}

        with pytest.raises(KeyError):
            check_dict(data_dict, self.mandatory_keys)


class TestBasicFeatures:
    """Test cases for _BasicFeatures class."""

    class MockParent:
        """Artificial parent class for _BasicFeatures."""

        def __init__(self):
            self.transitivity = 5
            self.window = 10

    @pytest.fixture()
    def parent(self):
        """Fixture for creating a mock parent object."""
        return self.MockParent()

    def test_initialization_default_values(self):
        """Test _BasicFeatures initialization with default values."""
        bf = _BasicFeatures()

        assert bf.name == ""
        assert bf.transitivity == 0
        assert bf.window == 1

    def test_initialization_custom_values(self):
        """Test _BasicFeatures initialization with custom values."""
        bf = _BasicFeatures(name="test", transitivity=2, window=5)

        assert bf.name == "test"
        assert bf.transitivity == 2
        assert bf.window == 5

    def test_string_representations(self):
        """Test __repr__ and __str__ methods."""
        bf = _BasicFeatures(name="test_feature")

        expected = "_BasicFeatures name: test_feature"
        assert repr(bf) == expected
        assert str(bf) == expected

    def test_read_parameters_with_values(self):
        """Test read_parameters method with parameter values present."""
        bf = _BasicFeatures()
        parameters = {"transitivity": 3, "window": 7}
        param_keys = {"transitivity", "window"}

        bf.read_parameters(parameters, param_keys)

        assert bf.transitivity == 3
        assert bf.window == 7

    def test_read_parameters_with_parent(self, parent):
        """Test read_parameters method with parent fallback."""

        bf = _BasicFeatures()
        bf._parent = parent
        parameters = {}
        param_keys = {"transitivity", "window"}

        bf.read_parameters(parameters, param_keys)

        assert bf.transitivity == 5
        assert bf.window == 10

    def test_read_parameters_with_parent_and_values(self, parent):
        """Test read_parameters method with parent fallback and parameter values present."""
        bf = _BasicFeatures()
        bf._parent = parent
        parameters = {"transitivity": 3}
        param_keys = {"transitivity", "window"}

        bf.read_parameters(parameters, param_keys)

        assert bf.transitivity == 3
        assert bf.window == 10


class TestFuFamFeatures:
    """Test cases for _FuFamFeatures class."""

    def test_initialization_default_values(self):
        """Test _FuFamFeatures initialization with default values."""
        fff = _FuFamFeatures()

        assert fff.presence == ""
        assert fff.duplicate == 0
        assert fff.exchangeable == set()
        assert fff.multi_system is False
        assert fff.multi_model is True
        assert fff._parent is None

    def test_initialization_custom_values(self):
        """Test _FuFamFeatures initialization with custom values."""
        exchangeable_set = {"fam1", "fam2"}
        fff = _FuFamFeatures(
            presence="mandatory",
            duplicate=2,
            exchangeable=exchangeable_set,
            multi_system=True,
            multi_model=False,
        )

        assert fff.presence == "mandatory"
        assert fff.duplicate == 2
        assert fff.exchangeable == exchangeable_set
        assert fff.multi_system is True
        assert fff.multi_model is False


class TestModFuFeatures:
    """Test cases for _ModFuFeatures class."""

    class MockChild:
        """Artificial child class for _ModFuFeatures."""

        def __init__(self, name, presence="", duplicate=0):
            self.name = name
            self.presence = presence
            self.duplicate = duplicate

    @pytest.fixture()
    def child(self):
        """Fixture for creating a mock child object."""
        return self.MockChild("test_child")

    @pytest.fixture()
    def children(self):
        """Fixture for creating a list of mock child objects."""
        yield [
            self.MockChild("child1"),
            self.MockChild("child2"),
            self.MockChild("child3"),
        ]

    def test_initialization_default_values(self):
        """Test _ModFuFeatures initialization with default values."""
        mff = _ModFuFeatures()

        assert mff.mandatory == set()
        assert mff.accessory == set()
        assert mff.forbidden == set()
        assert mff.neutral == set()
        assert mff.min_mandatory == 1
        assert mff.min_total == 1
        assert mff.same_strand is False

    def test_children_property(self, children):
        """Test _children property returns all child elements."""
        child1, child2, child3 = children

        mff = _ModFuFeatures(mandatory={child1}, accessory={child2}, forbidden={child3})

        children_list = list(mff._children)
        assert len(children_list) == 3
        assert child1 in children_list
        assert child2 in children_list
        assert child3 in children_list

    def test_child_names_all(self, children):
        """Test _child_names method returns all child names."""
        child1, child2, child3 = children

        mff = _ModFuFeatures(mandatory={child1}, accessory={child2}, neutral={child3})

        names = mff._child_names()
        assert names == {"child1", "child2", "child3"}

    def test_child_names_filtered(self, child):
        """Test _child_names method with presence filter."""

        mff = _ModFuFeatures(mandatory={child})
        child.presence = "mandatory"
        mandatory_names = mff._child_names(presence="mandatory")
        assert mandatory_names == {"test_child"}

    def test_get_child_type(self, child):
        """Tests the retrieval of the child type from the `_ModFuFeatures` class."""
        mff = _ModFuFeatures(mandatory={child})
        child.presence = "mandatory"

        assert mff.child_type == type(child).__name__

    def test_inconsistent_child_type(self, child):
        """Test that an TypeError is raised when the child type is inconsistent."""

        class DumbChild:
            """Dumb child class for testing."""

            presence = "mandatory"

        dumb_child = DumbChild()
        mff = _ModFuFeatures(mandatory={child, dumb_child})
        child.presence = "mandatory"

        with pytest.raises(
            TypeError,
            match=f"The child type is inconsistent. It contains {type(child).__name__} and {type(dumb_child).__name__}",
        ):
            _ = mff.child_type

    def test_check_validation_success(self, children):
        """Test _check method with valid configuration."""
        child1, child2, _ = children
        child1.presence = "mandatory"
        child2.presence = "accessory"

        mff = _ModFuFeatures(
            mandatory={child1}, accessory={child2}, min_mandatory=1, min_total=2
        )

        # Should not raise any exception
        mff._check()

    def test_check_validation_failure_no_mandatory(self, child):
        """Test _check method fails with no mandatory children."""
        mff = _ModFuFeatures(min_mandatory=1, min_total=1)

        with pytest.raises(
            Exception,
            match=f"There are no mandatory {mff.child_type}. "
            f"You should have at least one mandatory {mff.child_type} with mandatory presence.",
        ):
            mff._check()

        mff.accessory = {child}
        with pytest.raises(
            Exception,
            match=f"There are no mandatory {mff.child_type}. "
            f"You should have at least one mandatory {mff.child_type} with mandatory presence.",
        ):
            mff._check()

    def test_check_min_mandatory_less_than_required(self, child):
        """Test exception when mandatory elements are less than the required minimum."""
        mff = _ModFuFeatures(mandatory={child}, min_mandatory=2, min_total=2)
        with pytest.raises(
            Exception,
            match=f"There are less mandatory {mff.child_type} than the minimum mandatory",
        ):
            mff._check()

    def test_check_min_mandatory_success_with_duplicate(self, child):
        """Test success when mandatory elements are less than the required minimum but are duplicated."""
        mff = _ModFuFeatures(mandatory={child}, min_mandatory=2, min_total=2)
        child.duplicate = 1
        mff._check()

    def test_check_min_mandatory_less_than_required_with_duplicate(self, child):
        """Test exception when mandatory elements are less than the required minimum."""

        mff = _ModFuFeatures(mandatory={child}, min_mandatory=3, min_total=3)
        child.duplicate = 1
        with pytest.raises(
            Exception,
            match=f"There are less mandatory {mff.child_type} than the minimum mandatory",
        ):
            mff._check()

    def test_check_min_total_less_than_required(self, children):
        """Test exception when mandatory elements are less than the required minimum."""
        child1, child2, child3 = children
        mff = _ModFuFeatures(
            mandatory={child1},
            accessory={child2},
            neutral={child3},
            min_mandatory=1,
            min_total=3,
        )
        with pytest.raises(
            Exception,
            match=f"There are less {mff.child_type} than the minimum total",
        ):
            mff._check()

    def test_check_min_total_success_with_duplicate(self, children):
        """Test success when mandatory elements are less than the required minimum but are duplicated."""
        child1, child2, child3 = children
        mff = _ModFuFeatures(
            mandatory={child1},
            accessory={child2},
            neutral={child3},
            min_mandatory=1,
            min_total=4,
        )
        child1.duplicate = 1
        child2.duplicate = 1
        mff._check()

    def test_check_min_total_less_than_required_with_duplicate(self, children):
        """Test exception when mandatory elements are less than the required minimum."""

        child1, child2, child3 = children
        mff = _ModFuFeatures(
            mandatory={child1},
            accessory={child2},
            neutral={child3},
            min_mandatory=1,
            min_total=5,
        )
        child1.duplicate = 1
        child2.duplicate = 1
        child3.duplicate = 2
        with pytest.raises(
            Exception,
            match=f"There are less {mff.child_type} than the minimum total",
        ):
            mff._check()

    def test_check_min_mandatory_greater_than_min_total(self, child):
        """
        Test exception when minimum mandatory count is greater than minimum total count.
        """
        mff = _ModFuFeatures(mandatory={child}, min_mandatory=3, min_total=1)
        child.presence = "mandatory"
        child.duplicate = 2
        with pytest.raises(
            Exception,
            match=f"Minimum mandatory {mff.child_type} value is greater than minimum total",
        ):
            mff._check()

    @pytest.mark.parametrize(
        "presence", ["mandatory", "accessory", "forbidden", "neutral"]
    )
    def test_add_child_success(self, child, presence):
        """Test adding a child with mandatory presence."""

        child.presence = presence
        mff = _ModFuFeatures()

        mff.add(child)

        assert child in mff.__getattribute__(presence)

    def test_add_child_failure_inconsistent_child_type(self, child):
        """Test that an TypeError is raised when the child type is inconsistent."""

        class DumbChild:
            """Dumb child class for testing."""

            presence = "mandatory"

        dumb_child = DumbChild()
        mff = _ModFuFeatures(mandatory={child})
        with pytest.raises(
            TypeError,
            match=f"The child type is inconsistent. Expected {type(child).__name__} but found {type(dumb_child).__name__}",
        ):
            mff.add(dumb_child)

    def test_add_child_failure_unexepected_presence_name(self, child):
        """Test that an exception is raised when the child type is inconsistent."""
        mff = _ModFuFeatures()
        child.presence = "unexpected"
        with pytest.raises(
            ValueError,
            match=re.escape(
                f"The child {child.name} does not have a valid presence attribute ({child.presence})."
            ),
        ):
            mff.add(child)

    def test_mk_child_getter_success(self, children):
        """Test that the child getter returns the expected dictionary."""
        mff = _ModFuFeatures(mandatory=set(children))
        mff._mk_child_getter()

        assert mff._child_getter == {child.name: child for child in children}

    def test_get_child_success(self, children):
        child1, child2, child3 = children
        mff = _ModFuFeatures(mandatory=set(children))

        assert mff.get("child1") == child1
        assert mff.get("child2") == child2
        assert mff.get("child3") == child3

    def test_get_child_failure(self, children):
        mff = _ModFuFeatures(mandatory=set(children))
        name = "notin_child"
        with pytest.raises(
            KeyError, match=f"No such {mff.child_type} with name {name} in {type(mff)}"
        ):
            mff.get(name)


class TestFamily:
    """Test cases for Family class."""

    class MockModel:
        """Represents a mock model class."""

        pass

    class MockFuncUnit:
        """Represents a mock func_unit class."""

        def __init__(self, name, presence):
            self.name = name
            self.presence = presence
            self.families = []

    @pytest.fixture
    def model(self):
        """Creates a fixture that yields a mock model instance."""
        yield self.MockModel()

    @pytest.fixture
    def func_unit(self, model):
        """Creates a fixture that yields a mock func_unit instance."""
        func_unit = self.MockFuncUnit("test_fu", "mandatory")
        func_unit.model = model
        yield func_unit

    def test_initialization_default_values(self):
        """Test Family initialization with default values."""
        family = Family()

        assert family.name == ""
        assert family.presence == ""
        assert family.transitivity == 0
        assert family.window == 1
        assert family.duplicate == 0
        assert family.exchangeable == set()
        assert family._parent is None
        assert family.multi_system is False
        assert family.multi_model is False

    def test_initialization_custom_values(self):
        """Test Family initialization with custom values."""
        family = Family(
            name="test_family",
            presence="mandatory",
            transitivity=2,
            window=5,
            duplicate=1,
        )

        assert family.name == "test_family"
        assert family.presence == "mandatory"
        assert family.transitivity == 2
        assert family.window == 5
        assert family.duplicate == 1

    def test_func_unit_property(self, func_unit):
        """Test func_unit property getter and setter."""
        family = Family()

        family.func_unit = func_unit
        assert family._parent == func_unit
        assert family.func_unit == func_unit

        del family.func_unit
        assert family.func_unit is None

    def test_model_property(self, func_unit, model):
        """Test model property returns func_unit's model."""
        family = Family()
        family.func_unit = func_unit

        assert family.model == model

    def test_read_method(self):
        """Test read method processes family data correctly."""
        family = Family()
        data_fam = {
            "name": "test_family",
            "presence": "mandatory",
            "parameters": {"transitivity": 2, "duplicate": 1},
            "exchangeable": ["fam1", "fam2"],
        }

        family.read(data_fam)

        assert family.name == "test_family"
        assert family.presence == "mandatory"
        assert family.transitivity == 2
        assert family.duplicate == 1
        assert family.exchangeable == {"fam1", "fam2"}


class TestFuncUnit:
    """Test cases for FuncUnit class."""

    class MockModel:
        """Mock model class for testing."""

        def __init__(self, name=""):
            self.name = name

    @pytest.fixture
    def model(self):
        """Creates a fixture that yields a mock model instance."""
        yield self.MockModel("test_model")

    @pytest.fixture
    def family(self):
        """Creates a fixture that yields a family instance."""
        yield Family(name="test_family", presence="mandatory")

    @pytest.fixture
    def families(self):
        """Creates a fixture that yields a list of families."""
        yield [
            Family(name="fam1", presence="mandatory"),
            Family(name="fam2", presence="accessory"),
        ]

    def test_initialization_default_values(self):
        """Test FuncUnit initialization with default values."""
        func_unit = FuncUnit()

        assert func_unit.name == ""
        assert func_unit.presence == ""
        assert func_unit.mandatory == set()
        assert func_unit.accessory == set()
        assert func_unit.forbidden == set()
        assert func_unit.neutral == set()
        assert func_unit.min_mandatory == 1
        assert func_unit.min_total == 1
        assert func_unit.transitivity == 0
        assert func_unit.window == 1
        assert func_unit.duplicate == 0
        assert func_unit.same_strand is False
        assert func_unit._parent is None
        assert func_unit.exchangeable == set()
        assert func_unit.multi_system is False
        assert func_unit.multi_model is False

    def test_initialization_custom_values(self, families):
        """Test FuncUnit initialization with custom values."""
        family1, family2 = families
        func_unit = FuncUnit(
            name="test_fu",
            presence="mandatory",
            mandatory={family1},
            accessory={family2},
            min_mandatory=1,
            min_total=3,
            transitivity=2,
            window=5,
            duplicate=1,
            same_strand=True,
        )

        assert func_unit.name == "test_fu"
        assert func_unit.presence == "mandatory"
        assert func_unit.min_mandatory == 1
        assert func_unit.min_total == 3
        assert func_unit.transitivity == 2
        assert func_unit.window == 5
        assert func_unit.duplicate == 1
        assert func_unit.same_strand is True

    def test_model_property(self, model):
        """Test model property getter and setter."""
        func_unit = FuncUnit()

        func_unit.model = model
        assert func_unit._parent == model
        assert func_unit.model == model

        del func_unit.model
        assert func_unit.model is None

    def test_families_property(self, families):
        """Test families property returns child families."""
        family1, family2 = families
        func_unit = FuncUnit(mandatory={family1}, accessory={family2})

        families_list = set(func_unit.families)
        assert families_list == set(families)

    def test_size_property(self, families):
        """Test size property returns number of families."""
        family1, family2 = families
        func_unit = FuncUnit(mandatory={family1}, accessory={family2})

        assert func_unit.size == 2

    def test_add_family(self, family):
        """Test add_family method adds a family to the FuncUnit."""
        func_unit = FuncUnit()
        func_unit.add(family)

        assert list(func_unit.families) == [family]

    def test_get_families(self, family):
        """Test get_families method returns the list of child fus."""
        func_unit = FuncUnit()
        func_unit.add(family)

        assert func_unit.get("test_family") == family

    @pytest.mark.parametrize(
        "presence", ["mandatory", "accessory", "forbidden", "neutral"]
    )
    def test_list_families_names(self, families, presence):
        """Test list_families method returns the list of child fus."""
        func_unit = FuncUnit()
        fam1, fam2 = families
        func_unit.add(fam1)
        func_unit.add(fam2)
        fam3 = Family(name="test_family", presence=presence)
        func_unit.add(fam3)

        assert func_unit.families_names() == {fam1.name, fam2.name, fam3.name}
        if fam1.presence == fam3.presence:
            assert func_unit.families_names(fam3.presence) == {fam1.name, fam3.name}
        elif fam2.presence == fam3.presence:
            assert func_unit.families_names(fam3.presence) == {fam2.name, fam3.name}
        else:
            assert func_unit.families_names(fam3.presence) == {fam3.name}

    @pytest.mark.parametrize(
        "presence", ["mandatory", "accessory", "forbidden", "neutral"]
    )
    def test_get_duplicate_families(self, families, presence):
        """Test get_duplicate_families method returns the list of child fus with duplicate data."""
        func_unit = FuncUnit()
        fam1, fam2 = families
        func_unit.add(fam1)
        func_unit.add(fam2)
        fam3 = Family(name="duplicated_fam", presence=presence, duplicate=1)
        func_unit.add(fam3)

        assert set(func_unit.duplicate_fam()) == {fam3}
        assert set(func_unit.duplicate_fam(presence)) == {fam3}

    def test_read_method(self):
        """Test read method processes functional unit data correctly."""
        func_unit = FuncUnit()
        data_fu = {
            "name": "test_fu",
            "presence": "mandatory",
            "families": [
                {"name": "fam1", "presence": "mandatory"},
                {"name": "fam2", "presence": "accessory"},
            ],
            "parameters": {"duplicate": 1, "min_total": 2},
        }

        func_unit.read(data_fu)

        assert func_unit.name == "test_fu"
        assert func_unit.presence == "mandatory"
        mandatory_fam = list(func_unit.mandatory)[0]
        assert mandatory_fam.name == "fam1"
        accessory_fam = list(func_unit.accessory)[0]
        assert accessory_fam.name == "fam2"
        assert func_unit.duplicate == 1
        assert func_unit.min_total == 2
        assert func_unit.size == 2  # Two families added


class TestModel:
    """Test cases for Model class."""

    @pytest.fixture
    def func_unit(self):
        """Creates a fixture that yields a FuncUnit instance."""
        func_unit = FuncUnit(name="test_fu", presence="mandatory")
        yield func_unit

    @pytest.fixture
    def families(self):
        """Creates a fixture that yields a list of families."""
        family1 = Family(name="fam1", presence="mandatory")
        family2 = Family(name="fam2", presence="accessory")
        family3 = Family(name="fam3", presence="mandatory")
        family4 = Family(name="fam4", presence="mandatory")
        family5 = Family(name="fam5", presence="accessory")
        family6 = Family(name="fam6", presence="mandatory")
        family7 = Family(name="fam7", presence="accessory")
        family8 = Family(name="fam8", presence="accessory")
        family9 = Family(name="fam9", presence="mandatory")
        family10 = Family(name="fam10", presence="accessory")
        yield [
            family1,
            family2,
            family3,
            family4,
            family5,
            family6,
            family7,
            family8,
            family9,
            family10,
        ]

    @pytest.fixture
    def func_units(self, families):
        """Creates a fixture that yields a list of FuncUnit instances."""
        fam1, fam2, fam3, fam4, fam5, fam6, fam7, fam8, fam9, fam10 = families
        func_unit1 = FuncUnit(
            name="test_fu1", presence="mandatory", mandatory={fam1}, accessory={fam2}
        )
        func_unit2 = FuncUnit(name="test_fu2", presence="mandatory", mandatory={fam3})
        func_unit3 = FuncUnit(
            name="test_fu3", presence="forbidden", mandatory={fam4}, accessory={fam5}
        )
        func_unit4 = FuncUnit(
            name="test_fu4",
            presence="accessory",
            mandatory={fam6},
            accessory={fam7, fam8},
        )
        func_unit5 = FuncUnit(
            name="test_fu5", presence="neutral", mandatory={fam9}, accessory={fam10}
        )

        yield [func_unit1, func_unit2, func_unit3, func_unit4, func_unit5]

    def test_initialization_default_values(self):
        """Test Model initialization with default values."""
        model = Model()

        assert model.name == ""
        assert model.mandatory == set()
        assert model.accessory == set()
        assert model.forbidden == set()
        assert model.neutral == set()
        assert model.min_mandatory == 1
        assert model.min_total == 1
        assert model.transitivity == 0
        assert model.window == 1
        assert model.same_strand is False
        assert model.canonical == []

    def test_initialization_custom_value(self):
        """Test Model initialization with custom values."""
        func_unit1 = FuncUnit(name="fun1")
        func_unit2 = FuncUnit(name="fun2")
        func_unit3 = FuncUnit(name="fun3")
        func_unit4 = FuncUnit(name="fun4")
        model = Model(
            name="test_model",
            mandatory={func_unit1},
            accessory={func_unit2},
            forbidden={func_unit3},
            neutral={func_unit4},
            min_mandatory=1,
            min_total=3,
            transitivity=2,
            window=5,
            same_strand=True,
            canonical=["canonical_model"],
        )

        assert model.name == "test_model"
        assert model.mandatory == {func_unit1}
        assert model.accessory == {func_unit2}
        assert model.forbidden == {func_unit3}
        assert model.neutral == {func_unit4}
        assert model.min_mandatory == 1
        assert model.min_total == 3
        assert model.transitivity == 2
        assert model.window == 5
        assert model.same_strand is True
        assert model.canonical == ["canonical_model"]

    def test_func_units_property(self, func_units):
        """Test func_units property returns child functional units."""
        model = Model()
        fu1, fu2, fu3, fu4, fu5 = func_units

        model.add(fu1)
        model.add(fu2)
        model.add(fu3)
        model.add(fu4)
        model.add(fu5)

        func_units_set = set(model.func_units)
        assert func_units_set == set(func_units)

    @pytest.mark.parametrize(
        "presence", [None, "mandatory", "accessory", "forbidden", "neutral"]
    )
    def test_get_func_units_names(self, func_units, presence):
        """Test get_func_units_names method returns the list of functional unit names."""
        model = Model()
        fu1, fu2, fu3, fu4, fu5 = func_units

        model.add(fu1)
        model.add(fu2)
        model.add(fu3)
        model.add(fu4)
        model.add(fu5)

        func_units_names = set(model.func_units_names(presence))
        assert func_units_names == {
            fu.name for fu in func_units if fu.presence == presence or presence is None
        }

    def test_families_property(self, families, func_units):
        """Test families property returns all families from functional units."""
        model = Model()
        fu1, fu2, fu3, fu4, fu5 = func_units
        fam1, fam2, fam3, fam4, fam5, fam6, fam7, fam8, fam9, fam10 = families
        model.add(fu1)
        model.add(fu2)
        model.add(fu3)
        model.add(fu4)
        model.add(fu5)

        families_set = set(model.families)
        assert families_set == {
            fam1,
            fam2,
            fam3,
            fam4,
            fam5,
            fam6,
            fam7,
            fam8,
            fam9,
            fam10,
        }

    def test_size_property(self, func_units, families):
        """Test size property returns tuple of (func_units_count, families_count)."""
        model = Model()
        fu1, fu2, fu3, fu4, fu5 = func_units
        fam1, fam2, fam3, fam4, fam5, fam6, fam7, fam8, fam9, fam10 = families
        model.add(fu1)
        model.add(fu2)
        model.add(fu3)
        model.add(fu4)
        model.add(fu5)

        assert model.size == (len(set(func_units)), len(set(families)))  # 1 func_unit, 2 families

    def test_get_duplicate_fu(self, func_units):
        """Test get_duplicate_fu method returns the list of functional units with duplicate data."""
        model = Model()
        fu1, fu2, fu3, fu4, fu5 = func_units
        fu1.duplicate = 1
        fu3.duplicate = 2
        model.add(fu1)
        model.add(fu2)
        model.add(fu3)
        model.add(fu4)
        model.add(fu5)

        assert set(model.duplicate_fu()) == {fu1, fu3}
        assert set(model.duplicate_fu("mandatory")) == {fu1}

    def test_get_functional_units(self, func_units):
        """Test get_functional_units method returns the list of functional units."""
        model = Model()
        fu1, fu2, fu3, fu4, fu5 = func_units
        model.add(fu1)
        model.add(fu2)
        model.add(fu3)
        model.add(fu4)
        model.add(fu5)
        assert model.get("test_fu1") == fu1
        assert model.get("test_fu2") == fu2
        assert model.get("test_fu3") == fu3
        assert model.get("test_fu4") == fu4
        assert model.get("test_fu5") == fu5

    def test_read_method(self):
        """Test read method processes model data correctly."""
        model = Model()
        data_model = {
            "name": "test_model",
            "parameters": {"transitivity": 2, "min_mandatory": 1, "min_total": 1},
            "func_units": [
                {
                    "name": "fu1",
                    "presence": "mandatory",
                    "families": [{"name": "fam1", "presence": "mandatory"}],
                    "parameters": {}
                }
            ],
            "canonical": ["element1", "element2"],
        }

        model.read(data_model)

        assert model.name == "test_model"
        assert model.transitivity == 2
        assert model.min_mandatory == 1
        assert model.min_total == 1
        assert model.canonical == ["element1", "element2"]
        assert model.size == (1, 1)  # 1 func_unit, 1 family

    def test_read_model_static_method(self):
        """Test read_model static method creates and returns Model instance."""
        data_model = {
            "name": "test_model",
            "parameters": {"transitivity": 2, "min_mandatory": 1, "min_total": 1},
            "func_units": [
                {
                    "name": "fu1",
                    "presence": "mandatory",
                    "families": [{"name": "fam1", "presence": "mandatory"}],
                    "parameters": {}
                }
            ],
            "canonical": ["element1", "element2"],
        }

        model = Model.read_model(data_model)

        assert isinstance(model, Model)
        assert model.name == "test_model"
        assert model.transitivity == 2
        assert model.min_mandatory == 1
        assert model.min_total == 1
        assert model.canonical == ["element1", "element2"]
        assert model.size == (1, 1)  # 1 func_unit, 1 family


class TestModels:
    """Test cases for Models class."""

    def test_initialization_default(self):
        """Test Models initialization with default values."""
        models = Models()

        assert models.size == 0
        assert list(models.value) == []

    def test_initialization_with_models(self):
        """Test Models initialization with provided models."""
        model1 = Model(name="model1")
        model2 = Model(name="model2")
        model_dict = {"model1": model1, "model2": model2}

        models = Models()
        models._model_getter = model_dict

        assert models.size == 2

        model_names = [model.name for model in models.value]
        assert "model1" in model_names
        assert "model2" in model_names

    def test_add_model_success(self):
        """Test adding a new model successfully."""
        models = Models()
        model = Model(name="test_model")
        # Add a mandatory functional unit to make the model valid
        fu = FuncUnit(name="test_fu", presence="mandatory")
        model.add(fu)

        models.add_model(model)

        assert models.size == 1
        assert models.get_model("test_model") == model

    def test_add_model_duplicate_name(self):
        """Test adding a model with duplicate name raises exception."""
        models = Models()
        model1 = Model(name="duplicate_name")
        model2 = Model(name="duplicate_name")

        # Add mandatory functional units to make models valid
        fu1 = FuncUnit(name="fu1", presence="mandatory")
        fu2 = FuncUnit(name="fu2", presence="mandatory")
        model1.add(fu1)
        model2.add(fu2)

        models.add_model(model1)

        with pytest.raises(Exception) as exc_info:
            models.add_model(model2)

        assert "already in set of value" in str(exc_info.value)

    def test_get_model_success(self):
        """Test getting an existing model by name."""
        models = Models()
        model = Model(name="get_test_model")
        fu = FuncUnit(name="test_fu", presence="mandatory")
        model.add(fu)

        models.add_model(model)
        retrieved_model = models.get_model("get_test_model")

        assert retrieved_model == model

    def test_get_model_not_found(self):
        """Test getting a non-existent model raises KeyError."""
        models = Models()

        with pytest.raises(KeyError) as exc_info:
            models.get_model("non_existent_model")

        assert "Model not present in set of value" in str(exc_info.value)

    @patch("builtins.open", new_callable=mock_open)
    @patch("json.load")
    def test_read_method_success(self, mock_json_load, mock_file):
        """Test read method successfully loads model from JSON file."""
        # Mock JSON data
        mock_json_data = {
            "name": "file_test_model",
            "parameters": {"transitivity": 1, "min_mandatory": 1, "min_total": 1},
            "func_units": [
                {
                    "name": "fu1",
                    "presence": "mandatory",
                    "families": [{"name": "fam1", "presence": "mandatory"}],
                }
            ],
        }
        mock_json_load.return_value = mock_json_data

        models = Models()
        test_path = Path("test_model.json")

        models.read(test_path)

        assert models.size == 1
        model = models.get_model("file_test_model")
        assert model.name == "file_test_model"

    def test_func_units_property(self):
        """Test func_units property returns all functional units from all models."""
        models = Models()

        # Create first model with functional unit
        model1 = Model(name="model1")
        fu1 = FuncUnit(name="fu1", presence="mandatory")
        model1.add(fu1)

        # Create second model with functional unit
        model2 = Model(name="model2")
        fu2 = FuncUnit(name="fu2", presence="mandatory")
        model2.add(fu2)

        models.add_model(model1)
        models.add_model(model2)

        func_units_list = list(models.func_units)
        func_unit_names = [fu.name for fu in func_units_list]

        assert len(func_units_list) == 2
        assert "fu1" in func_unit_names
        assert "fu2" in func_unit_names

    def test_families_property(self):
        """Test families property returns all families from all models."""
        models = Models()

        # Create model with functional unit and families
        model = Model(name="model1")
        fu = FuncUnit(name="fu1", presence="mandatory")
        family1 = Family(name="fam1", presence="mandatory")
        family2 = Family(name="fam2", presence="accessory")

        fu.add(family1)
        fu.add(family2)
        model.add(fu)
        models.add_model(model)

        families_list = list(models.families)
        family_names = [fam.name for fam in families_list]

        assert len(families_list) == 2
        assert "fam1" in family_names
        assert "fam2" in family_names


class TestIntegration:
    """Integration tests for the complete system."""

    def test_complete_model_creation_and_validation(self):
        """Test creating a complete model hierarchy and validating it."""
        # Create a complete model with functional units and families
        model = Model(
            name="integration_test_model",
            min_mandatory=1,
            min_total=2,
            transitivity=1,
            window=2,
        )

        # Create functional unit
        func_unit = FuncUnit(
            name="test_functional_unit",
            presence="mandatory",
            min_mandatory=1,
            min_total=1,
        )
        func_unit.model = model

        # Create families
        family1 = Family(name="family1", presence="mandatory", duplicate=0)
        family1.func_unit = func_unit

        family2 = Family(
            name="family2",
            presence="accessory",
            duplicate=1,
            exchangeable={"alt_family"},
        )
        family2.func_unit = func_unit

        # Build the hierarchy
        func_unit.add(family1)
        func_unit.add(family2)
        model.add(func_unit)

        # Validate the complete model
        model.check_model()  # Should not raise exception

        # Test properties and relationships
        assert model.size == (1, 2)  # 1 func_unit, 2 families
        assert func_unit.size == 2  # 2 families
        assert family1.model == model
        assert family2.model == model
        assert func_unit.model == model

        # Test retrieving elements
        retrieved_fu = model.get("test_functional_unit")
        assert retrieved_fu == func_unit

        retrieved_family = func_unit.get("family1")
        assert retrieved_family == family1

    def test_models_collection_management(self):
        """Test managing multiple models in Models collection."""
        models = Models()

        # Create multiple models
        for i in range(3):
            model = Model(name=f"model_{i}")
            func_unit = FuncUnit(name=f"fu_{i}", presence="mandatory")
            family = Family(name=f"fam_{i}", presence="mandatory")

            func_unit.model = model
            family.func_unit = func_unit

            func_unit.add(family)
            model.add(func_unit)
            models.add_model(model)

        # Test collection properties
        assert models.size == 3

        func_units_list = list(models.func_units)
        families_list = list(models.families)

        assert len(func_units_list) == 3
        assert len(families_list) == 3

        # Test mappings
        fu_to_model_map = models.func_units_to_model()
        fam_to_model_map = models.families_to_model()

        assert len(fu_to_model_map) == 3
        assert len(fam_to_model_map) == 3

        # Verify mappings are correct
        for fu in func_units_list:
            assert fu_to_model_map[fu] == fu.model

        for fam in families_list:
            assert fam_to_model_map[fam] == fam.model
