#! /usr/bin/env python3

import pytest
from unittest import mock
import pickle
from random import randint
from typing import Set, List

from panorama.systems.models import _BasicFeatures, _FuFamFeatures, _ModFuFeatures, Model, FuncUnit, Family


class TestBasicFeatures:
    @pytest.fixture
    def basic(self):
        """Fixture to initialize a _BasicFeatures object.

        Returns:
            _BasicFeatures: A _BasicFeatures object initialized with test data.
        """
        return _BasicFeatures(name="Test", transitivity=2, window=3)

    def test_initialization(self, basic):
        """Test the initialization of _BasicFeatures attributes.

        Args:
            basic (_BasicFeatures): Fixture providing a basic _BasicFeatures object.
        """
        assert basic.name == "Test"
        assert basic.transitivity == 2
        assert basic.window == 3

    def test_repr(self, basic):
        """Test the __repr__ method of _BasicFeatures.

        Args:
            basic (_BasicFeatures): Fixture providing a basic _BasicFeatures object.
        """
        assert repr(basic) == "_BasicFeatures name: Test"

    def test_str(self, basic):
        """Test the __str__ method of _BasicFeatures.

        Args:
            basic (_BasicFeatures): Fixture providing a basic _BasicFeatures object.
        """
        assert str(basic) == "_BasicFeatures name: Test"


class TestFuFamFeatures:
    @pytest.fixture
    def fufam(self):
        """Fixture to initialize a _FuFamFeatures object.

        Returns:
            _FuFamFeatures: A _FuFamFeatures object initialized with test data.
        """
        return _FuFamFeatures(presence="mandatory", duplicate=2, multi_system=True)

    def test_initialization(self, fufam):
        """Test the initialization of _FuFamFeatures attributes.

        Args:
            fufam (_FuFamFeatures): Fixture providing a basic _FuFamFeatures object.
        """
        assert fufam.presence == "mandatory"
        assert fufam.duplicate == 2
        assert fufam.multi_system is True
        assert fufam.multi_model is True  # default value


class TestModFuFeatures:
    @pytest.fixture
    def modfu(self):
        """Fixture to initialize a _ModFuFeatures object.

        Returns:
            _ModFuFeatures: A _ModFuFeatures object initialized with test data.
        """
        return _ModFuFeatures(min_mandatory=1, min_total=3, same_strand=True)

    def test_initialization(self, modfu):
        """Test the initialization of _ModFuFeatures attributes.

        Args:
            modfu (_ModFuFeatures): Fixture providing a basic _ModFuFeatures object.
        """
        assert modfu.min_mandatory == 1
        assert modfu.min_total == 3
        assert modfu.same_strand is True
        assert modfu.mandatory == set()  # default value
        assert modfu.accessory == set()  # default value
        assert modfu.forbidden == set()  # default value
        assert modfu.neutral == set()  # default value

    def test_add(self, modfu):
        """Test the add method of _ModFuFeatures.

        Args:
            modfu (_ModFuFeatures): Fixture providing a _ModFuFeatures object.
        """
        mock_child = mock.Mock()
        mock_child.presence = "mandatory"
        modfu.add(mock_child)
        assert mock_child in modfu.mandatory

    def test_check_mandatory_exception(self, modfu):
        """Test the _check method for mandatory items exception.

        Args:
            modfu (_ModFuFeatures): Fixture providing a _ModFuFeatures object.

        Raises:
            Exception: If mandatory items are missing.
        """
        child_accessory = mock.Mock()
        child_accessory.presence = "accessory"
        modfu.add(child_accessory)
        with pytest.raises(Exception, match=f"There are less mandatory {modfu.child_type}"):
            modfu._check()

    def test_check_min_total_exception(self, modfu):
        """Test the _check method for total items exception.

        Args:
            modfu (_ModFuFeatures): Fixture providing a _ModFuFeatures object.

        Raises:
            Exception: If the total number of items is less than the minimum total.
        """
        mock_child = mock.Mock()
        mock_child.presence = "mandatory"
        mock_child.duplicate = 0
        modfu.add(mock_child)
        modfu.min_mandatory = 1
        modfu.min_total = 3
        with pytest.raises(Exception, match=f"There are less {modfu.child_type} than the minimum total"):
            modfu._check()


class TestFamily:
    @pytest.fixture
    def family(self):
        """Fixture to create a Family object.

        Returns:
            Family: A Family object initialized with test data.
        """
        return Family(name="TestFamily", transitivity=1, window=2, presence="mandatory")

    def test_initialization(self, family):
        """Test the initialization of Family attributes.

        Args:
            family (Family): Fixture providing a Family object.
        """
        assert family.name == "TestFamily"
        assert family.transitivity == 1
        assert family.window == 2
        assert family.presence == "mandatory"
        assert family.duplicate == 0
        assert family.multi_system is False
        assert family.multi_model is False
        assert family.exchangeable == set()
        assert family.func_unit is None

    def test_pickling(self, family):
        """Test the serialization and deserialization of a Family object.

        Args:
            family (Family): Fixture providing a Family object.
        """
        serialized = pickle.dumps(family)
        deserialized = pickle.loads(serialized)
        self.test_initialization(deserialized)

    @pytest.fixture
    def func_unit(self):
        """Fixture to create a FuncUnit object.

        Yields:
            FuncUnit: A FuncUnit object initialized with test data.
        """
        yield FuncUnit(name="TestFuncUnit", presence="mandatory", min_mandatory=1, min_total=2, same_strand=True)

    def test_set_func_unit(self, family, func_unit):
        """Test setting the func_unit for a Family object.

        Args:
            family (Family): Fixture providing a Family object.
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        family.func_unit = func_unit
        assert family._parent == func_unit

    @pytest.mark.parametrize("no_fu", [0, "a", 0.5, [], {}, set()])
    def test_set_not_expected_type_raise_type_error(self, family, no_fu):
        """Test that setting a non-FuncUnit raises a TypeError.

        Args:
            family (Family): Fixture providing a Family object.
            no_fu (Any): Non-FuncUnit object to test.

        Raises:
            TypeError: When trying to set a non-FuncUnit object.
        """
        with pytest.raises(TypeError):
            family.func_unit = no_fu

    def test_get_func_unit(self, family, func_unit):
        """Test getting the parent FuncUnit of a Family object.

        Args:
            family (Family): Fixture providing a Family object.
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        family.func_unit = func_unit
        assert family.func_unit == func_unit


class TestFuncUnit:
    @pytest.fixture
    def func_unit(self):
        """Fixture to create a FuncUnit object.

        Returns:
            FuncUnit: A FuncUnit object initialized with test data.
        """
        return FuncUnit(name="TestFuncUnit", transitivity=2, window=3, min_mandatory=1, min_total=3)

    def test_initialization(self, func_unit):
        """Test the initialization of FuncUnit attributes.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        assert func_unit.name == "TestFuncUnit"
        assert func_unit.transitivity == 2
        assert func_unit.window == 3
        assert func_unit.min_mandatory == 1
        assert func_unit.min_total == 3
        assert func_unit.multi_system is False
        assert func_unit.multi_model is False
        assert func_unit.exchangeable == set()
        assert func_unit.mandatory == set()
        assert func_unit.accessory == set()
        assert func_unit.neutral == set()
        assert func_unit.forbidden == set()
        assert func_unit.model is None
        assert func_unit.child_type == 'Family'
        assert func_unit.duplicate == 0
        assert func_unit.same_strand is False

    def test_pickling(self, func_unit):
        """Test the serialization and deserialization of a FuncUnit object.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        serialized = pickle.dumps(func_unit)
        deserialized = pickle.loads(serialized)
        self.test_initialization(deserialized)

    @pytest.fixture
    def model(self) -> Model:
        """
        Generate a Model object.

        Returns:
            Fixture providing a Model object.
        """
        return Model(name="TestModel", transitivity=3, window=2, min_mandatory=1)

    def test_set_model(self, model, func_unit):
        """Test setting the model for a FuncUnit object.

        Args:
            model (Model): Fixture providing a Model object.
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        func_unit.model = model
        assert func_unit._parent == model

    @pytest.mark.parametrize("no_model", [0, "a", 0.5, [], {}, set()])
    def test_set_model_not_isinstance_model_raise_type_error(self, func_unit, no_model):
        """Test that setting a non-Model object raises a TypeError.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
            no_model (Any): Non-Model object to test.

        Raises:
            TypeError: When trying to set a non-Model object.
        """
        with pytest.raises(TypeError):
            func_unit.model = no_model

    def test_get_model(self, func_unit, model):
        """Test getting the parent model of a FuncUnit object.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
            model (Model): Fixture providing a Model object.
        """
        func_unit.model = model
        assert func_unit.model == model

    def test_check_func_unit(self, func_unit):
        """Test the check_func_unit method of FuncUnit.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.

        Raises:
            Exception: If any validation fails during the check.
        """
        with pytest.raises(Exception):
            func_unit.check_func_unit()

    @pytest.fixture
    def families(self) -> Set[Family]:
        """
        Generate a Set of Family objects.

        Returns:
            A random set of Family objects.
        """
        families = set()
        for i in range(randint(2, 10)):
            t = randint(1, 5)
            family = Family(name=f"TestFamily{i}", transitivity=t, window=t + 1,
                            presence="mandatory")
            families.add(family)
        return families

    def test_add_families(self, func_unit, families):
        """Test adding families to a FuncUnit object.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
            families (set): Fixture providing a set of Family objects.
        """
        for family in families:
            func_unit.add(family)
        assert set(func_unit.families) == families

    def test_equal_hash(self, func_unit):
        """Test equality of hash of FuncUnit objects.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        fu2 = FuncUnit(name="TestFuncUnit", transitivity=2, window=3, min_mandatory=1, min_total=3)
        assert hash(fu2) == hash(func_unit)

    def test_func_unit_equality(self, families, func_unit):
        """Test equality of FuncUnit objects with families added.

        Args:
            families (set): Fixture providing a set of Family objects.
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
        """
        fu2 = FuncUnit(name="TestFuncUnit", transitivity=2, window=3, min_mandatory=1, min_total=3)
        for fam in families:
            func_unit.add(fam)
            fu2.add(fam)
        assert fu2 == func_unit

    def test_pickling_with_families(self, func_unit, families):
        """Test serialization and deserialization of a FuncUnit with families.

        Args:
            func_unit (FuncUnit): Fixture providing a FuncUnit object.
            families (set): Fixture providing a set of Family objects.
        """
        for family in families:
            func_unit.add(family)
        serialized = pickle.dumps(func_unit)
        deserialized = pickle.loads(serialized)
        assert set(deserialized.gene_families) == families


class TestModel:
    @pytest.fixture
    def model(self):
        """Fixture to create a Model object.

        Returns:
            Model: A Model object initialized with test data.
        """
        return Model(name="TestModel", transitivity=3, window=4, min_mandatory=1, min_total=2, canonical=["Canon1", "Canon2"])

    def test_initialization(self, model):
        """Test the initialization of Model attributes.

        Args:
            model (Model): Fixture providing a Model object.
        """
        assert model.name == "TestModel"
        assert model.transitivity == 3
        assert model.window == 4
        assert model.min_mandatory == 1
        assert model.min_total == 2
        assert model.mandatory == set()
        assert model.accessory == set()
        assert model.neutral == set()
        assert model.forbidden == set()
        assert model.child_type == 'FuncUnit'
        assert model.same_strand is False
        assert model.canonical == ["Canon1", "Canon2"]

    def test_pickling(self, model):
        """Test the serialization and deserialization of a Model object.

        Args:
            model (Model): Fixture providing a Model object.
        """
        serialized = pickle.dumps(model)
        deserialized = pickle.loads(serialized)
        self.test_initialization(deserialized)

    @pytest.fixture
    def families(self) -> List[Family]:
        """
        Generate a list of Family objects.

        Returns:
            List[Family]: Fixture providing a set of Family objects.
        """
        families = []
        for i in range(randint(2, 10)):
            t = randint(1, 5)
            family = Family(name=f"TestFamily{i}", transitivity=t, window=t + 1,
                            presence="mandatory")
            families.append(family)
        return families

    @pytest.fixture
    def func_units(self, families) -> Set[FuncUnit]:
        """
        Generate a set of FuncUnit objects.

        Args:
            families: Fixture providing a set of Family objects.

        Returns:
            Set[FuncUnit]: Fixture providing a set of FuncUnit objects.
        """
        func_units = set()
        for i in range(1, len(families) + 1):
            t = randint(1, 5)
            fu = FuncUnit(name=f"TestFu{i}", transitivity=t, window=t + 1,
                          presence="mandatory")
            fu.add(families[i - 1])
            func_units.add(fu)
        return func_units

    def test_add_func_units(self, model, func_units):
        """Test adding func_units to a Model object.

        Args:
            model (Model): Fixture providing a Model object.
            func_units (set): Fixture providing a set of FuncUnit objects.
        """
        for fu in func_units:
            model.add(fu)
        assert set(model.func_units) == func_units

    def test_pickling_with_func_units(self, model, func_units):
        """Test serialization and deserialization of a Model with func_units.

        Args:
            model (Model): Fixture providing a Model object.
            func_units (set): Fixture providing a set of FuncUnit objects.
        """
        for fu in func_units:
            model.add(fu)
        serialized = pickle.dumps(model)
        deserialized = pickle.loads(serialized)
        assert set(deserialized.func_units) == func_units
