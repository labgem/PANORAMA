#! /usr/bin/env python3

import pytest
from random import randint
from typing import Set
import pickle

from panorama.geneFamily import GeneFamily
from panorama.systems.system import SystemUnit, System
from panorama.systems.models import Model, FuncUnit, Family


def gen_families(min_mandatory: int, min_accessory, nb_forbidden: int, nb_neutral: int) -> Set[Family]:
    families = set()

    def gen_fam(nb, presence, **kwargs):
        for i in range(nb):
            families.add(Family(name=f"Family_{presence}_{i}", presence=presence, **kwargs))

    gen_fam(min_mandatory, "mandatory")
    gen_fam(min_accessory, "accessory")
    gen_fam(nb_forbidden, "forbidden")
    gen_fam(nb_neutral, "mandatory")

    return families


@pytest.fixture
def func_unit():
    min_mandatory = randint(1, 3)
    min_accessory = randint(0, 2)
    nb_forbidden = randint(0, 1)
    nb_neutral = randint(0, 1)
    t = randint(1, 5)
    fu = FuncUnit(name="TestFuncUnit", presence="mandatory", min_mandatory=min_mandatory,
                  min_total=min_mandatory + min_accessory, transitivity=t, window=t + 1)
    for fam in gen_families(min_mandatory, min_accessory, nb_forbidden, nb_neutral):
        fu.add(fam)

    return fu


@pytest.fixture
def model(func_unit):
    return Model(name='TestModel', mandatory={func_unit})


@pytest.fixture
def gene_families_to_families(func_unit: FuncUnit):
    families = list(func_unit.families)
    gf2fam = {}
    for i in range(len(families)):
        gf = GeneFamily(family_id=0, name=f"GF_{i}")
        gf2fam[gf] = ("test", randint(1, 5))
    return gf2fam


@pytest.fixture
def gene_families(gene_families_to_families):
    gfs = set(gene_families_to_families.keys())
    for i in range(len(gfs), len(gfs) + randint(2, 5)):
        gfs.add(GeneFamily(family_id=0, name=f"GF_{i}"))
    return gfs


class TestSystemUnit:
    @pytest.fixture
    def unit(self, func_unit, gene_families, gene_families_to_families):
        return SystemUnit(functional_unit=func_unit, source="test", gene_families=gene_families,
                          gene_families_to_metainfo=gene_families_to_families)

    def test_init(self, unit, func_unit, gene_families, gene_families_to_families):
        assert unit.ID == SystemUnit._id_counter == 1
        assert unit.functional_unit == func_unit
        assert unit.source == "test"
        assert unit._families_getter == {gf.name: gf for gf in gene_families}
        assert unit._families2metainfo == {gf: gene_families_to_families.get(gf, ('', 0)) for gf in gene_families}
        assert unit._regions_getter == {}
        assert unit._spots_getter == {}
        assert unit._modules_getter == {}
        assert unit._models_families is None
        assert unit._system is None

    def test_pickable(self, unit, func_unit, gene_families, gene_families_to_families):
        serialized = pickle.dumps(unit)
        deserialized = pickle.loads(serialized)
        assert unit.ID == deserialized.ID
        assert unit.functional_unit == deserialized.functional_unit
        assert unit.source == deserialized.source
        assert unit._families_getter == deserialized._families_getter
        assert unit._families2metainfo == deserialized._families2metainfo
        assert unit._regions_getter == deserialized._regions_getter
        assert unit._spots_getter == deserialized._spots_getter
        assert unit._modules_getter == deserialized._modules_getter
        assert unit._models_families == deserialized._models_families
        assert unit._system == deserialized._system


class TestSystem:
    @pytest.fixture
    def unit(self, func_unit, gene_families, gene_families_to_families):
        return SystemUnit(functional_unit=func_unit, source="test", gene_families=gene_families,
                          gene_families_to_metainfo=gene_families_to_families)

    @pytest.fixture
    def system(self, unit, model):
        return System(model=model, source="test", units={unit})

    def test_init(self, system, unit, model):
        assert system.ID == str(System._id_counter) == str(1)
        assert system.model == model
        assert system.source == "test"
        assert system._unit_getter == {unit.name: unit}
        assert system.canonical == set()

    def test_pickable(self, system, unit, model):
        serialized = pickle.dumps(system)
        deserialized = pickle.loads(serialized)
        self.test_init(deserialized, unit, model)
