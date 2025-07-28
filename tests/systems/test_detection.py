from collections import defaultdict

import networkx as nx

from panorama.geneFamily import GeneFamily
from panorama.systems.system import SystemUnit
from panorama.systems.models import Models
from panorama.systems.detection import (
    get_functional_unit_gene_families,
    search_unit_in_cc,
    search_unit_in_context,
    search_unit_in_combination,
    search_system_units,
    search_system,
    search_systems,
    search_for_system,
    search_system_units,
    check_for_needed_units,
    get_system_unit_combinations,
)


def test_get_functional_unit_gene_families_simple_case(multi_unit_model, simple_gf2fam):
    """Tests that the function correctly retrieves the (mandatory | accessory) and neutral families within a functional unit."""
    # Get the (mandatory | accessory) and neutral family sets per functional unit with at least one occurence in the pangenome
    gfs = list(simple_gf2fam.keys())
    fu1 = multi_unit_model.get("fu1")  # fu1 = {GF0, GF1, GF2, GF3, GF4, GF5} | {GF9}
    assert get_functional_unit_gene_families(fu1, gfs, simple_gf2fam) == (
        set(gfs[:6]),
        {gfs[-1]},
    )


def test_get_functional_unit_gene_families_exclude_non_pangenome_fams(
    multi_unit_model, simple_gf2fam
):
    """Tests that the function excludes families in the functional unit with no association to any GF in the pangenome."""
    gfs = list(simple_gf2fam.keys())
    fu2 = multi_unit_model.get("fu2")  # fu2 = {GF6, GF7, GF8, GF10} | {GF11}
    assert get_functional_unit_gene_families(fu2, gfs, simple_gf2fam) == (
        set(gfs[6:9]),
        set(),
    )  # GF10 and 11 excluded as they are not in the pangenome (not in gf2fam)


def test_search_unit_in_cc(
    multi_unit_model, simple_gf2fam, simple_fam2source, simple_pangenome
):
    """Tests that the function correctly identifies units in a connected component of the filtered context graph."""
    # simple_pangenome argument needed to assign metadata to GFs
    detected_units = set()
    gfs = list(simple_gf2fam.keys())
    comb_families = set(
        gfs[:6]
    )  # a combination that is subset of (mandatory | accessory) families of fu1 AND found in a connected component of the context graph
    fu = multi_unit_model.get("fu1")

    filtered_cc = nx.Graph()  # a filtered connected component of the context graph
    fu_families_list = sorted(
        comb_families, key=lambda gf: gf.ID
    )  # sort to ensure consistent order
    for i in range(len(fu_families_list) - 1):
        if i == 3:  # skip to create a disconnected component
            continue
        filtered_cc.add_edge(fu_families_list[i], fu_families_list[i + 1])
    context_gf = GeneFamily(
        family_id=11, name="GF11"
    )  # add context family node (not in the model)
    filtered_cc.add_edge(fu_families_list[0], context_gf)
    # filtered_cc = GF0 -- GF1 -- GF2 -- GF3 -- GF11
    #               GF4 -- GF5

    # returns detected model families combination as frozenset
    assert search_unit_in_cc(
        filtered_cc,
        comb_families,
        fu,
        "source1",
        simple_gf2fam,
        simple_fam2source,
        detected_units,
    ) == {frozenset(gfs[:4])}

    # the argument detected_units is updated with the complete identified unit, including context family
    unit_families = set(next(iter(detected_units)).families)
    assert unit_families == set(gfs[:4]) | {context_gf}


def test_search_unit_in_combination(
    multi_unit_model, simple_gf2fam, simple_fam2source, simple_matrix, simple_pangenome
):
    """Tests that the function correctly identifies units in a connected component of the context graph starting a set of potential combinations."""
    fu = multi_unit_model.get("fu1")
    gfs = list(simple_gf2fam.keys())
    families_in_cc = set(
        gfs[:6]
    )  # (mandatory | accessory) families of fu1 found in a connected component of the context graph

    cc = nx.Graph()  # a connected component of the context graph
    for i, _ in enumerate(gfs[:6]):
        cc.add_edge(gfs[i], gfs[i + 1])
    # cc = GF0 -- GF1 -- GF2 -- GF3 -- GF4 -- GF5
    combinations_in_cc = list(
        {frozenset(gfs[:4])}
    )  # list of all combinations subset of families_in_cc
    combinations2orgs = defaultdict(
        set,
        {
            frozenset(
                gfs[:4]
            ): {  # dict of all comb-to-orgs associations returned by compute_gene_context_graph of ppanggolin
                "org1",
                "org2",
            }
        },
    )  # used for local filtering
    detected = search_unit_in_combination(
        cc,
        families_in_cc,
        simple_gf2fam,
        simple_fam2source,
        fu,
        "source1",
        simple_matrix,
        combinations_in_cc.copy(),
        combinations2orgs,
    )  # a copy of combinations_in_cc is passed since it is modified in-place
    assert set(next(iter(detected)).families) == set(gfs[:4])

    # Add context family node
    context_gf = GeneFamily(family_id=11, name="GF11")
    cc.add_edge(gfs[0], context_gf)
    detected = search_unit_in_combination(
        cc,
        families_in_cc,
        simple_gf2fam,
        simple_fam2source,
        fu,
        "source1",
        simple_matrix,
        combinations_in_cc.copy(),
        combinations2orgs,
    )
    assert set(next(iter(detected)).families) == set(gfs[:4]) | {context_gf}


def test_search_unit_in_context(
    multi_unit_model, simple_gf2fam, simple_fam2source, simple_matrix, simple_pangenome
):
    """Tests that the function correctly identifies units given the full context graph"""
    detected = set()
    gfs = list(simple_gf2fam.keys())
    fu_families = set(
        gfs[:6]
    )  # (mandatory | accessory) families of fu1 that are present in the pangenome

    combinations2orgs = defaultdict(
        set,
        {
            frozenset(
                gfs[:4]
            ): {  # dict of all comb-to-orgs associations returned by compute_gene_context_graph of ppanggolin
                "org1",
                "org2",
            }
        },
    )  # used for local filtering
    func_unit = multi_unit_model.get("fu1")
    context_graph = (
        nx.Graph()
    )  # context graph returned by compute_gene_context_graph of ppanggolin
    for i, _ in enumerate(gfs[:-1]):
        if i == 5:  # skip to create a disconnected component
            continue
        context_graph.add_edge(gfs[i], gfs[i + 1])  # genomes = {'org1'}
    # context_graph = GF0 -- GF1 -- GF2 -- GF3 -- GF4 -- GF5
    #                 GF6 -- GF7 -- GF8
    detected = search_unit_in_context(
        context_graph,
        fu_families,
        simple_gf2fam,
        simple_fam2source,
        simple_matrix,
        func_unit,
        "source1",
        combinations2orgs,
        local=True,
    )
    assert set(next(iter(detected)).families) == set(
        gfs[:4]
    )  # detected unit corresponding to the combination in combinations2orgs


def test_search_system_units(
    multi_unit_model, simple_gf2fam, simple_fam2source, simple_pangenome
):
    """Tests that the context graph is properly constructed, the potential combinations are identified, and the units are detected."""
    gfs = list(simple_gf2fam.keys())
    model_families = set(gfs[:6])  # all GFs corresponding to model families
    detected_systems = search_system_units(
        multi_unit_model,
        model_families,
        simple_gf2fam,
        simple_fam2source,
        source="source1",
    )
    # This method executes ppanggolin context graph construction -> context_graph = GF0 -- GF1 -- GF2 -- GF3 -- GF4 -- GF5 -- GF6 -- GF7 -- GF8
    # since all GFs correspond to sequential genes of the same contig of the same organism as defined in simple_gfs fixture

    assert detected_systems.keys() == {"fu1"}  # only fu1 should be detected
    assert set(next(iter(detected_systems["fu1"])).families) == set(
        gfs[:7]
    )  # GF6 detected as context family (since transitivity = window = 1)
    assert set(next(iter(detected_systems["fu1"])).models_families) == set(
        gfs[:6]
    )  # model families of fu1

    # Increment window size of fu1
    multi_unit_model.get("fu1").window = 2
    detected_systems = search_system_units(
        multi_unit_model,
        model_families,
        simple_gf2fam,
        simple_fam2source,
        source="source1",
    )
    assert set(next(iter(detected_systems["fu1"])).families) == set(
        gfs[:8]
    )  # GF7 is now detected as context family
    assert set(next(iter(detected_systems["fu1"])).models_families) == set(
        gfs[:6]
    )  # model families of fu1 remain the same


def test_check_for_needed_units_satisfied(simple_gfs, multi_unit_model):
    """Tests that the function correctly identifies when detected functional units satisfy the model requirements."""
    gf2meta_info = {
        gf: ("source1", 1) for gf in simple_gfs[:6]
    }  # dict of GFs-to-metadata for the GFs of one detected unit
    su = SystemUnit(
        functional_unit=multi_unit_model.get("fu1"),
        source="source1",  # a detected system unit
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo=gf2meta_info,
    )
    detected_units = {
        "fu1": su
    }  # dict of all detected units corresponding to each functional unit in the model

    # multi_unit_model has fu1 mandatory and fu2 accessory => satisfied
    assert check_for_needed_units(detected_units, multi_unit_model) is True


def test_check_for_needed_units_not_satisfied(simple_gfs, multi_unit_model):
    """Tests that the function correctly identifies when detected functional units do not satisfy the model requirements."""
    gf2meta_info = {
        gf: ("source1", 1) for gf in simple_gfs[:6]
    }  # dict of GFs-to-metadata for the GFs of one detected unit
    su = SystemUnit(
        functional_unit=multi_unit_model.get("fu1"),
        source="source1",  # a detected system unit
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo=gf2meta_info,
    )
    detected_units = {
        "fu1": su
    }  # dict of all detected units corresponding to each functional unit in the model

    # Add fu2 as mandatory functional unit -> model no longer satisfied
    multi_unit_model.mandatory.add(multi_unit_model.get("fu2"))
    multi_unit_model.min_mandatory = 2
    assert check_for_needed_units(detected_units, multi_unit_model) is False


def test_get_system_unit_combinations(simple_gfs, multi_unit_model):
    """Tests that the function correctly retrieves all unit combinations satisfying the model requirements."""
    gf2meta_info = {
        gf: ("source1", 1) for gf in simple_gfs[:6]
    }  # dict of GFs-to-metadata for the GFs of one detected unit
    su1 = SystemUnit(
        functional_unit=multi_unit_model.get("fu1"),
        source="source1",
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo=gf2meta_info,
    )
    detected_units = {
        "fu1": {su1}
    }  # dict of all detected units corresponding to each functional unit in the model

    # Return all system unit combinations satisifying the model requirements
    assert get_system_unit_combinations(detected_units, multi_unit_model) == [[su1]]

    # Let fu2 be detected
    gfs = simple_gfs[6:9] + [
        GeneFamily(family_id=10, name="GF10")
    ]  # GF10 is the mandatory family of fu2
    gf2meta_info = {gf: ("source1", 1) for gf in gfs}
    su2 = SystemUnit(
        functional_unit=multi_unit_model.get("fu2"),
        source="source1",
        gene_families=set(simple_gfs),
        families_to_metainfo=gf2meta_info,
    )
    detected_units = {"fu1": {su1}, "fu2": {su2}}

    # fu2 is accessory -> model satisfied with and without su2
    assert [
        set(c) for c in get_system_unit_combinations(detected_units, multi_unit_model)
    ] == [
        {su1},
        {su1, su2},
    ]  # convert combs to sets to avoid order discrepancies

    # Let fu2 be mandatory => model only satisfied with both su and su2
    multi_unit_model.mandatory.add(multi_unit_model.get("fu2"))
    multi_unit_model.accessory.remove(multi_unit_model.get("fu2"))
    multi_unit_model.min_mandatory = 2
    assert [
        set(c) for c in get_system_unit_combinations(detected_units, multi_unit_model)
    ] == [{su1, su2}]


def test_search_system(
    multi_unit_model, simple_gf2fam, simple_fam2source, simple_pangenome
):
    """Tests that the function correctly identifies systems in the pangnomes after units detection."""
    gfs = list(simple_gf2fam.keys())
    model_fams = set(gfs[:6])  # all GFs corresponding to model families (fu1 only)

    detected_system = next(
        iter(
            search_system(
                multi_unit_model,
                model_fams,
                simple_gf2fam,
                simple_fam2source,
                "source1",
            )
        )
    )  # extract the first detected system
    expected_unit = SystemUnit(
        functional_unit=multi_unit_model.get("fu1"),
        source="source1",
        gene_families=set(
            gfs[:7]
        ),  # GF6 does not have corresponding metainfo => context family
        families_to_metainfo={gf: ("source1", 1) for gf in gfs[:6]},
    )

    # the first detected unit of the detected system corresponds to the expected unit (__eq__ based on models_families attribute only)
    assert next(iter((set(detected_system.units)))) == expected_unit
    assert set(expected_unit.models_families) == set(
        gfs[:6]
    )  # ensure models_families is as expected
    assert (
        set(next(iter((set(detected_system.units)))).families)
        == set(expected_unit.families)
        == set(gfs[:7])
    )  # unit families include GF6 as well


def test_search_systems(simple_gfs, multi_unit_model, simple_pangenome):
    """Tests that the function correctly identifies systems in the pangenome across different models."""
    # Executes all the previous methods, and all the methods in utils, given the model and pangenome as arguments
    models = Models(models={multi_unit_model.name: multi_unit_model})

    # No return value; systems added to the pangenome
    assert (
        search_systems(
            models,
            simple_pangenome,
            "source1",
            metadata_sources=["source1"],
            sensitivity=3,
        )
        is None
    )
    expected_units = {
        SystemUnit(
            functional_unit=multi_unit_model.get("fu1"),
            source="source1",
            gene_families=set(simple_gfs[:7]),
            families_to_metainfo={gf: ("source1", 1) for gf in simple_gfs[:6]},
        )
    }
    assert set(next(simple_pangenome.systems).units) == expected_units


def test_search_for_system():  # broken logic; to be fixed
    pass
