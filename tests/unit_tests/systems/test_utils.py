from collections import defaultdict

import networkx as nx
import pandas as pd
from natsort import natsorted
from ppanggolin.meta.meta import assign_metadata

from panorama.geneFamily import GeneFamily
from panorama.systems.models import Family
from panorama.systems.utils import (
    check_for_families,
    check_needed_families,
    dict_families_context,
    filter_global_context,
    filter_local_context,
    get_gfs_matrix_combination,
    get_metadata_to_families,
)
from tests.conftest import DummyGeneFamily


def test_get_metadata_to_families_simple_case(simple_pangenome):
    """Tests that metadata for each source is properly mapped to gene families."""
    annot2fam = get_metadata_to_families(simple_pangenome, sources=["source1"])
    gfs = [gf for gf in simple_pangenome.gene_families]

    assert annot2fam == {
        "source1": defaultdict(set, {f"protein{i}": {gfs[i]} for i in range(len(gfs))})
    }  # mapping from protein annotations to GFs


def test_get_metadata_to_families_multiple_annotations_case(simple_pangenome):
    """Tests that many-to-one associations in both directions are properly handled."""
    extra_gfs = [GeneFamily(family_id=i, name=f"GF{i}") for i in range(10, 18)]
    gfs = list(simple_pangenome.gene_families) + extra_gfs
    # Add multiple GFs for some protein annotations: GF_i -> protein_i
    #                                                GF_i+10 -> protein_i
    families = [gf.name for gf in gfs[:10]] + [gf.name for gf in gfs[10:18]]
    protein_names = [f"protein{i}" for i in range(10)] + [
        f"protein{i}" for i in range(8)
    ]
    # Add mutiple annotations for some GFs: GF0 -> protein0, protein8
    #                                       GF1 -> protein1, protein9
    families += [gfs[0].name, gfs[1].name]
    protein_names += ["protein8", "protein9"]

    for gf in extra_gfs:
        simple_pangenome.add_gene_family(gf)
    metadata = pd.DataFrame(
        {
            "families": families,
            "protein_name": protein_names,
        }
    )
    # assign updated metadata as source2
    assign_metadata(metadata, simple_pangenome, source="source2", metatype="families")

    assert get_metadata_to_families(simple_pangenome, sources=["source2"]) == {
        "source2": defaultdict(
            set,
            {
                f"protein{i}": {gfs[i], gfs[i + 10]}  # protein_i -> {GF_i, GF_i+10}
                for i in range(8)
            }
            | {
                f"protein{i}": {gfs[i], gfs[i - 8]}  # protein8 -> {GF0, GF8}
                for i in range(8, 10)
            },
        )
    }  # protein9 -> {GF1, GF9}


def test_dict_families_context_simple_case(simple_gfs, single_unit_model):
    """
    Tests that the function returns the correct set of model GFs,
    mapping to protein annotations, and source per annotation.
    """
    annot2fam = {
        "source1": defaultdict(
            set, {f"protein{i}": {simple_gfs[i]} for i in range(len(simple_gfs))}
        )
    }  # protein_i -> {GF_i}
    gf2fam, fam2source = dict_families_context(single_unit_model, annot2fam)
    # mandatory: protein0, protein1, protein2
    # accessory: protein3, protein4, protein5
    assert gf2fam == defaultdict(
        set,
        {
            simple_gfs[i]: {
                [f for f in single_unit_model.families if f.name == f"protein{i}"][0]
            }
            for i in range(6)
        },
    )  # mapping from GFs to model families
    assert fam2source == {
        f"protein{i}": "source1" for i in range(6)
    }  # annotation source mapping


def test_dict_families_context_many_to_many_case(simple_gfs, single_unit_model):
    """Tests that the function properly handles many-to-many associations between GFs and protein annotations."""
    # Add many-to-many GF-protein associations
    simple_gfs.extend([GeneFamily(family_id=i, name=f"GF{i}") for i in range(10, 18)])
    annot2fam = {
        "source2": defaultdict(
            set,
            {
                f"protein{i}": {
                    simple_gfs[i],
                    simple_gfs[i + 10],
                }  # protein_i -> {GF_i, GF_i+10}
                for i in range(8)
            }
            | {f"protein{i}": {simple_gfs[i], simple_gfs[i - 8]} for i in range(8, 10)},
        )
    }  # protein8 -> {GF0, GF8}; protein9 -> {GF1, GF9}
    # Add new family to model
    new_family = Family(name="protein8", presence="mandatory")
    fu = next(single_unit_model.func_units)
    fu.mandatory.add(new_family)

    gf2fam, fam2source = dict_families_context(single_unit_model, annot2fam)

    family_lookup = {f.name: f for f in single_unit_model.families}
    assert gf2fam == defaultdict(
        set, {simple_gfs[i]: {family_lookup[f"protein{i}"]} for i in range(1, 6)}
    ) | defaultdict(
        set, {simple_gfs[i + 10]: {family_lookup[f"protein{i}"]} for i in range(6)}
    ) | {simple_gfs[0]: {family_lookup["protein0"], new_family}} | {
        simple_gfs[8]: {new_family}
    }
    assert fam2source == {f"protein{i}": "source2" for i in range(6)} | {
        "protein8": "source2"
    }


def test_get_gfs_matrix_combination_simple_case(simple_gfs):
    """Tests that the function correctly returns the binary matrix of GF-to-protein annotation associations."""
    gf2fam = defaultdict(
        set,
        {
            simple_gfs[i]: {
                Family(
                    name=f"protein{i}", presence="mandatory" if i < 3 else "accessory"
                )
            }
            for i in range(6)
        },
    )
    fu_families = set(simple_gfs[:6])
    matrix = get_gfs_matrix_combination(fu_families, gf2fam)
    # binary matrix with protein annotations as rows and GFs as columns
    assert (
        matrix.sort_index(axis=0)
        .sort_index(axis=1)
        .equals(  # sort rows and columns to ensure consistency
            pd.DataFrame(
                [[int(i == j) for j in range(6)] for i in range(6)],
                index=[f"protein{i}" for i in range(6)],
                columns=[f"GF{i}" for i in range(6)],
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
    )  # Ensure the matrix matches the expected simple diagonal matrix


def test_get_gfs_matrix_combination_many_to_many_case(simple_gfs):
    """Tests that the function correctly handles many-to-many associations between GFs and protein annotations."""
    # Consider many-to-many GF-protein associations
    simple_gfs.extend([GeneFamily(family_id=i, name=f"GF{i}") for i in range(10, 18)])
    fu_families = set(simple_gfs[:6]) | set(simple_gfs[10:16]) | {simple_gfs[8]}
    gf2fam = (
        defaultdict(
            set,
            {
                simple_gfs[i]: {
                    Family(
                        name=f"protein{i}",
                        presence="mandatory" if i < 3 else "accessory",
                    )
                }
                for i in range(1, 6)
            },
        )
        | defaultdict(
            set,
            {
                simple_gfs[i + 10]: {
                    Family(
                        name=f"protein{i}",
                        presence="accessory" if i < 3 else "accessory",
                    )
                }
                for i in range(6)
            },
        )
        | {
            simple_gfs[0]: {
                Family(name="protein0", presence="mandatory"),
                Family(name="protein8", presence="mandatory"),
            }
        }
        | {simple_gfs[8]: {Family(name="protein8", presence="mandatory")}}
    )
    matrix = get_gfs_matrix_combination(fu_families, gf2fam)
    matrix = matrix.reindex(
        index=sorted(matrix.index), columns=natsorted(matrix.columns)
    )  # natural sort for columns to ensure consistent order

    annotations = sorted(
        {family.name for gf in fu_families for family in gf2fam[gf]}
    )  # rows
    model_fams = natsorted(gf.name for gf in fu_families)  # columns
    expected_matrix = pd.DataFrame(
        [
            [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],  # protein0: GF0, GF10
            [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],  # protein1: GF1, GF11
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],  # protein2: GF2, GF12
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],  # protein3: GF3, GF13
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],  # protein4: GF4, GF14
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],  # protein5: GF5, GF15
            [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],  # protein8: GF0, GF8
        ],
        index=annotations,
        columns=model_fams,
    )

    assert matrix.equals(
        expected_matrix
    )  # Ensure the matrix matches the expected output


def test_check_needed_families(single_unit_model, simple_matrix):
    """Tests that the function correctly checks if the model is satisfied given the matrix of GFs to annotations matrix
    """
    fu = next(single_unit_model.func_units)  # the first functional unit of the model
    # min_total = 4, min_mandatory = 2
    # 3 mandatory and 3 accessory families covered => model satisfied
    assert check_needed_families(simple_matrix, fu) is True

    # Updating min_total to 7 -> model unsatisfied
    fu.min_total = 7
    assert (
        check_needed_families(simple_matrix, next(single_unit_model.func_units))
        is False
    )


def test_check_for_families(
    simple_gf2fam, simple_fam2source, single_unit_model, simple_pangenome
):
    """
    Tests that the function correctly checks if the model is satisfied given a set of GFs and the metadata dictionary
    AND returns the correct GFs-to-source mapping dictionary.
    """
    # simple_pangenome is needed to assign metadata to families
    gfs_in_cc = list(simple_gf2fam.keys())[
        :4
    ]  # assume first 4 GFs in a connected component
    fu = next(single_unit_model.func_units)
    # 3 mandatory and 1 accessory families covered => model satisfied
    check, gf2meta_info = check_for_families(
        gfs_in_cc, simple_gf2fam, simple_fam2source, fu
    )
    assert check
    assert gf2meta_info == {
        gf: ("source1", 1) for gf in gfs_in_cc
    }  # mapping from GFs to (source, best metadata ID)

    # Incrementing min_total to 5
    fu.min_total = 5
    assert check_for_families(gfs_in_cc, simple_gf2fam, simple_fam2source, fu) == (
        False,
        {},
    )  # model unsatisfied


def test_check_for_families_forbidden_case(
    simple_gf2fam, simple_fam2source, single_unit_model, simple_pangenome
):
    """Tests that the function properly handles forbidden families in the model."""
    gfs_in_cc = list(simple_gf2fam.keys())[
        :7
    ]  # assume first 7 GFs in a connected component

    forbidden_fam = Family(
        name="protein6", presence="forbidden"
    )  # let GF6 be forbidden
    modified_gf2fam = simple_gf2fam.copy()
    modified_gf2fam[gfs_in_cc[6]].add(
        forbidden_fam
    )  # add forbidden annotation to GF6; GF6 -> 1 accessory family + 1 forbidden family

    fu = next(single_unit_model.func_units)
    fu.forbidden = {forbidden_fam}

    assert check_for_families(gfs_in_cc, simple_gf2fam, simple_fam2source, fu) == (
        False,
        {},
    )  # model unsatisfied due to forbidden family


def test_filter_global_context():
    """Tests the global context filtering function."""
    gf1 = DummyGeneFamily("GF1", ["org1", "org2"])
    gf2 = DummyGeneFamily("GF2", ["org1", "org2"])
    gf3 = DummyGeneFamily("GF3", ["org1", "org2", "org3"])

    G = nx.Graph()
    G.add_edge(gf1, gf2, genomes={"org1", "org2"})  # proportion = 1.0 for both
    G.add_edge(gf1, gf3, genomes={"org1"})  # proportion = 0.5 for gf1, 0.333 for gf3

    filtered = filter_global_context(G, jaccard_threshold=0.8)

    # Only edge (gf1, gf2) remains
    assert set(filtered.edges()) == {(gf1, gf2)}
    edge_data = filtered.get_edge_data(gf1, gf2)
    assert edge_data["f1"] == "GF1"
    assert edge_data["f2"] == "GF2"
    assert edge_data["f1_jaccard_gene"] == edge_data["f2_jaccard_gene"] == 1.0


def test_filter_local_context():
    """Tests the local context filtering function."""
    gf1 = DummyGeneFamily("GF1", ["org1", "org2"])
    gf2 = DummyGeneFamily("GF2", ["org1", "org2"])
    gf3 = DummyGeneFamily("GF3", ["org1", "org2", "org3"])

    orgs = {"org1"}

    G = nx.Graph()
    G.add_edge(
        gf1, gf2, genomes={"org1", "org2"}
    )  # local proportion = 1.0 for both (since only org1 considered in orgs)
    G.add_edge(gf1, gf3, genomes={"org1"})  # local proportion = 1.0 for both

    filtered = filter_local_context(G, orgs, jaccard_threshold=0.8)

    # both edges remain
    assert set(filtered.edges()) == {(gf1, gf2), (gf1, gf3)}
    edge_data1 = filtered.get_edge_data(gf1, gf2)
    edge_data2 = filtered.get_edge_data(gf1, gf3)
    assert edge_data1["f1_jaccard_gene"] == edge_data1["f2_jaccard_gene"] == 1.0
    assert edge_data2["f1_jaccard_gene"] == edge_data2["f2_jaccard_gene"] == 1.0

    orgs = {
        "org2",
        "org3",
    }  # Assume combination in org2 and org3 instead -> local proportion = 0 for both sides of edge (gf1, gf3)

    G.add_edge(
        gf2, gf3, genomes={"org2"}
    )  # local proportion = 1.0 for gf2, 0.5 for gf3
    filtered = filter_local_context(G, orgs, jaccard_threshold=0.65)
    assert set(filtered.edges()) == {(gf1, gf2)}  # only edge (gf1, gf2) remains
