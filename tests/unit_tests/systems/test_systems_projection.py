from collections import namedtuple

import networkx as nx
import pandas as pd
import pandas.testing as pdt

from panorama.systems.system import System, SystemUnit
from panorama.systems.systems_projection import (
    compute_gene_components,
    compute_genes_graph,
    get_org_df,
    has_short_path,
    project_pangenome_systems,
    project_unit_on_organisms,
    system_projection,
    unit_projection,
    write_projection_systems,
)


def test_has_short_path():
    """
    Tests if the function correctly identifies whether a path of length <= n exists between
    a set of given nodes in a graph.
    """
    G = nx.Graph()
    G.add_edges_from([(1, 2), (2, 3), (3, 4)])  # simple graph: 1-2-3-4
    assert has_short_path(G, [1, 2, 3, 4], 1)
    assert not has_short_path(G, [1, 2, 4, 3], 1)  # function is order-dependent


def test_compute_genes_graph(simple_gfs, simple_fu, simple_orgs, simple_pangenome):
    """
    Tests that the function correctly computes the genes graph given a set of genes corresponding
    to model GFs in one organism.
    """
    gf2meta_info = {
        gf: ("source1", 1) for gf in simple_gfs[:4]
    }  # dict of GFs-to-metadata for the GFs of one detected unit
    unit = SystemUnit(
        functional_unit=simple_fu,
        source="source1",
        gene_families=set(simple_gfs[:4]),
        families_to_metainfo=gf2meta_info,
    )
    org_unit_genes = {
        gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:4])
    }  # genes of the unit model families present in the organism

    genes_graph = compute_genes_graph(org_unit_genes, unit)

    assert org_unit_genes == {
        gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:4])
    }  # ensure correct genes are included
    assert set(genes_graph.nodes) == {gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:6])}

    G = nx.Graph()
    expected_edges = [
        (
            simple_gfs[0]["gene_0_0"],
            simple_gfs[1]["gene_1_0"],
            {"transitivity": 0},
        ),  # 0 -> 1
        (
            simple_gfs[0]["gene_0_0"],
            simple_gfs[2]["gene_2_0"],
            {"transitivity": 1},
        ),  # 0 -> 2
        (
            simple_gfs[1]["gene_1_0"],
            simple_gfs[2]["gene_2_0"],
            {"transitivity": 0},
        ),  # 1 -> 2
        (
            simple_gfs[1]["gene_1_0"],
            simple_gfs[3]["gene_3_0"],
            {"transitivity": 1},
        ),  # 1 -> 3
        (
            simple_gfs[2]["gene_2_0"],
            simple_gfs[3]["gene_3_0"],
            {"transitivity": 0},
        ),  # 2 -> 3
        (
            simple_gfs[2]["gene_2_0"],
            simple_gfs[4]["gene_4_0"],
            {"transitivity": 1},
        ),  # 2 -> 4
        (
            simple_gfs[3]["gene_3_0"],
            simple_gfs[4]["gene_4_0"],
            {"transitivity": 0},
        ),  # 3 -> 4
        (
            simple_gfs[3]["gene_3_0"],
            simple_gfs[5]["gene_5_0"],
            {"transitivity": 1},
        ),  # 3 -> 5
        (
            simple_gfs[4]["gene_4_0"],
            simple_gfs[5]["gene_5_0"],
            {"transitivity": 0},
        ),  # 4 -> 5
    ]
    G.add_edges_from(expected_edges)

    assert G.edges == genes_graph.edges
    assert G.nodes == genes_graph.nodes
    assert all(G[u][v] == genes_graph[u][v] for u, v in G.edges)  # ensure the transitivity attributes match


def test_compute_gene_components(simple_gfs):
    """
    Tests that the function correctly identifies the components of genes that exist
    within a given window size distance of each other.
    """
    model_genes = {
        gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:4])
    }  # genes of the unit model families present in the organism
    window_size = 2
    assert compute_gene_components(model_genes, window_size) == [
        [gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:6])]
    ]


def test_project_unit_on_organisms(simple_gfs, simple_fu, simple_pangenome):
    """Tests that the function correctly projects a system unit onto a set of organisms"""
    unit = SystemUnit(
        functional_unit=simple_fu,
        source="source1",
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo={gf: ("source1", 1) for gf in simple_gfs[:4]},
    )
    model_genes = {gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:4])}
    components = [[gf[f"gene_{i}_0"] for i, gf in enumerate(simple_gfs[:6])]]
    org_proj_name = [
        "unit_name",
        "unit_ID",
        "organism_name",
        "gf_name",
        "partition",
        "annotation",
        "secondary_annotation",
        "gene_ID",
        "local_ID",
        "contig",
        "start",
        "stop",
        "strand",
        "fragment",
        "category",
        "genomic_organization",
        "completeness",
        "product",
    ]
    org_proj_line = namedtuple("OrgLine", org_proj_name)
    expected = [
        [
            "fu",
            "1",
            "org_0",
            "GF0",
            "persistent",
            "protein0",
            "",
            "gene_0_0",
            "",
            "contig_0",
            "100",
            "200",
            "+",
            "False",
            "model",
            "strict",
            "1.0",
            "",
        ],
        [
            "fu",
            "1",
            "org_0",
            "GF1",
            "persistent",
            "protein1",
            "",
            "gene_1_0",
            "",
            "contig_0",
            "200",
            "300",
            "+",
            "False",
            "model",
            "strict",
            "1.0",
            "",
        ],
        [
            "fu",
            "1",
            "org_0",
            "GF2",
            "persistent",
            "protein2",
            "",
            "gene_2_0",
            "",
            "contig_0",
            "300",
            "400",
            "+",
            "False",
            "model",
            "strict",
            "1.0",
            "",
        ],
        [
            "fu",
            "1",
            "org_0",
            "GF3",
            "persistent",
            "protein3",
            "",
            "gene_3_0",
            "",
            "contig_0",
            "400",
            "500",
            "+",
            "False",
            "model",
            "strict",
            "1.0",
            "",
        ],
        [
            "fu",
            "1",
            "org_0",
            "GF4",
            "persistent",
            "",
            "",
            "gene_4_0",
            "",
            "contig_0",
            "500",
            "600",
            "+",
            "False",
            "context",
            "strict",
            "1.0",
            "",
        ],
        [
            "fu",
            "1",
            "org_0",
            "GF5",
            "persistent",
            "",
            "",
            "gene_5_0",
            "",
            "contig_0",
            "600",
            "700",
            "+",
            "False",
            "context",
            "strict",
            "1.0",
            "",
        ],
    ]
    expected_res = [org_proj_line(*map(str, exp)) for exp in expected]
    assert project_unit_on_organisms(components, unit, model_genes) == expected_res


def test_unit_projection(simple_gfs, simple_fu, simple_gf2fam, simple_pangenome):
    """Tests that the function correctly projects a system unit onto the pangenome and organisms."""
    unit = SystemUnit(
        functional_unit=simple_fu,
        source="source1",
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo={gf: ("source1", 1) for gf in simple_gfs[:4]},
    )
    indexes = (
        simple_pangenome.compute_org_bitarrays()
    )  # equivalent to: indexes = {gf: i for i, gf in enumerate(simple_gfs)}
    #                org.mk_bitarray(indexes)
    pangenome_proj, org_proj = unit_projection(unit, simple_gf2fam, indexes)

    # Expected pangenome projection DataFrame
    expected_pangenome = pd.DataFrame(
        {
            0: ["fu"] * 3,  # Functional unit name
            1: ["org_0", "org_1", "org_2"],
            2: ["GF0,GF1,GF2,GF3"] * 3,  # Model families
            3: ["GF4,GF5"] * 3,  # Context families
            4: ["persistent"] * 3,  # Partition
            5: [1.0] * 3,  # Completeness
            6: [6] * 3,  # Strict count
            7: [0] * 3,  # Split count
            8: [0] * 3,  # Extended count
        }
    )

    # Expected organisms projection DataFrame
    expected_organisms = pd.DataFrame(
        {
            "unit_name": ["fu"] * 18,
            "unit_ID": ["1"] * 18,
            "organism_name": ["org_0"] * 6 + ["org_1"] * 6 + ["org_2"] * 6,
            "gf_name": ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"] * 3,
            "partition": ["persistent"] * 18,
            "annotation": ["protein0", "protein1", "protein2", "protein3", "", ""] * 3,
            "secondary_annotation": [""] * 18,
            "gene_ID": [f"gene_{i}_{org}" for org in range(3) for i in range(6)],
            "local_ID": [""] * 18,
            "contig": [f"contig_{org}" for org in range(3) for i in range(6)],
            "start": ["100", "200", "300", "400", "500", "600"] * 3,
            "stop": ["200", "300", "400", "500", "600", "700"] * 3,
            "strand": ["+"] * 18,
            "fragment": ["False"] * 18,
            "category": ["model", "model", "model", "model", "context", "context"] * 3,
            "genomic_organization": ["strict"] * 18,
            "completeness": ["1.0"] * 18,
            "product": [""] * 18,
        }
    )

    pangenome_proj = pangenome_proj.sort_values(by=[1]).reset_index(drop=True)  # Sort by organism
    org_proj = org_proj.sort_values(by=["organism_name", "gf_name"]).reset_index(
        drop=True
    )  # Sort by organism, then gene family
    print(org_proj.columns, expected_organisms)
    pdt.assert_frame_equal(pangenome_proj, expected_pangenome)
    pdt.assert_frame_equal(org_proj, expected_organisms)


def test_system_projection(simple_gfs, single_unit_model, simple_gf2fam, simple_pangenome):
    """Tests that the function correctly projects a system onto the pangenome and organisms."""
    system = System(model=single_unit_model, source="source1", system_id=1)
    indexes = simple_pangenome.compute_org_bitarrays()
    pangenome_proj, org_proj = system_projection(system, indexes, simple_gf2fam)

    assert pangenome_proj.empty
    assert org_proj.empty

    unit = SystemUnit(
        functional_unit=single_unit_model.get("fu"),
        source="source1",
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo={gf: ("source1", 1) for gf in simple_gfs[:4]},
    )
    system.add_unit(unit)
    pangenome_proj, org_proj = system_projection(system, indexes, simple_gf2fam)

    expected_pangenome = pd.DataFrame(
        {
            0: ["1"] * 3,  # System ID
            1: ["TestModel"] * 3,  # Model name
            2: ["fu"] * 3,  # Functional unit name
            3: ["org_0", "org_1", "org_2"],  # Organism names
            4: ["GF0,GF1,GF2,GF3"] * 3,  # Model families
            5: ["GF4,GF5"] * 3,  # Context families
            6: ["persistent"] * 3,  # Partition
            7: [1.0] * 3,  # Completeness
            8: [6] * 3,  # Strict count
            9: [0] * 3,  # Split count
            10: [0] * 3,  # Extended count
        }
    )

    expected_organisms = pd.DataFrame(
        {
            0: ["1"] * 18,  # System ID
            1: ["TestModel"] * 18,  # Model name
            2: ["fu"] * 18,  # Functional unit name
            3: ["1"] * 18,  # System number
            4: ["org_0"] * 6 + ["org_1"] * 6 + ["org_2"] * 6,  # Organism names
            5: ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"] * 3,  # Gene families
            6: ["persistent"] * 18,  # Partition
            7: ["protein0", "protein1", "protein2", "protein3", "", ""] * 3,  # Protein names
            8: [""] * 18,  # Secondary names
            9: [f"gene_{i}_{org}" for org in range(3) for i in range(6)],  # Gene IDs
            10: [""] * 18,  # Gene names
            11: [f"contig_{org}" for org in range(3) for i in range(6)],  # Contig names
            12: ["100", "200", "300", "400", "500", "600"] * 3,  # Start positions
            13: ["200", "300", "400", "500", "600", "700"] * 3,  # Stop positions
            14: ["+"] * 18,  # Strand
            15: ["False"] * 18,  # Is fragment
            16: ["model", "model", "model", "model", "context", "context"] * 3,  # Category
            17: ["strict"] * 18,  # Genomic organization
            18: ["1.0"] * 18,  # Completeness
            19: [""] * 18,  # Product
        }
    )
    pangenome_proj_sorted = pangenome_proj.sort_values(by=[3]).reset_index(drop=True)  # Sort by organism
    org_proj_sorted = org_proj.sort_values(by=[4, 5]).reset_index(drop=True)  # Sort by organism, then gene family

    pdt.assert_frame_equal(pangenome_proj_sorted, expected_pangenome)
    pdt.assert_frame_equal(org_proj_sorted, expected_organisms)


def test_project_pangenome_systems(simple_gfs, single_unit_model, simple_pangenome):
    """Tests that the function correctly projects all systems within a pangenome."""
    unit = SystemUnit(
        functional_unit=single_unit_model.get("fu"),
        source="source1",
        gene_families=set(simple_gfs[:6]),
        families_to_metainfo={gf: ("source1", 1) for gf in simple_gfs[:4]},
    )
    system = System(model=single_unit_model, source="source1", units={unit}, system_id=1)
    simple_pangenome.add_system(system)
    pangenome_proj, org_proj = project_pangenome_systems(simple_pangenome, "source1")

    expected_pangenome = pd.DataFrame(
        {
            "system number": ["1"] * 3,  # System ID
            "system name": ["TestModel"] * 3,  # Model name
            "functional unit name": ["fu"] * 3,  # Functional unit name
            "organism": ["org_0", "org_1", "org_2"],  # Organism names
            "model_GF": ["GF0,GF1,GF2,GF3"] * 3,  # Model families
            "context_GF": ["GF4,GF5"] * 3,  # Context families
            "partition": ["persistent"] * 3,  # Partition
            "completeness": [1.0] * 3,  # Completeness
            "strict": [6] * 3,  # Strict count
            "split": [0] * 3,  # Split count
            "extended": [0] * 3,  # Extended count
        }
    )

    expected_organisms = pd.DataFrame(
        {
            "system number": ["1"] * 18,
            "system name": ["TestModel"] * 18,
            "functional unit name": ["fu"] * 18,
            "subsystem number": ["1"] * 18,
            "organism": ["org_0"] * 6 + ["org_1"] * 6 + ["org_2"] * 6,
            "gene family": ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"] * 3,
            "partition": ["persistent"] * 18,
            "annotation": ["protein0", "protein1", "protein2", "protein3", "", ""] * 3,  # Changed from 'protein'
            "secondary_names": [""] * 18,  # Changed from 'secondary_name'
            "gene.ID": [f"gene_{i}_{org}" for org in range(3) for i in range(6)],  # Changed from 'gene_name'
            "gene.name": [""] * 18,  # Changed from 'gene_id'
            "contig": [f"contig_{org}" for org in range(3) for i in range(6)],
            "start": ["100", "200", "300", "400", "500", "600"] * 3,
            "stop": ["200", "300", "400", "500", "600", "700"] * 3,
            "strand": ["+"] * 18,
            "is_fragment": ["False"] * 18,
            "category": ["model", "model", "model", "model", "context", "context"] * 3,
            "genomic organization": ["strict"] * 18,
            "completeness": ["1.0"] * 18,
            "product": [""] * 18,
        }
    )

    pangenome_proj = pangenome_proj.sort_values(by=["organism"]).reset_index(drop=True)  # Sort by organism
    org_proj = org_proj.sort_values(by=["organism", "gene family"]).reset_index(
        drop=True
    )  # Sort by organism, then gene family

    pdt.assert_frame_equal(pangenome_proj, expected_pangenome)
    pdt.assert_frame_equal(org_proj, expected_organisms)


def test_get_org_df():
    """Tests the org_df is modified as expected"""
    organisms_df = pd.DataFrame(
        {
            "system number": ["1"] * 6 + ["2"],
            "system name": ["TestModel"] * 6 + ["TestModel"],
            "functional unit name": ["fu"] * 6 + ["fu"],
            "subsystem number": ["1"] * 6 + ["2"],
            "organism": ["org_0"] * 7,
            "gene family": ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"] + ["GF5"],
            "partition": ["persistent"] * 6 + ["persistent"],
            "annotation": ["protein0", "protein1", "protein2", "protein3", "", ""] + [""],
            "secondary_names": [""] * 7,
            "gene.ID": [f"gene_{i}_0" for i in range(6)] + ["gene_5_0"],
            "gene.name": [""] * 7,
            "contig": ["contig_0"] * 7,
            "start": ["100", "200", "300", "400", "500", "600"] + ["600"],
            "stop": ["200", "300", "400", "500", "600", "700"] + ["700"],
            "strand": ["+"] * 6 + ["+"],
            "is_fragment": ["False"] * 6 + ["False"],
            "category": ["model", "model", "model", "model", "context", "context"] + ["context"],
            "genomic organization": ["strict"] * 6 + ["strict"],
            "completeness": ["1.0"] * 6 + ["0.8"],
            "product": [""] * 7,
        }
    )

    df, org = get_org_df(organisms_df)

    # ensure that GF5 corresponding to both systems 1 and 2 is properly handled
    expected_df = pd.DataFrame(
        {
            "system number": ["1", "1", "1", "1", "1", "1, 2"],
            "system name": ["TestModel"] * 6,
            "functional unit name": ["fu"] * 6,
            "subsystem number": ["1", "1", "1", "1", "1", "1, 2"],
            "gene family": ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"],
            "partition": ["persistent"] * 6,
            "annotation": ["protein0", "protein1", "protein2", "protein3", "", ""],
            "secondary_names": [""] * 6,
            "gene.ID": [
                "gene_0_0",
                "gene_1_0",
                "gene_2_0",
                "gene_3_0",
                "gene_4_0",
                "gene_5_0",
            ],
            "gene.name": [""] * 6,
            "contig": ["contig_0"] * 6,
            "start": ["100", "200", "300", "400", "500", "600"],
            "stop": ["200", "300", "400", "500", "600", "700"],
            "strand": ["+", "+", "+", "+", "+", "+"],
            "is_fragment": ["False", "False", "False", "False", "False", "False"],
            "category": [
                "model",
                "model",
                "model",
                "model",
                "context",
                "context, context",
            ],
            "genomic organization": ["strict"] * 5 + ["strict, strict"],
            "completeness": ["1.0"] * 5 + ["0.8, 1.0"],
            "product": [""] * 6,
        }
    )

    assert org == "org_0"
    pdt.assert_frame_equal(df, expected_df)


def test_write_projection_systems(simple_orgs):
    """
    Tests that the function correctly writes the pangenome and organism projections to expected files.
    Creates a temporary directory where projection files are written and verified, then removes it.
    """
    pangenome_proj = pd.DataFrame(
        {
            "system number": ["1"] * 3,  # System ID
            "system name": ["TestModel"] * 3,  # Model name
            "functional unit name": ["fu"] * 3,  # Functional unit name
            "organism": ["org_0", "org_1", "org_2"],  # Organism names
            "model_GF": ["GF0,GF1,GF2,GF3"] * 3,  # Model families
            "context_GF": ["GF4,GF5"] * 3,  # Context families
            "partition": ["persistent"] * 3,  # Partition
            "completeness": [1.0] * 3,  # Completeness
            "strict": [6] * 3,  # Strict count
            "split": [0] * 3,  # Split count
            "extended": [0] * 3,  # Extended count
        }
    )

    org_proj = pd.DataFrame(
        {
            "system number": ["1"] * 18,
            "system name": ["TestModel"] * 18,
            "functional unit name": ["fu"] * 18,
            "subsystem number": ["1"] * 18,
            "organism": ["org_0"] * 6 + ["org_1"] * 6 + ["org_2"] * 6,
            "gene family": ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"] * 3,
            "partition": ["persistent"] * 18,
            "annotation": ["protein0", "protein1", "protein2", "protein3", "", ""] * 3,  # Changed from 'protein'
            "secondary_names": [""] * 18,  # Changed from 'secondary_name'
            "gene.ID": [f"gene_{i}_{org}" for org in range(3) for i in range(6)],  # Changed from 'gene_name'
            "gene.name": [""] * 18,  # Changed from 'gene_id'
            "contig": [f"contig_{org}" for org in range(3) for i in range(6)],
            "start": ["100", "200", "300", "400", "500", "600"] * 3,
            "stop": ["200", "300", "400", "500", "600", "700"] * 3,
            "strand": ["+"] * 18,
            "is_fragment": ["False"] * 18,
            "category": ["model", "model", "model", "model", "context", "context"] * 3,
            "genomic organization": ["strict"] * 18,
            "completeness": ["1.0"] * 18,
            "product": [""] * 18,
        }
    )

    import shutil
    from pathlib import Path

    path = Path("pytest-projection-test/")

    try:
        write_projection_systems(path, pangenome_proj, org_proj, simple_orgs)

        # Check that the expected files were created
        assert (path / "systems.tsv").exists()
        assert (path / "projection").exists()

        # Read and verify systems.tsv
        pangenome_file = pd.read_csv(path / "systems.tsv", sep="\t").fillna("")

        assert len(pangenome_file) == 1
        assert pangenome_file["system number"].iloc[0] == 1
        assert pangenome_file["system name"].iloc[0] == "TestModel"
        assert pangenome_file["functional unit name"].iloc[0] == "fu"
        organism_value = str(pangenome_file["organism"].iloc[0])
        assert organism_value == "org_0, org_1, org_2"

        # Read and verify individual organism projection files: org_0.tsv, org_1.tsv, org_2.tsv
        for i, org_name in enumerate(["org_0", "org_1", "org_2"]):
            org_file_path = path / "projection" / f"{org_name}.tsv"
            assert org_file_path.exists(), f"File {org_name}.tsv should exist"

            org_file = pd.read_csv(org_file_path, sep="\t").fillna("")

            # Each organism file should have 6 rows (one per GF)
            assert len(org_file) == 6, f"Expected 6 rows for {org_name}, got {len(org_file)}"

            # Check basic structure
            assert all(org_file["system number"] == 1)
            assert all(org_file["system name"] == "TestModel")
            assert all(org_file["functional unit name"] == "fu")
            assert all(org_file["subsystem number"] == 1)

            # Check gene families are present
            expected_families = ["GF0", "GF1", "GF2", "GF3", "GF4", "GF5"]
            assert set(org_file["gene family"]) == set(expected_families)

            # Check gene IDs match expected pattern
            expected_gene_ids = [f"gene_{j}_{i}" for j in range(6)]
            assert set(org_file["gene.ID"]) == set(expected_gene_ids)

            # Check contigs
            assert all(org_file["contig"] == f"contig_{i}")

            # Check categories
            model_rows = org_file[org_file["category"] == "model"]
            context_rows = org_file[org_file["category"] == "context"]
            assert len(model_rows) == 4  # GF0-GF3
            assert len(context_rows) == 2  # GF4-GF5

            # Check annotations for model genes
            model_annotations = model_rows["annotation"].tolist()
            expected_annotations = ["protein0", "protein1", "protein2", "protein3"]
            assert set(model_annotations) == set(expected_annotations)

            # Check context genes have empty annotations
            context_annotations = context_rows["annotation"].tolist()
            assert all(ann == "" for ann in context_annotations)

    finally:
        # Remove the test directory after the test
        if path.exists():
            shutil.rmtree(path)
