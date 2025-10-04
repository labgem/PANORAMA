import networkx as nx
import pytest

from panorama.compare import context


class MockGeneFamily:
    """Mock gene family object with akin attribute for testing."""

    def __init__(self, name: str, akin: int):
        self.name = name
        self.akin = akin

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"MockGeneFamily(name='{self.name}', akin={self.akin})"

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if isinstance(other, MockGeneFamily):
            return self.name == other.name
        return False


def get_gene_families_by_names(families, names):
    """Helper function to extract multiple gene families by names.

    Args:
        families: Set or collection of gene family objects
        names: List of gene family names to extract

    Returns:
        Dictionary mapping names to gene family objects
    """
    families_dict = {gf.name: gf for gf in families}
    return {name: families_dict[name] for name in names if name in families_dict}


def compare_ccc_results(actual_results, expected_results):
    """Helper function to compare CCC results ignoring order."""
    # Convert to sets for comparison (each element is (set, set))
    actual_set = {(frozenset(r), frozenset(g)) for r, g in actual_results}
    expected_set = {(frozenset(r), frozenset(g)) for r, g in expected_results}
    return actual_set == expected_set


def test_create_metanodes_simple():
    # Family cluster 100 is found only in g1
    g1_node_2_family_cluster = {"G1": 5, "G3": 3, "G4": 1, "G5": 2, "G6": 4, "G2": 100}
    # Family cluster 200 is found only in g2
    g2_node_2_family_cluster = {"R1": 1, "R2": 2, "R3": 3, "R5": 5, "R7": 4, "R4": 200}

    meta_nodes_2_attr, _, _ = context.create_metanodes(g1_node_2_family_cluster, g2_node_2_family_cluster)

    # family cluster 100 and 200 are only found in one graph
    # so it is not keep in the metanodes
    assert set(meta_nodes_2_attr) == {"1", "2", "3", "4", "5"}


def test_create_metanodes_cluster_found_twice():
    # gene family G5 and G3 have the same cluster
    g1_node_2_family_cluster = {"G1": 5, "G3": 5, "G4": 1, "G5": 2, "G6": 4}

    g2_node_2_family_cluster = {"R1": 1, "R2": 2, "R3": 3, "R5": 5, "R7": 4}

    meta_nodes_2_attr, _, _ = context.create_metanodes(g1_node_2_family_cluster, g2_node_2_family_cluster)

    # family cluster is represented twice in metanodes
    # because it includes two gene families in g1
    assert set(meta_nodes_2_attr) == {"1", "2", "4", "5", "5´"}


@pytest.fixture
def simple_graph():
    G = nx.Graph()
    G.add_edges_from([("G1", "G2"), ("G1", "G3"), ("G4", "G7"), ("G2", "G4")])
    return G


def test_get_multigraph_edges_simple(simple_graph):
    """
    Test case for translating edges into metanodes edges.
    """

    node2metanode = {"G1": ["A"], "G2": ["B"], "G3": ["C"], "G4": ["D"]}

    g_multigraph_edges = context.get_multigraph_edges(simple_graph, node2metanode)

    assert sorted(g_multigraph_edges) == sorted(
        [
            ("A", "B"),
            ("A", "C"),
            ("B", "D"),
        ]
    )


def test_get_multigraph_edges_multiple_metanodes(simple_graph):
    """
    Test case for translating edges when a node is associated to multiple metanodes.

    This test verifies that edges from nodes to metanodes are correctly translated,
    including supplementary edges for additional metanodes associated with a node.
    """

    node2metanode = {"G1": ["A"], "G2": ["B"], "G3": ["C"], "G4": ["D", "D_prime"]}

    g_multigraph_edges = context.get_multigraph_edges(simple_graph, node2metanode)

    assert sorted(g_multigraph_edges) == sorted(
        [
            ("A", "B"),
            ("A", "C"),
            ("B", "D"),
            ("B", "D_prime"),
        ]
    )


def create_r_graph_with_akin_mapping(akin_mapping):
    """Helper to create r_graph with custom akin mappings."""
    families = {}
    for name, akin in akin_mapping.items():
        if name.startswith("R"):
            families[name] = MockGeneFamily(name, akin)

    r_edges = [
        (families["R1"], families["R2"]),
        (families["R2"], families["R3"]),
        (
            families["R2"],
            families.get("R4", MockGeneFamily("R4", akin_mapping.get("R4", 200))),
        ),
        (
            families.get("R4", MockGeneFamily("R4", akin_mapping.get("R4", 200))),
            families.get("R6", MockGeneFamily("R6", akin_mapping.get("R6", 6))),
        ),
        (families["R3"], families["R5"]),
        (
            families.get("R6", MockGeneFamily("R6", akin_mapping.get("R6", 6))),
            families["R5"],
        ),
        (families["R5"], families["R7"]),
    ]

    r_graph = nx.Graph()
    r_graph.add_edges_from(r_edges)

    # Add required edge attributes
    for u, v in r_graph.edges():
        r_graph[u][v]["mean_transitivity"] = 0.5  # Default value for testing

    return r_graph


def create_g_graph_with_akin_mapping(akin_mapping):
    """Helper to create g_graph with custom akin mappings."""
    families = {}
    for name, akin in akin_mapping.items():
        if name.startswith("G"):
            families[name] = MockGeneFamily(name, akin)

    g_edges = [
        (
            families["G1"],
            families.get("G2", MockGeneFamily("G2", akin_mapping.get("G2", 100))),
        ),
        (
            families.get("G2", MockGeneFamily("G2", akin_mapping.get("G2", 100))),
            families["G3"],
        ),
        (families["G3"], families["G4"]),
        (families["G4"], families["G5"]),
        (families["G5"], families["G6"]),
    ]

    g_graph = nx.Graph()
    g_graph.add_edges_from(g_edges)

    # Add required edge attributes
    for u, v in g_graph.edges():
        g_graph[u][v]["mean_transitivity"] = 0.5  # Default value for testing

    return g_graph


@pytest.fixture
def r_graph():
    """Original r_graph fixture with standard akin mappings."""
    akin_mapping = {"R1": 1, "R2": 2, "R3": 3, "R4": 200, "R5": 5, "R6": 6, "R7": 4}
    return create_r_graph_with_akin_mapping(akin_mapping)


@pytest.fixture
def g_graph():
    """Original g_graph fixture with standard akin mappings."""
    akin_mapping = {"G1": 5, "G2": 100, "G3": 3, "G4": 1, "G5": 2, "G6": 4}
    return create_g_graph_with_akin_mapping(akin_mapping)


def test_get_conserved_genomics_contexts(r_graph, g_graph):
    """
    Test case to get CCC from two graphs.

    This simple test dataset is taken from Boyer et al. article
    (https://doi.org/10.1093/bioinformatics/bti711).
    """

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph, g_graph, min_cgc_size=1)

    # Extract the gene family objects from the graphs for comparison
    r_families = set(r_graph.nodes())
    g_families = set(g_graph.nodes())

    # Get gene families by names for expected results
    r = get_gene_families_by_names(r_families, ["R1", "R2", "R3", "R5", "R7"])
    g = get_gene_families_by_names(g_families, ["G1", "G3", "G4", "G5", "G6"])

    # Convert results to tuples format for easier comparison
    result_tuples = [(info["graphA_nodes"], info["graphB_nodes"]) for info in ccc_results]
    expected_results = [
        ({r["R1"], r["R2"], r["R3"]}, {g["G4"], g["G5"], g["G3"]}),
        ({r["R7"]}, {g["G6"]}),
        ({r["R5"]}, {g["G1"]}),
    ]

    assert compare_ccc_results(result_tuples, expected_results)


def test_get_conserved_genomics_contexts_duplicated_clstr():
    """
    Test case where two nodes of the same graph belong to the same cluster.
    """

    # R1 and R2 belong to the same cluster (akin=1)
    r_akin_mapping = {"R1": 1, "R2": 1, "R3": 3, "R4": 200, "R5": 5, "R6": 6, "R7": 4}
    g_akin_mapping = {"G1": 5, "G2": 100, "G3": 3, "G4": 1, "G5": 2, "G6": 4}

    r_graph = create_r_graph_with_akin_mapping(r_akin_mapping)
    g_graph = create_g_graph_with_akin_mapping(g_akin_mapping)

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph, g_graph, min_cgc_size=1)

    # Extract families for comparison
    r_families = set(r_graph.nodes())
    g_families = set(g_graph.nodes())

    # Get gene families by names for expected results
    r = get_gene_families_by_names(r_families, ["R1", "R2", "R3", "R5", "R7"])
    g = get_gene_families_by_names(g_families, ["G1", "G3", "G4", "G6"])

    # Convert results to tuples format for easier comparison
    result_tuples = [(info["graphA_nodes"], info["graphB_nodes"]) for info in ccc_results]
    # G5 is no more found in conserved connected components
    # as no R node belong to cluster 2
    expected_results = [
        ({r["R1"], r["R2"], r["R3"]}, {g["G4"], g["G3"]}),
        ({r["R7"]}, {g["G6"]}),
        ({r["R5"]}, {g["G1"]}),
    ]

    assert compare_ccc_results(result_tuples, expected_results)

    ccc_results_min_size_3, _ = context.get_conserved_genomics_contexts(r_graph, g_graph, min_cgc_size=3)

    # Convert results to tuples format for easier comparison
    result_tuples_min_size_3 = [(info["graphA_nodes"], info["graphB_nodes"]) for info in ccc_results_min_size_3]
    # Only large contexts are kept
    expected_results = [({r["R1"], r["R2"], r["R3"]}, {g["G4"], g["G3"]})]

    assert compare_ccc_results(result_tuples_min_size_3, expected_results)


def test_get_conserved_genomics_contexts_extreme_clstr():
    """
    Test case where all nodes of the same graph belong to the same cluster.
    """

    # All R nodes belong to cluster 1
    r_akin_mapping = {"R1": 1, "R2": 1, "R3": 1, "R4": 1, "R5": 1, "R6": 1, "R7": 1}
    g_akin_mapping = {"G1": 5, "G2": 100, "G3": 3, "G4": 1, "G5": 2, "G6": 4}

    r_graph = create_r_graph_with_akin_mapping(r_akin_mapping)
    g_graph = create_g_graph_with_akin_mapping(g_akin_mapping)

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph, g_graph, min_cgc_size=1)

    # Extract families for comparison
    r_families = set(r_graph.nodes())
    g_families = set(g_graph.nodes())

    # Get gene families by names for expected results (get all R nodes that are in the graph)
    r_names = [f"R{i}" for i in [1, 2, 3, 4, 5, 6, 7]]  # All R nodes that should be in graph
    r = get_gene_families_by_names(r_families, r_names)
    g = get_gene_families_by_names(g_families, ["G4"])

    # Convert results to tuples format for easier comparison
    result_tuples = [(info["graphA_nodes"], info["graphB_nodes"]) for info in ccc_results]
    # All R belong to cluster 1. Only G4 belong to cluster 1
    # Include all R nodes that are actually in the graph
    expected_results = [
        ({r["R1"]}, {g["G4"]}),
        ({r["R2"]}, {g["G4"]}),
        ({r["R3"]}, {g["G4"]}),
        ({r["R4"]}, {g["G4"]}),
        ({r["R5"]}, {g["G4"]}),
        ({r["R6"]}, {g["G4"]}),
        ({r["R7"]}, {g["G4"]}),
    ]

    assert compare_ccc_results(result_tuples, expected_results)


def test_get_conserved_genomics_no_shared_clstr():
    """
    Test case where graph share no cluster family.
    """

    # No shared clusters
    r_akin_mapping = {"R1": 1, "R2": 2, "R3": 3, "R4": 4, "R5": 5, "R6": 6, "R7": 7}
    g_akin_mapping = {"G1": 15, "G2": 16, "G3": 13, "G4": 11, "G5": 12, "G6": 14}

    r_graph = create_r_graph_with_akin_mapping(r_akin_mapping)
    g_graph = create_g_graph_with_akin_mapping(g_akin_mapping)

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph, g_graph)

    assert ccc_results == []


def test_get_conserved_genomics_with_min_cgc():
    """
    Test case where share cluster does not reach minimum size of context
    """

    # Share only two clusters
    r_akin_mapping = {"R1": 1, "R2": 2, "R3": 3, "R4": 200, "R5": 5, "R6": 6, "R7": 4}
    g_akin_mapping = {"G1": 15, "G2": 100, "G3": 13, "G4": 1, "G5": 2, "G6": 14}

    r_graph = create_r_graph_with_akin_mapping(r_akin_mapping)
    g_graph = create_g_graph_with_akin_mapping(g_akin_mapping)

    ccc_results_min_size_2, _ = context.get_conserved_genomics_contexts(r_graph, g_graph, min_cgc_size=2)

    # Extract families for comparison
    r_families = set(r_graph.nodes())
    g_families = set(g_graph.nodes())

    # Get gene families by names for expected results
    r = get_gene_families_by_names(r_families, ["R1", "R2"])
    g = get_gene_families_by_names(g_families, ["G4", "G5"])

    # Convert results to tuples format for easier comparison
    result_tuples_min_size_2 = [(info["graphA_nodes"], info["graphB_nodes"]) for info in ccc_results_min_size_2]
    expected_results_min_size_2 = [({r["R1"], r["R2"]}, {g["G4"], g["G5"]})]
    assert compare_ccc_results(result_tuples_min_size_2, expected_results_min_size_2)

    ccc_results_min_size_3, _ = context.get_conserved_genomics_contexts(r_graph, g_graph, min_cgc_size=3)

    assert ccc_results_min_size_3 == []


def test_get_connected_components():
    nodes = ["A", "B", "C"]
    edges = [("A", "C")]
    cc_result = list(context.get_connected_components(nodes, edges))

    assert sorted(cc_result) == [{"A", "C"}, {"B"}]
