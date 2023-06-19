from panorama.compare import context
import pytest
import networkx as nx


def test_create_metanodes_simple():
    # Family cluster 100 is found only in g1
    g1_node_2_family_cluster = {'G1': 5,  'G3': 3, 'G4': 1, 
                                'G5': 2, 'G6': 4, 
                                'G2':100} 
    # Family cluster 200 is found only in g2 
    g2_node_2_family_cluster = {'R1': 1, 'R2': 2, 'R3': 3,
                                'R5': 5, 'R7': 4, "R4":200}

    
    meta_nodes_2_attr, _, _ = context.create_metanodes(g1_node_2_family_cluster,
                                          g2_node_2_family_cluster)
    

    # family cluster 100 and 200 are only found in one graph
    # so it is not keep in the metanodes 
    assert set(meta_nodes_2_attr) == {"1","2","3","4","5"}


def test_create_metanodes_cluster_found_twice():


    # gene family G5 and G3 have the same cluster
    g1_node_2_family_cluster = {'G1': 5,  'G3': 5, 'G4': 1, 
                                'G5': 2, 'G6': 4}
    
    g2_node_2_family_cluster = {'R1': 1, 'R2': 2, 'R3': 3,
                                'R5': 5, 'R7': 4}

    
    meta_nodes_2_attr, _, _ = context.create_metanodes(g1_node_2_family_cluster,
                                          g2_node_2_family_cluster)
    
    # family cluster is represented twice in metanodes 
    # because it includes two gene families in g1
    assert set(meta_nodes_2_attr) == {"1", "2", "4", "5", "5´"}

@pytest.fixture
def simple_graph():
    G = nx.Graph()
    G.add_edges_from([("G1", "G2"),
                      ("G1", "G3"), 
                      ("G4", "G7"),
                      ("G2", "G4")])
    return G


def test_get_multigraph_edges_simple(simple_graph):
    """
    Test case for translating edges into metanodes edges.
    """

    node2metanode = {"G1":["A"],
                     "G2":["B"],
                     "G3":["C"],
                     "G4":["D"]}

    g_multigraph_edges = context.get_multigraph_edges(simple_graph, node2metanode)

    assert sorted(g_multigraph_edges) == sorted([("A", "B"),
                                                 ("A", "C"),
                                                 ("B", "D"),
                                                 ])


def test_get_multigraph_edges_multiple_metanodes(simple_graph):
    """
    Test case for translating edges when a node is associated to multiple metanodes.

    This test verifies that edges from nodes to metanodes are correctly translated, 
    including supplementary edges for additional metanodes associated with a node.
    """
    
    node2metanode = {"G1":["A"],
                     "G2":["B"],
                     "G3":["C"],
                     "G4":["D", "D_prime"]}

    g_multigraph_edges = context.get_multigraph_edges(simple_graph, node2metanode)

    assert sorted(g_multigraph_edges) == sorted([("A", "B"),
                                                 ("A", "C"),
                                                 ("B", "D"),
                                                 ("B", "D_prime"),
                                                 ])

@pytest.fixture
def r_graph():

    r_edges = [('R1', 'R2'), 
                ('R2', 'R3'), 
                ('R2', 'R4'), 
                ('R4', 'R6'),
                ('R3', 'R5'), 
                ('R6', 'R5'), 
                ('R5', 'R7')]

    r_graph = nx.Graph()
    r_graph.add_edges_from(r_edges)
    return r_graph


@pytest.fixture
def g_graph():

    g_edges = [('G1', 'G2'), 
               ('G2', 'G3'),
               ('G3', 'G4'),
               ('G4', 'G5'),
               ('G5', 'G6')]

    g_graph = nx.Graph()
    g_graph.add_edges_from(g_edges)
    return g_graph


def test_get_conserved_genomics_contexts(r_graph, g_graph):
    """
    Test case to get CCC from two graphs.

    This simple test dataset is taken from Boyer et al. article
    (https://doi.org/10.1093/bioinformatics/bti711).
    """

    node_2_cluster_node = {'R1': 1, 'R2': 2, 'R3': 3, 'R5': 5, 'R7': 4,
                              'G1': 5, 'G3': 3, 'G4': 1, 'G5': 2, 'G6': 4}

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph,
                                            g_graph,
                                            node_2_cluster_node, min_cgc_size=1)

    expected_results = [({"R1", "R2", "R3"}, {"G4", "G5", "G3"}), 
                        ({"R7"}, {"G6"}),
                        ({"R5"}, {"G1"})]

    assert sorted(ccc_results, key=lambda x : len(x[0]) ) == sorted(expected_results, key=lambda x : len(x[0]) )

def test_get_conserved_genomics_contexts_duplicated_clstr(r_graph, g_graph):
    """
    Test case where two nodes of the same graph belong to the same cluster.
    """

    # R1 and R2 belong to the same cluster
    node_2_cluster_node = {'R1': 1, 'R2': 1, 'R3': 3, 'R5': 5, 'R7': 4,
                              'G1': 5, 'G3': 3, 'G4': 1, 'G5': 2, 'G6': 4}

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph,
                                            g_graph,
                                            node_2_cluster_node, min_cgc_size=1)
    # G5 is no more found in conserved connected components
    #  as no R node belong to cluster 2

    expected_results = [({"R1", "R2", "R3"}, {"G4", "G3"}), 
                        ({"R7"}, {"G6"}),
                        ({"R5"}, {"G1"})]
                        
    assert sorted(ccc_results, key=lambda x : len(x[0]) ) == sorted(expected_results, key=lambda x : len(x[0]) )
    
    ccc_results_min_size_3, _ = context.get_conserved_genomics_contexts(r_graph,
                                            g_graph,
                                            node_2_cluster_node, min_cgc_size=3)
    # G5 is no more found in conserved connected components
    #  as no R node belong to cluster 2

    expected_results = [({"R1", "R2", "R3"}, {"G4", "G3"})]
                        
    assert sorted(ccc_results_min_size_3, key=lambda x : len(x[0]) ) == sorted(expected_results, key=lambda x : len(x[0]) )

def test_get_conserved_genomics_contexts_extreme_clstr(r_graph, g_graph):
    """
    Test case where all nodes of the same graph belong to the same cluster.
    """

    # R1 and R2 belong to the same cluster
    node_2_cluster_node = {'R1': 1, 'R2': 1, 'R3': 1, 'R5': 1, 'R7': 1,
                              'G1': 5, 'G3': 3, 'G4': 1, 'G5': 2, 'G6': 4}

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph,
                                            g_graph,
                                            node_2_cluster_node, min_cgc_size=1)
    # All R belong to cluster 1. Only G4 belong to cluster 1
    expected_results = [({"R1"}, {"G4"}),
                        ({"R2"}, {"G4"}),
                        ({"R3"}, {"G4"}),
                        ({"R5"}, {"G4"}),
                        ({"R7"}, {"G4"}),
                            ]
                        
    assert sorted(ccc_results, key=lambda x : len(x[0]) ) == sorted(expected_results, key=lambda x : len(x[0]) )


def test_get_conserved_genomics_no_shared_clstr(r_graph, g_graph):
    """
    Test case where graph share no cluster family.
    """

    # not shared clusters
    node_2_cluster_node = {'R1': 1, 'R2': 2, 'R3': 3, 'R5': 4, 'R7': 5,
                              'G1': 15, 'G3': 13, 'G4': 11, 'G5': 12, 'G6': 14}

    ccc_results, _ = context.get_conserved_genomics_contexts(r_graph,
                                            g_graph,
                                            node_2_cluster_node)
                        
    assert ccc_results ==  []


def test_get_conserved_genomics_with_min_cgc(r_graph, g_graph):
    """
    Test case where share cluster does not reach minimum size of context
    """

    # share only two cluster

    node_2_cluster_node = {'R1': 1, 'R2': 2, 'R3': 3, 'R5': 5, 'R7': 4,
                           'G1': 15, 'G3': 13, 'G4': 1, 'G5': 2, 'G6': 14}

    ccc_results_min_size_2, _ = context.get_conserved_genomics_contexts(r_graph,
                                            g_graph,
                                            node_2_cluster_node, min_cgc_size=2)

                        
    assert ccc_results_min_size_2 == [({"R1", "R2"}, {"G4", "G5"})]

    ccc_results_min_size_3, _ = context.get_conserved_genomics_contexts(r_graph,
                                        g_graph,
                                        node_2_cluster_node, min_cgc_size=3)

    assert ccc_results_min_size_3 == []

def test_get_connected_components():
    
    nodes = ["A", "B", "C"]
    edges = [("A", "C")]
    cc_result = list(context.get_connected_components(nodes, edges))

    assert sorted(cc_result) == [{"A", "C"}, {"B"}]