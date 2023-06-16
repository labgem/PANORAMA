from panorama.compare import context
import pytest

def test_create_metanodes_simple():
    # Family cluster 100 is found only in g1
    g1_node_2_family_cluster = {'G1': 5,  'G3': 3, 'G4': 1, 
                                'G5': 2, 'G6': 4, 
                                'G2':100} 
    # Family cluster 200 is found only in g2 
    g2_node_2_family_cluster = {'R1': 1, 'R2': 2, 'R3': 3,
                                'R5': 5, 'R7': 4, "R4":200}

    
    meta_nodes_and_attr, _, _ = context.create_metanodes(g1_node_2_family_cluster,
                                          g2_node_2_family_cluster)
    

    meta_nodes = [m_n for m_n, _ in meta_nodes_and_attr]

    # family cluster 100 and 200 are only found in one graph
    # so it is not keep in the metanodes 
    assert set(meta_nodes) == {"1","2","3","4","5"}


def test_create_metanodes_cluster_found_twice():


    # gene family G5 and G3 have the same cluster
    g1_node_2_family_cluster = {'G1': 5,  'G3': 5, 'G4': 1, 
                                'G5': 2, 'G6': 4}
    
    g2_node_2_family_cluster = {'R1': 1, 'R2': 2, 'R3': 3,
                                'R5': 5, 'R7': 4}

    
    meta_nodes_and_attr, _, _ = context.create_metanodes(g1_node_2_family_cluster,
                                          g2_node_2_family_cluster)
    

    meta_nodes = [m_n for m_n, _ in meta_nodes_and_attr]
    # family cluster is represented twice in metanodes 
    # because it includes two gene families in g1
    assert set(meta_nodes) == {"1", "2", "4", "5", "5´"}