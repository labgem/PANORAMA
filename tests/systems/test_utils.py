import pytest

from ppanggolin.meta.meta import assign_metadata
from ppanggolin.pangenome import Pangenome

from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family

from panorama.systems.utils import (get_metadata_to_families, dict_families_context, 
                                    get_gfs_matrix_combination, filter_global_context, 
                                    filter_local_context, check_for_families, 
                                    check_needed_families)

from collections import defaultdict
import pandas as pd
import networkx as nx



@pytest.fixture()
def simple_gfs():
    return [GeneFamily(family_id=i, name=f'GF{i}') for i in range(10)]

@pytest.fixture()
def simple_pangenome(simple_gfs):
    pangenome = Pangenome()
    for gf in simple_gfs:
        pangenome.add_gene_family(gf)
    # GF_i -> protein_i
    metadata = pd.DataFrame({
        'families': [gf.name for gf in simple_gfs],
        'protein_name': [f'protein{i}' for i in range(10)],
        'score': [1.0 for _ in range(10)],
        })
    assign_metadata(metadata, pangenome, source='source1', metatype='families')
    return pangenome

@pytest.fixture()
def simple_model():
    mandatory_gfs = {Family(name=f'protein{i}', presence='mandatory') for i in range(3)} 
    accessory_gfs = {Family(name=f'protein{i+3}', presence='accessory') for i in range(3)}
    fu = FuncUnit(mandatory=mandatory_gfs, accessory=accessory_gfs, min_mandatory=2, min_total=4,transitivity=1)
    model = Model(name='TestModel', mandatory={fu})
    return model



def test_get_metadata_to_families(simple_pangenome):
    annot2fam = get_metadata_to_families(simple_pangenome, sources=["source1"])
    gfs = [gf for gf in simple_pangenome.gene_families]
    
    assert annot2fam == {'source1': defaultdict(set, 
                        {f'protein{i}': {gfs[i]} for i in range(len(gfs))})} # mapping from protein annotations to GFs
    
    
    extra_gfs = [GeneFamily(family_id=i, name=f'GF{i}') for i in range(10,18)]
    gfs.extend(extra_gfs)
    # Add multiple annotations per protein: GF_i -> protein_i 
    #                                       GF_i+10 -> protein_i
    families = [gf.name for gf in gfs[:10]] + [gf.name for gf in gfs[10:18]]  
    protein_names = [f'protein{i}' for i in range(10)] + [f'protein{i}' for i in range(8)]
    # Add mutiple annotations for some families: GF0 -> protein0, protein8
    #                                            GF1 -> protein1, protein9                                             
    families += [gfs[0].name, gfs[1].name]
    protein_names += ['protein8', 'protein9']
    
    for gf in extra_gfs:
        simple_pangenome.add_gene_family(gf)
    metadata = pd.DataFrame({
        'families': families,
        'protein_name': protein_names,
        })
    metadata['score'] = 1.0
    # assign updated metadata as source2
    assign_metadata(metadata, simple_pangenome, source='source2', metatype='families')
    
    assert get_metadata_to_families(simple_pangenome, sources=["source2"]) == {'source2': defaultdict(set, 
                                                                              {f'protein{i}': {gfs[i], gfs[i+10]} 
                                                                              for i in range(8)} | {f'protein{i}': 
                                                                              {gfs[i], gfs[i-8]} for i in range(8,10)})}


def test_dict_families_context(simple_model, simple_gfs):
    annot2fam = {'source1': defaultdict(set, {f'protein{i}': {simple_gfs[i]} for i in range(len(simple_gfs))})} # protein_i -> {GF_i}
    fu = next(simple_model.func_units)
    result = dict_families_context(simple_model, annot2fam)
    # mandatory: protein0, protein1, protein2
    # accessory: protein3, protein4, protein5
    assert result[0] == set(simple_gfs[:6]) # all GFs corresponding to model families
    assert result[1] == defaultdict(set, {simple_gfs[i]: 
                                    {[f for f in simple_model.families if f.name == f'protein{i}'][0]} 
                                    for i in range(6)}) # mapping from GFs to model families
    assert result[2] == {f'protein{i}': 'source1' for i in range(6)} # annotation source mapping
    
    
    # Add many-to-many GF-protein associations
    simple_gfs.extend([GeneFamily(family_id=i, name=f'GF{i}') for i in range(10,18)])
    annot2fam = {'source2': defaultdict(set, {f'protein{i}': {simple_gfs[i], simple_gfs[i+10]}  # protein_i -> {GF_i, GF_i+10}
                                              for i in range(8)} | {f'protein{i}': {simple_gfs[i], simple_gfs[i-8]} # protein9 -> {GF1, GF9}; protein8 -> {GF0, GF8}
                                              for i in range(8,10)})}                                               # protein1 -> {GF1, GF11}; protein0 -> {GF0, GF10}        
    # Add new family to model
    new_family = Family(name='protein8', presence='mandatory')
    fu.mandatory.add(new_family)
    
    model_fams, gf2fam, fam2source = dict_families_context(simple_model, annot2fam)
    
    assert model_fams == set(simple_gfs[:6]) | set(simple_gfs[10:16]) | {simple_gfs[8]}
    family_lookup = {f.name: f for f in simple_model.families}
    assert gf2fam == defaultdict(set, 
                                 {simple_gfs[i]: {family_lookup[f'protein{i}']} 
                                 for i in range(1,6)}) | defaultdict(set, 
                                 {simple_gfs[i+10]: {family_lookup[f'protein{i}']} 
                                 for i in range(6)}) | {simple_gfs[0]: 
                                 {family_lookup['protein0'], new_family}} | {simple_gfs[8]: {new_family}}
    assert fam2source == {f'protein{i}': 'source2' for i in range(6)} | {'protein8': 'source2'}
    
    
def test_get_gfs_matrix_combination(simple_gfs):
    gf2fam = defaultdict(set, {simple_gfs[i]: {Family(name=f'protein{i}', presence='mandatory' 
                                                      if i < 3 else 'accessory')} for i in range(6)})
    fu_families = set(simple_gfs[:6])
    result = get_gfs_matrix_combination(fu_families, gf2fam)
    # binary matrix with protein annotations as rows and GFs as columns
    assert result.sort_index(axis=0).sort_index(axis=1).equals(
        pd.DataFrame(
            [[int(i == j) for j in range(6)] for i in range(6)],
            index=[f'protein{i}' for i in range(6)],
            columns=[f'GF{i}' for i in range(6)]
            ).sort_index(axis=0).sort_index(axis=1)
        )
    
    
    # Consider many-to-many GF-protein associations
    simple_gfs.extend([GeneFamily(family_id=i, name=f'GF{i}') for i in range(10,18)])
    fu_families = set(simple_gfs[:6]) | set(simple_gfs[10:16]) | {simple_gfs[8]}
    gf2fam = defaultdict(set, 
                            {simple_gfs[i]: {Family(name=f'protein{i}', presence='mandatory' if i < 3 else 'accessory')} 
                            for i in range(1,6)}) | defaultdict(set, 
                            {simple_gfs[i+10]: {Family(name=f'protein{i}', presence='accessory' if i < 3 else 'accessory')} 
                            for i in range(6)}) | {simple_gfs[0]: {Family(name='protein0', presence='mandatory'), 
                                                    Family(name='protein8', presence='mandatory')}} | {simple_gfs[8]: 
                                                    {Family(name='protein8', presence='mandatory')}}
    result = get_gfs_matrix_combination(fu_families, gf2fam)
    
    annotations = list({family.name for gf in fu_families for family in gf2fam[gf]}) # rows
    model_fams = [gf.name for gf in fu_families] # columns
    gf_lookup = {gf.name: gf for gf in fu_families}
    matrix = []
    for protein in annotations:
        row = []
        for gf_name in model_fams:
            gf = gf_lookup[gf_name]
            found = any(fam.name == protein for fam in gf2fam.get(gf, [])) # Check if gf has protein annotation
            row.append(1 if found else 0)
        matrix.append(row)

    assert result.sort_index(axis=0).sort_index(axis=1).equals(
        pd.DataFrame(matrix, index=annotations, columns=model_fams).sort_index(axis=0).sort_index(axis=1))
    

def test_check_needed_families(simple_model):
    matrix = pd.DataFrame([[int(i == j) for j in range(6)] for i in range(6)],
                            index=[f'protein{i}' for i in range(6)],
                            columns=[f'GF{i}' for i in range(6)])
    fu = next(simple_model.func_units) # the first functional unit of the model
    # min_total = 4, min_mandatory = 2
    # 3 mandatory and 3 accessory families covered => model satisfied
    assert check_needed_families(matrix, fu) is True
    
    # Updating min_total to 7 -> model unsatisfied
    fu.min_total = 7
    assert check_needed_families(matrix, next(simple_model.func_units)) is False


class DummyGeneFamily:
    def __init__(self, name, organisms):
        self.name = name
        self.organisms = organisms

def test_filter_global_context():
    gf1 = DummyGeneFamily("GF1", ["org1", "org2"])
    gf2 = DummyGeneFamily("GF2", ["org1", "org2"])
    gf3 = DummyGeneFamily("GF3", ["org1", "org2", "org3"])

    G = nx.Graph()
    G.add_edge(gf1, gf2, genomes={"org1", "org2"})  # proportion = 1.0 for both
    G.add_edge(gf1, gf3, genomes={"org1"})          # proportion = 0.5 for gf1, 0.333 for gf3

    filtered = filter_global_context(G, jaccard_threshold=0.8)

    # Only edge (gf1, gf2) remains
    assert set(filtered.edges()) == {(gf1, gf2)}
    edge_data = filtered.get_edge_data(gf1, gf2)
    assert edge_data['f1'] == "GF1"
    assert edge_data['f2'] == "GF2"
    assert edge_data['f1_jaccard_gene'] == edge_data['f2_jaccard_gene'] == 1.0
    
def test_filter_local_context():
    gf1 = DummyGeneFamily("GF1", ["org1", "org2"])
    gf2 = DummyGeneFamily("GF2", ["org1", "org2"])
    gf3 = DummyGeneFamily("GF3", ["org1", "org2", "org3"])

    orgs = {"org1"}
        
    G = nx.Graph()
    G.add_edge(gf1, gf2, genomes={"org1", "org2"})  # local proportion = 1.0 for both (when only org1 considered)
    G.add_edge(gf1, gf3, genomes={"org1"})          # local proportion = 1.0 for both

    filtered = filter_local_context(G, orgs, jaccard_threshold=0.8) 

    # both edges remain
    assert set(filtered.edges()) == {(gf1, gf2), (gf1, gf3)}
    edge_data1 = filtered.get_edge_data(gf1, gf2)
    edge_data2 = filtered.get_edge_data(gf1, gf3)
    assert edge_data1['f1_jaccard_gene'] == edge_data1['f2_jaccard_gene'] == 1.0
    assert edge_data2['f1_jaccard_gene'] == edge_data2['f2_jaccard_gene'] == 1.0
    
    orgs = {"org2", "org3"}  # Assume combination in org2 and org3 instead -> local proportion = 0 for both sides of edge (gf1, gf3)
    G.add_edge(gf2, gf3, genomes={"org2"})  # local proportion = 1.0 for gf2, 0.5 for gf3
    filtered = filter_local_context(G, orgs, jaccard_threshold=0.65)
    assert set(filtered.edges()) == {(gf1, gf2)} # only edge (gf1, gf2) remains
    

def test_check_for_families(simple_gfs, simple_model, simple_pangenome):
    # simple_pangenome is needed to assign metadata to families
    gfs_in_cc = simple_gfs[:4] # assume first 4 GFs in a connected component
    gf2fam = defaultdict(set, {simple_gfs[i]: {Family(name=f'protein{i}', presence='mandatory' 
                                                      if i < 3 else 'accessory')} for i in range(6)})
    fam2source = {f'protein{i}': 'source1' for i in range(6)}
    fu = next(simple_model.func_units)
    # 3 mandatory and 1 accessory families covered => model satisfied
    check, gf2meta_info = check_for_families(gfs_in_cc, gf2fam, fam2source, fu)
    assert check == True
    assert gf2meta_info == {gf: ('source1', 1) for gf in gfs_in_cc} # mapping from GFs to (source, best metadata ID)
    
    # Incrementing min_total to 5 
    fu.min_total = 5
    assert check_for_families(gfs_in_cc, gf2fam, fam2source, fu) == (False, {}) # model unsatisfied
    
    