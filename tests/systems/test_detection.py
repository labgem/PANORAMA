import pytest
from panorama.systems.detection import (get_functional_unit_gene_families, 
                                        search_unit_in_cc, search_unit_in_context, 
                                        search_unit_in_combination, search_system_units,
                                        search_system, search_systems, search_for_system, search_system_units,
                                        check_for_needed_units,check_for_forbidden_unit,
                                        get_system_unit_combinations)

from ppanggolin.genome import Gene, Organism, Contig

from ppanggolin.meta.meta import assign_metadata
from panorama.pangenomes import Pangenome 
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family
from panorama.systems.system import System, SystemUnit
from panorama.systems.models import Model, Models

from collections import defaultdict
import networkx as nx
import pandas as pd


@pytest.fixture()
def simple_gfs():
    gfs = [GeneFamily(family_id=i, name=f'GF{i}') for i in range(10)]
    org = Organism(name='test_org')
    contig = Contig(identifier=1, name='test_contig')
    # Add a gene to each GF and fill its attributes
    for i, gf in enumerate(gfs):
        gene = Gene(gene_id=f'gene_{i}')
        gene.family = gf
        gf[gene.ID] = gene
        gene.fill_parents(organism = org, contig = contig)
        gene.fill_annotations(position=i, strand='+', start=(i+1)*100, stop=(i+2)*100)
        contig.add(gene)
    return gfs

@pytest.fixture()
def simple_model():
    # Functional Unit 1
    mandatory_gfs = {Family(name=f'protein{i}', presence='mandatory') for i in range(3)} # fu1_mandatory: GF0, GF1, GF2
    accessory_gfs = {Family(name=f'protein{i+3}', presence='accessory') for i in range(3)} # fu1_accessory: GF3, GF4, GF5
    neutral_gfs = {Family(name='protein9', presence='neutral')} # fu1_neutral: GF9
    fu1 = FuncUnit(name='fu1', mandatory=mandatory_gfs, accessory=accessory_gfs, neutral=neutral_gfs, 
                   min_mandatory=2, min_total=4, transitivity=1)
    # Functional Unit 2
    extra_gfs = {Family(name=f'protein{i+6}', presence='accessory') for i in range(3)} # fu2_accessory: GF6, GF7, GF8
    fu2 = FuncUnit(name='fu2', mandatory={Family(name='protein10', presence='mandatory')}, accessory=extra_gfs, 
                   min_mandatory=1, min_total=3, transitivity=1) # fu2_mandatory: GF10
    # Model: fu1 mandatory; fu2 accessory
    model = Model(name='TestModel', mandatory={fu1}, accessory={fu2})
    return model

@pytest.fixture()
def simple_pangenome(simple_gfs):
    pangenome = Pangenome(name='TestPangenome')
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


def test_get_functional_unit_gene_families(simple_model, simple_gfs):
    family_lookup = {f.name: f for f in simple_model.families} # mapping protein names to model Family objects
    gf2fam = defaultdict(set, {simple_gfs[i]: {family_lookup[f'protein{i}']} for i in range(10)})
    
    # Get the (mandatory | accessory) and neutral family sets per functional unit with at least one occurence in the pangenome
    fu1 = simple_model.get('fu1') 
    assert get_functional_unit_gene_families(fu1, simple_gfs, gf2fam) == (set(simple_gfs[:6]), {simple_gfs[-1]}) # fu1 = {GF0, GF1, GF2, GF3, GF4, GF5} | {GF9}
    
    fu2 = simple_model.get('fu2') 
    assert get_functional_unit_gene_families(fu2, simple_gfs, gf2fam) == (set(simple_gfs[6:9]), set()) # GF10 excluded as it is not in the pangenome (not in gf2fam)


def test_search_unit_in_cc(simple_model, simple_gfs, simple_pangenome): 
    # simple_pangenome argument needed to assign metadata to GFs
    detected_units = set()
    comb_families = set(simple_gfs[:6]) # a combination that is subset of (mandatory | accessory) families of fu1 AND found in a connected component of the context graph
    fu = simple_model.get('fu1')
    gf2fam = defaultdict(set, {simple_gfs[i]: {Family(name=f'protein{i}', presence='mandatory' if i < 3 else 'accessory')} for i in range(6)})
    fam2source = {f'protein{i}': 'source1' for i in range(6)}
    
    filtered_cc = nx.Graph() # a filtered connected component of the context graph
    fu_families_list = sorted(comb_families, key=lambda gf: gf.ID) # sort to ensure consistent order
    for i in range(len(fu_families_list) - 1):
        if i==3:  # skip to create a disconnected component
            continue
        filtered_cc.add_edge(fu_families_list[i], fu_families_list[i + 1])
    context_gf = GeneFamily(family_id=11, name='GF11') # add context family node (not in the model)
    filtered_cc.add_edge(fu_families_list[0], context_gf) 
    # filtered_cc = GF0 -- GF1 -- GF2 -- GF3 -- GF11
    #               GF4 -- GF5
    
    # returns detected model families combination as frozenset
    assert search_unit_in_cc(filtered_cc, comb_families, fu, 'source1', gf2fam, fam2source, detected_units) == {frozenset(simple_gfs[:4])} 

    # argument detected_units is updated with the complete identified unit, including context family
    unit_families = set(next(iter(detected_units)).families)
    assert unit_families == set(simple_gfs[:4]) | {context_gf}


def test_search_unit_in_combination(simple_gfs, simple_model, simple_pangenome):
    fu = simple_model.get('fu1')
    gf2fam = defaultdict(set, {simple_gfs[i]: {Family(name=f'protein{i}', presence='mandatory' if i < 3 else 'accessory')} for i in range(6)})
    fam2source = {f'protein{i}': 'source1' for i in range(6)}
    families_in_cc = set(simple_gfs[:6]) # (mandatory | accessory) families of fu1 found in a connected component of the context graph
    matrix = pd.DataFrame([[float(i == j) for j in range(6)] for i in range(6)], # matrix of GF to model family associations 
                index=[f'protein{i}' for i in range(6)],                         # used to check if needed families are covered in each combination
                columns=[f'GF{i}' for i in range(6)])
    
    cc = nx.Graph() # a connected component of the context graph
    for i, gf in enumerate(simple_gfs[:6]):
        cc.add_edge(simple_gfs[i], simple_gfs[i + 1])
    # cc = GF0 -- GF1 -- GF2 -- GF3 -- GF4 -- GF5
    combinations_in_cc = list({frozenset(simple_gfs[:4])}) # list of all combinations subset of families_in_cc
    combinations2orgs = defaultdict(set, {frozenset(simple_gfs[:4]): # dict of all comb-to-orgs associations returned by compute_gene_context_graph of ppanggolin
                                           {'org1', 'org2'}})        # used for local filtering
    detected = search_unit_in_combination(cc, families_in_cc, gf2fam, fam2source, fu, 'source1', matrix, 
                                            combinations_in_cc.copy(), combinations2orgs) # a copy of combinations_in_cc is passed since it is modified in-place
    assert set(next(iter(detected)).families) == set(simple_gfs[:4]) 
    
    # Add context family node
    context_gf = GeneFamily(family_id=11, name='GF11') 
    cc.add_edge(simple_gfs[0], context_gf)
    detected = search_unit_in_combination(cc, families_in_cc, gf2fam, fam2source, fu, 'source1', matrix, 
                                    combinations_in_cc.copy(), combinations2orgs)
    assert set(next(iter(detected)).families) == set(simple_gfs[:4]) | {context_gf}


def test_search_unit_in_context(simple_model, simple_gfs, simple_pangenome):
    detected = set()
    fu_families = set(simple_gfs[:6]) # (mandatory | accessory) families of fu1 that are present in the pangenome
    family_lookup = {f.name: f for f in simple_model.families}
    gf2fam = defaultdict(set, {simple_gfs[i]: {family_lookup[f'protein{i}']} for i in range(10)})

    fam2source = {f'protein{i}': 'source1' for i in range(6)}
    matrix = pd.DataFrame([[float(i == j) for j in range(6)] for i in range(6)], # matrix of GF to model family associations 
                            index=[f'protein{i}' for i in range(6)],             # used to check if needed families are covered in each combination
                            columns=[f'GF{i}' for i in range(6)])
    combinations2orgs = defaultdict(set, {frozenset(simple_gfs[:4]): # dict of all comb-to-orgs associations returned by compute_gene_context_graph of ppanggolin
                                           {'org1', 'org2'}})        # used for local filtering
    func_unit = simple_model.get('fu1') 
    context_graph = nx.Graph() # context graph returned by compute_gene_context_graph of ppanggolin
    for i, _ in enumerate(simple_gfs[:-1]):
        if i == 5: # skip to create a disconnected component
            continue  
        context_graph.add_edge(simple_gfs[i], simple_gfs[i + 1]) # genomes = {'org1'}
    # context_graph = GF0 -- GF1 -- GF2 -- GF3 -- GF4 -- GF5
    #                 GF6 -- GF7 -- GF8
    detected = search_unit_in_context(context_graph, fu_families,  gf2fam, fam2source, matrix, func_unit, 'source1', combinations2orgs, local=True)
    assert set(next(iter(detected)).families) == set(simple_gfs[:4]) # detected unit corresponding to the combination in combinations2orgs


def test_search_system_units(simple_model, simple_gfs, simple_pangenome):
    family_lookup = {f.name: f for f in simple_model.families} # mapping protein names to model Family objects
    gf2fam = defaultdict(set, {simple_gfs[i]: {family_lookup[f'protein{i}']} for i in range(10)})
    fam2source = {f'protein{i}': 'source1' for i in range(6)}
    model_families = set(simple_gfs[:6]) # all GFs corresponding to model families
    detected_systems = search_system_units(simple_model, model_families, gf2fam, fam2source, source='source1')
    # This method executes ppanggolin context graph construction -> context_graph = GF0 -- GF1 -- GF2 -- GF3 -- GF4 -- GF5 -- GF6 -- GF7 -- GF8
    # since all GFs correspond to sequential genes of the same contig of the same organism as defined in simple_gfs fixture
    
    assert detected_systems.keys() == {'fu1'} # only fu1 should be detected
    assert set(next(iter(detected_systems['fu1'])).families) == set(simple_gfs[:7]) # GF6 detected as context family (since transitivity = window = 1)
    assert set(next(iter(detected_systems['fu1'])).models_families) == set(simple_gfs[:6]) # model families of fu1
    
    # Increment window size of fu1
    simple_model.get('fu1').window = 2
    detected_systems = search_system_units(simple_model, model_families, gf2fam, fam2source, source='source1')
    assert set(next(iter(detected_systems['fu1'])).families) == set(simple_gfs[:8]) # GF7 is now detected as context family
    assert set(next(iter(detected_systems['fu1'])).models_families) == set(simple_gfs[:6]) # model families of fu1 remain the same


def test_search_system(simple_model, simple_gfs, simple_pangenome):
    model_fams = set(simple_gfs[:6]) # all GFs corresponding to model families (fu1 only)
    gf2fam = defaultdict(set, {simple_gfs[i]: 
                        {[f for f in simple_model.families if f.name == f'protein{i}'][0]} 
                        for i in range(6)}) # mapping from GFs to model families
    fam2source = {f'protein{i}': 'source1' for i in range(6)}
    
    detected_system = next(iter(search_system(simple_model, model_fams, gf2fam, fam2source, 'source1'))) # extract the first detected system
    expected_unit = SystemUnit(functional_unit=simple_model.get('fu1'), source='source1', 
                               gene_families=set(simple_gfs[:7]), # GF6 does not have corresponding metainfo => context family
                               families_to_metainfo={gf: ('source1', 1) for gf in simple_gfs[:6]})
    
    # the first detected unit of the detected system corresponds to the expected unit (__eq__ based on models_families attribute only)
    assert next(iter((set(detected_system.units)))) == expected_unit
    assert set(expected_unit.models_families) == set(simple_gfs[:6]) # ensure models_families is as expected
    assert set(next(iter((set(detected_system.units)))).families) == set(expected_unit.families) == set(simple_gfs[:7]) # unit families include GF6 as well


def test_search_systems(simple_model, simple_gfs, simple_pangenome):
    # Executes all the previous methods, and all the methods in utils, given the model and pangenome as arguments
    models = Models(models = {simple_model.name: simple_model})
    
    # No return value; systems added to the pangenome
    assert search_systems(models, simple_pangenome, 'source1' , metadata_sources = ['source1'], sensitivity=3) is None
    expected_units = {SystemUnit(functional_unit=simple_model.get('fu1'), source='source1',
                          gene_families=set(simple_gfs[:7]), families_to_metainfo={gf: ('source1', 1) for gf in simple_gfs[:6]})}
    assert set(next(simple_pangenome.systems).units) == expected_units


def test_check_for_needed_units(simple_model, simple_gfs):
    gf2meta_info = {gf: ('source1', 1) for gf in simple_gfs[:6]} # dict of GFs-to-metadata for the GFs of one detected unit
    su = SystemUnit(functional_unit=simple_model.get('fu1'), source='source1', # a detected system unit
                    gene_families=set(simple_gfs[:6]), families_to_metainfo=gf2meta_info)
    detected_units = {'fu1': su}  # dict of all detected units corresponding to each functional unit in the model
    
    # Check if the detected units satisfy the model requirements
    # simple_model has fu1 mandatory and fu2 accessory => satisfied
    assert check_for_needed_units(detected_units, simple_model) is True
    
    # Add fu2 as mandatory functional unit -> model no longer satisfied
    simple_model.mandatory.add(simple_model.get('fu2'))
    simple_model.min_mandatory = 2
    assert check_for_needed_units(detected_units, simple_model) is False
    

def test_get_system_unit_combinations(simple_model, simple_gfs):
    gf2meta_info = {gf: ('source1', 1) for gf in simple_gfs[:6]} # dict of GFs-to-metadata for the GFs of one detected unit
    su1 = SystemUnit(functional_unit=simple_model.get('fu1'), source='source1',
                gene_families=set(simple_gfs[:6]), families_to_metainfo=gf2meta_info)
    detected_units = {'fu1': {su1}} # dict of all detected units corresponding to each functional unit in the model
    
    # Return all system unit combinations satisying the model requirements
    assert get_system_unit_combinations(detected_units, simple_model) == [[su1]]
    
    # Let fu2 be detected 
    gfs = simple_gfs[6:9] + [GeneFamily(family_id=10,name='GF10')] # GF10 is the mandatory family of fu2
    gf2meta_info = {gf: ('source1', 1) for gf in gfs}
    su2 = SystemUnit(functional_unit=simple_model.get('fu2'), source='source1',
                     gene_families=set(simple_gfs), families_to_metainfo=gf2meta_info)
    detected_units = {'fu1': {su1}, 'fu2': {su2}} 
    
    # fu2 is accessory -> model satisfied with and without su2
    assert [set(c) for c in get_system_unit_combinations(detected_units, simple_model)] == [{su1},{su1, su2}] # convert combs to sets to avoid order discrepancies
    
    # Let fu2 be mandatory => model only satisfied with both su and su2
    simple_model.mandatory.add(simple_model.get('fu2'))
    simple_model.accessory.remove(simple_model.get('fu2'))
    simple_model.min_mandatory = 2
    assert [set(c) for c in get_system_unit_combinations(detected_units, simple_model)] == [{su1, su2}] 


def test_search_for_system(): # broken logic; to be fixed
    pass
