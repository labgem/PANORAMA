from collections import defaultdict

import pytest
import pandas as pd
from ppanggolin.meta.meta import assign_metadata
from ppanggolin.genome import Gene, Organism, Contig

from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family


# Fixtures used across tests

@pytest.fixture
def simple_org():
    org = Organism(name='test_org')
    return org

@pytest.fixture
def simple_contig():
    contig = Contig(identifier=1, name='test_contig')
    return contig

@pytest.fixture
def simple_gfs(simple_org, simple_contig):
    gfs = [GeneFamily(family_id=i, name=f'GF{i}') for i in range(10)]
    # Add a gene to each GF and fill its attributes
    for i, gf in enumerate(gfs):
        gene = Gene(gene_id=f'gene_{i}')
        gene.family = gf
        gf[gene.ID] = gene
        gene.fill_parents(organism = simple_org, contig = simple_contig)
        gene.fill_annotations(position=i, strand='+', start=(i+1)*100, stop=(i+2)*100)
        simple_contig.add(gene)
    return gfs

@pytest.fixture
def simple_pangenome(simple_gfs):
    from panorama.pangenomes import Pangenome # Import here to avoid circular import issue; TODO fix it
    pangenome = Pangenome(name='test_pangenome')
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

@pytest.fixture
def single_unit_model():
    mandatory_gfs = {Family(name=f'protein{i}', presence='mandatory') for i in range(3)} 
    accessory_gfs = {Family(name=f'protein{i+3}', presence='accessory') for i in range(3)}
    fu = FuncUnit(mandatory=mandatory_gfs, accessory=accessory_gfs, min_mandatory=2, min_total=4,transitivity=1)
    model = Model(name='TestModel', mandatory={fu})
    return model

@pytest.fixture()
def multi_unit_model():
    # Functional Unit 1
    mandatory_gfs = {Family(name=f'protein{i}', presence='mandatory') for i in range(3)} # fu1_mandatory: GF0, GF1, GF2
    accessory_gfs = {Family(name=f'protein{i+3}', presence='accessory') for i in range(3)} # fu1_accessory: GF3, GF4, GF5
    neutral_gfs = {Family(name='protein9', presence='neutral')} # fu1_neutral: GF9
    fu1 = FuncUnit(name='fu1', mandatory=mandatory_gfs, accessory=accessory_gfs, neutral=neutral_gfs, 
                min_mandatory=2, min_total=4, transitivity=1)
    # Functional Unit 2
    extra_gfs = {Family(name=f'protein{i+6}', presence='accessory') for i in range(3)} # fu2_accessory: GF6, GF7, GF8
    fu2 = FuncUnit(name='fu2', mandatory={Family(name='protein10', presence='mandatory')}, accessory=extra_gfs,              # fu2_mandatory: GF10
                forbidden={Family(name='protein11', presence='forbidden')}, min_mandatory=1, min_total=3, transitivity=1) # fu2_forbidden: GF11
    # Model: fu1 mandatory; fu2 accessory
    model = Model(name='TestModel', mandatory={fu1}, accessory={fu2})
    fu1.model = fu2.model = model
    return model

@pytest.fixture
def simple_gf2fam(simple_gfs, multi_unit_model):
    family_lookup = {f.name: f for f in multi_unit_model.families}
    gf2fam = defaultdict(set, {simple_gfs[i]: {family_lookup[f'protein{i}']} for i in range(10)})
    return gf2fam

@pytest.fixture
def simple_fam2source():
    return {f'protein{i}': 'source1' for i in range(10)}

@pytest.fixture
def simple_matrix():
    # Creates a matrix of GF-to-model family associations
    matrix = pd.DataFrame([[float(i == j) for j in range(6)] for i in range(6)],
                        index=[f'protein{i}' for i in range(6)],
                        columns=[f'GF{i}' for i in range(6)])
    return matrix


# Helper classes

class DummyGeneFamily:
    """A dummy gene family class to test context filtering."""
    def __init__(self, name, organisms):
        self.name = name
        self.organisms = organisms


# Specifying the order of test files execution

def pytest_collection_modifyitems(config, items):
    """Reorder tests: utils first, then detection, then others."""
    utils_tests = []
    detection_tests = []
    other_tests = []
    
    for item in items:
        if "test_utils" in str(item.fspath):
            utils_tests.append(item)
        elif "test_detection" in str(item.fspath):
            detection_tests.append(item)
        else:
            other_tests.append(item)
    
    # Reorder: utils first, then detection, then others
    items[:] = utils_tests + detection_tests + other_tests


