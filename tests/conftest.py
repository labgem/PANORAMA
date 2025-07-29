from collections import defaultdict

import pytest
import pandas as pd
from ppanggolin.meta.meta import assign_metadata
from ppanggolin.genome import Gene, Organism, Contig

from panorama.geneFamily import GeneFamily
from panorama.systems.models import Model, FuncUnit, Family
from panorama.systems.system import System


# Fixtures used across tests

@pytest.fixture
def simple_contigs():
    return [Contig(identifier=num, name=f'contig_{num}') for num in range(3)]

@pytest.fixture
def simple_orgs(simple_contigs):
    # create 3 organisms with 1 contig each
    orgs = []
    for i, contig in enumerate(simple_contigs):
        org = Organism(name=f'org_{i}')
        org[contig.name] = contig
        orgs.append(org)
    return orgs

@pytest.fixture
def simple_gfs(simple_contigs, simple_orgs):
    gfs = [GeneFamily(family_id=i, name=f'GF{i}') for i in range(10)]
    # Add 3 genes to each GF; one gene for each GF per contig
    for num, contig in enumerate(simple_contigs):
        for i, gf in enumerate(gfs):
            gene = Gene(gene_id=f'gene_{i}_{num}')
            gf[gene.ID] = gene
            gene.family = gf
            gene.fill_parents(organism = simple_orgs[num], contig = contig)
            gene.fill_annotations(position=i, strand='+', start=(i+1)*100, stop=(i+2)*100)
            contig.add(gene)
    for gf in gfs: # add partition to each GF
        gf.partition = 'P'
    return gfs

@pytest.fixture
def simple_pangenome(simple_gfs, simple_orgs):
    from panorama.pangenomes import Pangenome # Import here to avoid circular import issue; TODO fix it
    pangenome = Pangenome(name='test_pangenome')
    for gf in simple_gfs:
        pangenome.add_gene_family(gf)
    for org in simple_orgs:
        pangenome.add_organism(org)
    metadata = pd.DataFrame({
        'families': [gf.name for gf in simple_gfs],
        'protein_name': [f'protein{i}' for i in range(10)],
        'score': [1.0 for _ in range(10)],
        })
    assign_metadata(metadata, pangenome, source='source1', metatype='families')
    return pangenome

@pytest.fixture
def simple_fu():
    mandatory_gfs = {Family(name=f'protein{i}', presence='mandatory') for i in range(3)}
    accessory_gfs = {Family(name=f'protein{i+3}', presence='accessory') for i in range(3)}
    fu = FuncUnit(name='fu', mandatory=mandatory_gfs, accessory=accessory_gfs,
                  min_mandatory=2, min_total=4, transitivity=1, window=2)
    return fu

@pytest.fixture
def single_unit_model(simple_fu):
    model = Model(name='TestModel', mandatory={simple_fu})
    simple_fu.model = model
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


