# default libraries
import hashlib
import logging
import os
import shutil
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Union

# install libraries
import pandas as pd
import pytest
from ppanggolin.genome import Contig, Gene, Organism
from ppanggolin.meta.meta import assign_metadata

# local libraries
from panorama.geneFamily import GeneFamily
from panorama.systems.models import Family, FuncUnit, Model

logger = logging.getLogger(__name__)


def calculate_md5(file_path: Path) -> str:
    """
    Calculate MD5 checksum of a file.

    Args:
        file_path: Path to the file

    Returns:
        MD5 checksum as hexadecimal string
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        # Read file in chunks to handle large files efficiently
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def copy_and_verify(source: Path, destination: Path) -> None:
    """
    Copy a file and verify integrity using MD5 checksum.

    Args:
        source: Source file path
        destination: Destination file path

    Raises:
        ValueError: If MD5 checksums don't match after copying
    """
    # Calculate MD5 of source file
    source_md5 = calculate_md5(source)
    logger.debug(f"Source file {source.name} MD5: {source_md5}")

    # Copy the file
    shutil.copy2(source, destination)

    # Calculate MD5 of copied file
    dest_md5 = calculate_md5(destination)
    logger.debug(f"Copied file {destination.name} MD5: {dest_md5}")

    # Verify checksums match
    if source_md5 != dest_md5:
        raise ValueError(f"MD5 checksum mismatch for {source.name}! Source: {source_md5}, Copied: {dest_md5}")

    logger.debug(f"✅ File {source.name} copied and verified successfully")


def validate_test_data_path(path_str: str = None) -> Union[Path, None]:
    """
    Validate and return the test data path.

    Args:
        path_str: Optional path string. If None, will check environment variable.

    Returns:
        Path object if valid, None otherwise

    Issues warnings for missing or invalid paths.
    """
    # If not provided, check environment variable
    if path_str is None:
        path_str = os.environ.get("PANORAMA_TEST_DATA_PATH")

    # If still not provided, issue a warning
    if path_str is None:
        warnings.warn(
            "Test data path not provided. "
            "Functional tests requiring datasets will be skipped. "
            "Clone https://github.com/labgem/PANORAMA_test and set via --test-data-path"
            " argument or PANORAMA_TEST_DATA_PATH environment variable.",
            UserWarning,
            stacklevel=3,
        )
        return None

    # Convert to Path object, expand user (~) and validate that it exists
    path = Path(path_str).expanduser().resolve()
    if not path.exists():
        warnings.warn(
            f"Test data path '{path}' does not exist. Functional tests requiring datasets will be skipped.",
            UserWarning,
            stacklevel=3,
        )
        return None

    return path


def pytest_addoption(parser):
    parser.addoption(
        "--cpu",
        action="store",
        default="1",
        help="Number of CPUs to use in functional tests",
    )
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Update golden hashes JSON instead of just testing.",
    )
    parser.addoption(
        "--test-data-path",
        action="store",
        default=None,
        help="Path to test dataset repository. Can also be set via "
        "PANORAMA_TEST_DATA_PATH environment variable. "
        "To get test data: git clone https://github.com/labgem/PANORAMA_test",
    )


@pytest.fixture(scope="session")
def num_cpus(request):
    return request.config.getoption("--cpu")


@pytest.fixture
def update_golden(request):
    return request.config.getoption("--update-golden")


@pytest.fixture(scope="session")
def golden_files_path():
    """
    Fixture to provide the path to golden files directory.
    Golden files are stored alongside the test code for version control.
    """
    return Path(__file__).parent / "golden"


@pytest.fixture(scope="session")
def test_data_path(request):
    """
    Fixture to provide the path to test datasets.

    Path can be provided via:
    1. --test-data-path command line argument
    2. PANORAMA_TEST_DATA_PATH environment variable
    3. If neither is provided, returns None and shows a warning
    """
    # First check command line argument
    path_str = request.config.getoption("--test-data-path")

    return validate_test_data_path(path_str)


@pytest.fixture(scope="session")
def pangenome_list_file(test_data_path, tmp_path_factory):
    """
    Create a temporary pangenome list file for testing.
    Copies pangenome files to temporary directory to avoid modifying originals.
    Uses tmp_path_factory for session scope compatibility.
    """
    pangenome_dir = test_data_path / "pangenomes"
    tmp_path = tmp_path_factory.mktemp("panorama_test")

    # Create a subdirectory for copied pangenome files
    tmp_pangenome_dir = tmp_path / "pangenomes"
    tmp_pangenome_dir.mkdir()

    pangenome_list_tsv = tmp_path / "pangenomes_list.tsv"

    with open(pangenome_list_tsv, "w") as f:
        for pangenome_file in pangenome_dir.glob("*.h5"):
            # Copy pangenome file to temporary directory with MD5 verification
            tmp_pangenome_file = tmp_pangenome_dir / pangenome_file.name

            logger.info(f"Copying and verifying pangenome file: {pangenome_file.name}")
            copy_and_verify(pangenome_file, tmp_pangenome_file)

            pangenome_name = pangenome_file.name.rsplit(".h5", 1)[0]
            # Write the path to the copied file in the list
            f.write(f"{pangenome_name}\t{tmp_pangenome_file}\n")

    logger.info(f"All pangenome files copied and verified. List created: {pangenome_list_tsv}")
    return pangenome_list_tsv


def pytest_collection_modifyitems(config, items):
    """
    Handle test collection:
    skip functional tests when no test data is available and reorder tests.
    """

    def get_test_priority(test_function):
        """Determines the priority of a test based on the test file name."""
        # Get the test_name from the test item
        test_name = str(test_function.fspath.relto(test_function.session.fspath))
        # Return the index if file is in our order list, otherwise put it at the end
        try:
            return test_order.index(test_name)
        except ValueError:
            # Files not in the list go to the end, maintain their relative order
            return len(test_order)

    # Get test data path from command line or environment
    test_data_path_str = config.getoption("--test-data-path")
    test_data_path_obj = validate_test_data_path(test_data_path_str)

    # Skip tests that require test data if no valid test data path is available
    if test_data_path_obj is None:
        logger.warning("No valid test data path available. Functional tests will be skipped.")
        skip_functional = pytest.mark.skip(
            reason="Test data not available."
            " Clone https://github.com/labgem/PANORAMA_test and set "
            "--test-data-path or PANORAMA_TEST_DATA_PATH environment variable."
        )

        for item in items:
            # Skip tests that specifically require test data
            if "requires_test_data" in item.keywords:
                item.add_marker(skip_functional)
    else:
        logger.info(f"Using test data path: '{test_data_path_obj}'")

    # Reorder tests: utils first, then detection, then others
    unit_test_system_order = [
        f"unit_tests/systems/{test}"
        for test in [
            "test_systems.py",
            "test_model.py",
            "test_utils.py",
            "test_detection.py",
            "test_systems_projection.py",
        ]
    ]
    unit_test_order = unit_test_system_order
    functional_test_order = [
        f"functional_tests/{test}"
        for test in [
            "test_info.py",
            "test_utils_cmd.py",
            "test_system_cmds.py",
            "test_pansystems.py",
            "test_compare_spots.py",
            "test_compare_systems.py",
        ]
    ]
    test_order = [f"tests/{test}" for test in unit_test_order + functional_test_order]
    items.sort(key=get_test_priority)


# Fixtures used across tests


@pytest.fixture
def simple_contigs():
    return [Contig(identifier=num, name=f"contig_{num}") for num in range(3)]


@pytest.fixture
def simple_orgs(simple_contigs):
    # create 3 organisms with 1 contig each
    orgs = []
    for i, contig in enumerate(simple_contigs):
        org = Organism(name=f"org_{i}")
        org[contig.name] = contig
        orgs.append(org)
    return orgs


@pytest.fixture
def simple_gfs(simple_contigs, simple_orgs):
    gfs = [GeneFamily(family_id=i, name=f"GF{i}") for i in range(10)]
    # Add 3 genes to each GF; one gene for each GF per contig
    for num, contig in enumerate(simple_contigs):
        for i, gf in enumerate(gfs):
            gene = Gene(gene_id=f"gene_{i}_{num}")
            gf[gene.ID] = gene
            gene.family = gf
            gene.fill_parents(organism=simple_orgs[num], contig=contig)
            gene.fill_annotations(position=i, strand="+", start=(i + 1) * 100, stop=(i + 2) * 100)
            contig.add(gene)
    for gf in gfs:  # add partition to each GF
        gf.partition = "P"
    return gfs


@pytest.fixture
def simple_pangenome(simple_gfs, simple_orgs):
    from panorama.pangenomes import Pangenome  # Import here to avoid circular import issue; TODO fix it

    pangenome = Pangenome(name="test_pangenome")
    for gf in simple_gfs:
        pangenome.add_gene_family(gf)
    for org in simple_orgs:
        pangenome.add_organism(org)
    metadata = pd.DataFrame(
        {
            "families": [gf.name for gf in simple_gfs],
            "protein_name": [f"protein{i}" for i in range(10)],
            "score": [1.0 for _ in range(10)],
        }
    )
    assign_metadata(metadata, pangenome, source="source1", metatype="families")
    return pangenome


@pytest.fixture
def simple_fu():
    mandatory_gfs = {Family(name=f"protein{i}", presence="mandatory") for i in range(3)}
    accessory_gfs = {Family(name=f"protein{i + 3}", presence="accessory") for i in range(3)}
    fu = FuncUnit(
        name="fu",
        mandatory=mandatory_gfs,
        accessory=accessory_gfs,
        min_mandatory=2,
        min_total=4,
        transitivity=1,
        window=2,
    )
    return fu


@pytest.fixture
def single_unit_model(simple_fu):
    model = Model(name="TestModel", mandatory={simple_fu})
    simple_fu.model = model
    return model


@pytest.fixture()
def multi_unit_model():
    # Functional Unit 1
    mandatory_gfs = {Family(name=f"protein{i}", presence="mandatory") for i in range(3)}  # fu1_mandatory: GF0, GF1, GF2
    accessory_gfs = {
        Family(name=f"protein{i + 3}", presence="accessory") for i in range(3)
    }  # fu1_accessory: GF3, GF4, GF5
    neutral_gfs = {Family(name="protein9", presence="neutral")}  # fu1_neutral: GF9
    fu1 = FuncUnit(
        name="fu1",
        mandatory=mandatory_gfs,
        accessory=accessory_gfs,
        neutral=neutral_gfs,
        min_mandatory=2,
        min_total=4,
        transitivity=1,
    )
    # Functional Unit 2
    extra_gfs = {Family(name=f"protein{i + 6}", presence="accessory") for i in range(3)}  # fu2_accessory: GF6, GF7, GF8
    fu2 = FuncUnit(
        name="fu2",
        mandatory={Family(name="protein10", presence="mandatory")},
        accessory=extra_gfs,  # fu2_mandatory: GF10
        forbidden={Family(name="protein11", presence="forbidden")},
        min_mandatory=1,
        min_total=3,
        transitivity=1,
    )  # fu2_forbidden: GF11
    # Model: fu1 mandatory; fu2 accessory
    model = Model(name="TestModel", mandatory={fu1}, accessory={fu2})
    fu1.model = fu2.model = model
    return model


@pytest.fixture
def simple_gf2fam(simple_gfs, multi_unit_model):
    family_lookup = {f.name: f for f in multi_unit_model.families}
    gf2fam = defaultdict(set, {simple_gfs[i]: {family_lookup[f"protein{i}"]} for i in range(10)})
    return gf2fam


@pytest.fixture
def simple_fam2source():
    return {f"protein{i}": "source1" for i in range(10)}


@pytest.fixture
def simple_matrix():
    # Creates a matrix of GF-to-model family associations
    matrix = pd.DataFrame(
        [[float(i == j) for j in range(6)] for i in range(6)],
        index=[f"protein{i}" for i in range(6)],
        columns=[f"GF{i}" for i in range(6)],
    )
    return matrix


# Helper classes


class DummyGeneFamily:
    """A dummy gene family class to test context filtering."""

    def __init__(self, name, organisms):
        self.name = name
        self.organisms = organisms
