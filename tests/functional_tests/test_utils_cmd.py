import pytest
from tests.utils.run_command import run_command
from pathlib import Path


@pytest.fixture(scope="session")
def utils_hmm_list(test_data_path: Path, tmp_path_factory):
    hmm_dir = test_data_path / "hmms"
    outdir = tmp_path_factory.mktemp("utils") / "hmm_utils" # Create temp dir using factory
    
    command = f"panorama utils --hmm {hmm_dir} -o {outdir}"

    run_command(command)

    return outdir / "hmm_list.tsv"

@pytest.mark.requires_test_data
def test_hmm_list(utils_hmm_list: Path):

    assert utils_hmm_list.exists(), f"Expected output file {utils_hmm_list} not found."

    with open(utils_hmm_list, "r") as f:
        lines = f.readlines()
        assert len(lines) > 1, "HMM list file is empty or has no entries."

        for line in lines[1:]:
            path = line.strip().split("\t")[2]
            assert Path(path).exists(), f"HMM file {path} in {utils_hmm_list} does not exist."

@pytest.fixture(scope="session")
def utils_model_list(test_data_path: Path, tmp_path_factory):
    model_dir = test_data_path / "models"
    outdir = tmp_path_factory.mktemp("utils") / "model_utils" # Create temp dir using factory
    model_files = " ".join((str(f) for f in model_dir.glob("*.json")))
    command = f"panorama utils --models {model_files} -o {outdir}"

    run_command(command)

    return outdir / "models_list.tsv"


@pytest.mark.requires_test_data
def test_model_list(utils_model_list: Path):

    assert utils_model_list.exists(), f"Expected output file {utils_model_list} not found."
    # assert file not empty
    assert utils_model_list.stat().st_size > 0, f"File {utils_model_list} is empty."
