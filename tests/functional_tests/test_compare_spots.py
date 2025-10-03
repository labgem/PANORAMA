import pytest

from pathlib import Path
from tests.utils.run_command import run_command
from tests.functional_tests.test_utils_cmd import utils_model_list, utils_clustering
from tests.utils.file_compare import assert_or_update_file


@pytest.mark.requires_test_data
def test_compare_spots_command(
    pangenome_list_file, utils_model_list, utils_clustering, num_cpus, update_golden
):
    # TODO Manage the problem with reproducibility for better testing

    outdir = pangenome_list_file.parent / "compare_spots_outdir"
    compare_spots_command = (
        "panorama compare_spots "
        f"--pangenomes {pangenome_list_file} "
        f"-o {outdir} "
        "--graph_formats gexf graphml "
        f"--cpus {num_cpus} "
        f"--cluster {utils_clustering} "
        "--systems "
        f"-s defensefinder "
        f"-m {utils_model_list} "
    )
    run_command(compare_spots_command)

    # Validate output structure and files
    assert outdir.exists(), f"Output directory {outdir} was not created"

    # Expected files in each source directory
    assert_or_update_files = [
        "all_conserved_spots.tsv",
        # "conserved_spots/conserved_spot_3.tsv",
        # "conserved_spots/conserved_spot_24.tsv",
        # "conserved_spots/conserved_spot_33.tsv",
        # "conserved_spots/conserved_spot_60.tsv",
    ]
    expected_files = map(
        Path,
        assert_or_update_files
        + [
            "conserved_spots.gexf",
            "conserved_spots.graphml",
            "systems_link_with_conserved_spots_louvain.gexf",
            "systems_link_with_conserved_spots_louvain.graphml",
        ],
    )
    # Check all expected files exist
    for expected_file in expected_files:
        file_path = outdir / expected_file
        assert file_path.exists(), f"Expected file {file_path} was not created"
        assert file_path.stat().st_size > 0, f"Expected file {file_path} is empty"
        # if expected_file in assert_or_update_files:
        #     assert_or_update_file(expected_file, file_path, update_golden)
