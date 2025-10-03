from pathlib import Path

import pytest

from tests.utils.run_command import run_command
from tests.functional_tests.test_utils_cmd import utils_model_list, utils_clustering


@pytest.mark.requires_test_data
def test_compare_systems_command(
    pangenome_list_file, utils_model_list, utils_clustering, num_cpus, update_golden
):
    # TODO Manage the problem with reproducibility for better testing

    outdir = pangenome_list_file.parent / "compare_systems_outdir"
    compare_systems_command = (
        "panorama compare_systems "
        f"--pangenomes {pangenome_list_file} "
        f"-o {outdir} "
        f"-s defensefinder "
        f"-m {utils_model_list} "
        "--graph_formats gexf graphml "
        "--gfrr_metrics min_gfrr_models "
        "--heatmap "
        f"--cpus {num_cpus} "
        f"--cluster {utils_clustering} "
    )
    run_command(compare_systems_command)

    # Validate output structure and files
    assert outdir.exists(), f"Output directory {outdir} was not created"

    # Expected files in each source directory
    assert_or_update_files = []
    expected_files = map(
        Path,
        assert_or_update_files
        + [
            "heatmap_number_systems.html",
            "heatmap_normalized_systems.html",
            "conserved_systems.gexf",
            "conserved_systems.graphml",
        ],
    )
    # Check all expected files exist
    for expected_file in expected_files:
        file_path = outdir / expected_file
        assert file_path.exists(), f"Expected file {file_path} was not created"
        assert file_path.stat().st_size > 0, f"Expected file {file_path} is empty"
        # if expected_file in assert_or_update_files:
        #     assert_or_update_file(expected_file, file_path, update_golden)
