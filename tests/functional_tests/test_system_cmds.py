from pathlib import Path
from time import sleep

import pytest

from tests.functional_tests.test_utils_cmd import (  # noqa: F401
    utils_hmm_list,
    utils_model_list,
)
from tests.utils.file_compare import assert_or_update_file
from tests.utils.run_command import run_command


@pytest.fixture(scope="session")
def annotation_and_systems_cmds(
    pangenome_list_file,
    utils_hmm_list,  # noqa: F811
    utils_model_list,  # noqa: F811
    num_cpus,
):
    annotation_command = (
        f"panorama annotation "
        f"--pangenomes {pangenome_list_file} "
        "--source defensefinder "
        f"--hmm {utils_hmm_list} "  # "--mode sensitive " this produce an error
        "--k_best_hit 3 "
        f"--threads {num_cpus}"
    )
    run_command(annotation_command)

    sleep(10)  # Give some time between the two commands

    systems_command = (
        "panorama systems "
        f"--pangenomes {pangenome_list_file} "
        f"--models {utils_model_list} "
        "--source defensefinder "
        "--annotation_sources defensefinder "
        f"--threads {num_cpus}"
    )

    run_command(systems_command)


@pytest.mark.requires_test_data
def test_write_systems(
    annotation_and_systems_cmds,
    pangenome_list_file,
    utils_model_list,  # noqa: F811
    num_cpus,
    update_golden,
):
    outdir = pangenome_list_file.parent / "write_systems_outdir"

    command = (
        f"panorama write_systems  --pangenomes {pangenome_list_file} "
        f"--output {outdir} --models {utils_model_list} "
        "--sources defensefinder "
        "--projection "
        f"--threads {num_cpus} "
        "--association all"
    )

    run_command(command)

    # Validate output structure and files
    assert outdir.exists(), f"Output directory {outdir} was not created"

    # Get list of pangenome names from the pangenome list file
    pangenome_names = []
    with open(pangenome_list_file, "r") as f:
        for line in f:
            if line.strip():
                pangenome_name = line.split("\t")[0]
                pangenome_names.append(pangenome_name)

    # Expected files in each source directory (without projection subdirectory)
    expected_files = [
        "association.tsv",
        "correlation_modules.html",
        "correlation_RGPs.html",
        "correlation_spots.html",
        "module_to_systems.tsv",
        "rgp_to_systems.tsv",
        "spot_to_systems.tsv",
        "systems.tsv",
    ]

    # Check each pangenome directory
    for pangenome_name in pangenome_names:
        pangenome_dir = outdir / pangenome_name
        assert pangenome_dir.exists(), f"Pangenome directory {pangenome_dir} was not created"

        # Check source directory
        source_dir = pangenome_dir / "defensefinder"
        assert source_dir.exists(), f"Source directory {source_dir} was not created"

        # Check all expected files exist
        for expected_file in expected_files:
            file_path = source_dir / expected_file
            assert file_path.exists(), f"Expected file {file_path} was not created"
            assert file_path.stat().st_size > 0, f"Expected file {file_path} is empty"
            if expected_file == "systems.tsv":
                golden_file_name = Path(f"test_write_systems.{pangenome_name}.systems.tsv")
                assert_or_update_file(golden_file_name, file_path, update_golden)

        # Check the projection directory exists (since --projection flag is used)
        projection_dir = source_dir / "projection"
        assert projection_dir.exists(), f"Projection directory {projection_dir} was not created"

        # Check that the projection directory contains TSV files
        projection_files = list(projection_dir.glob("*.tsv"))
        assert len(projection_files) > 0, f"No projection TSV files found in {projection_dir}"

        # Verify all projection files are non-empty
        for proj_file in projection_files:
            assert proj_file.stat().st_size > 0, f"Projection file {proj_file} is empty"
