from pathlib import Path

import pytest

from tests.utils.run_command import run_command


@pytest.mark.requires_test_data
def test_panorama_info(pangenome_list_file: Path, tmp_path: Path):
    output_dir = tmp_path / "output"

    command = (
        f"panorama info -p {pangenome_list_file} -o {output_dir} "
        "--status --content --disable_prog_bar --verbose 2"
    )
    # command = "panorama info -h"
    run_command(command)

    expected_files = [output_dir / "content_info.html", output_dir / "status_info.html"]
    for file in expected_files:
        assert file.exists(), f"Expected output file {file} not found."


@pytest.mark.requires_test_data
def test_panorama_info_content_only(pangenome_list_file: Path, tmp_path: Path):
    output_dir = tmp_path / "output"

    command = (
        f"panorama info -p {pangenome_list_file} -o {output_dir} "
        "--content --disable_prog_bar --verbose 2"
    )
    # command = "panorama info -h"
    run_command(command)

    expected_files = [output_dir / "content_info.html"]
    for file in expected_files:
        assert file.exists(), f"Expected output file {file} not found."
