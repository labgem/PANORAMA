from pathlib import Path
import shutil


GOLDEN_DIR = Path("tests/functional_tests/expected_outputs")


def assert_or_update_file(expected_file_name: Path, file_path: Path, update_golden: bool):
    """Compare file content with golden file or update it if --update-golden is set."""
    golden_file = GOLDEN_DIR / expected_file_name

    if update_golden:
        # Copy current file to golden
        shutil.copy(file_path, golden_file)
        print(f"[update-golden] Updated golden file for {file_path.name}")
    else:
        assert (
            golden_file.exists()
        ), f"No golden file for '{expected_file_name}' found in {GOLDEN_DIR}. Run pytest with --update-golden first."
        content_actual = file_path.read_text()
        content_expected = golden_file.read_text()
        print(content_actual)
        assert content_actual == content_expected, (
            f"Content mismatch for {file_path.name}. "
            f"Use --update-golden to update the reference."
        )

