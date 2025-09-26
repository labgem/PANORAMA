# Testing with PANORAMA

## Running Tests

PANORAMA uses pytest for testing with two types of tests:

### Unit Tests
Unit tests run by default and don't require any special setup:
```bash
pytest  # Runs all unit tests
```

### Functional Tests  
Functional tests require test datasets from the PANORAMA_test repository and run automatically when test data is available.

#### How to Run Functional Tests

**Step 1: Clone the test data repository**
```bash
git clone https://github.com/labgem/PANORAMA_test.git
```

**Step 2: Provide the test data path**
You can provide the path to the cloned repository in one of two ways:

**Option A: Command line argument**
```bash
pytest --test-data-path=/path/to/PANORAMA_test
```

**Option B: Environment variable**
```bash
export PANORAMA_TEST_DATA_PATH=/path/to/PANORAMA_test
pytest
```

#### Test Execution Behavior
- **No test data path**: Only unit tests run; functional tests are automatically skipped
- **Invalid test data path**: Only unit tests run; functional tests are skipped if path doesn't exist  
- **Valid test data path**: Both unit and functional tests run

#### What Happens Without Test Data
If you don't provide test data or the path doesn't exist, functional tests are automatically skipped. You'll see output like:

```
tests/test_functional.py::test_something SKIPPED (Test data not available. Clone https://github.com/labgem/PANORAMA_test and set --test-data-path or PANORAMA_TEST_DATA_PATH environment variable.)
```

## Available Pytest Options

- `--test-data-path=PATH`: Path to test dataset repository. When provided, functional tests will be executed automatically.
- `--cpu=N`: Number of CPUs to use in functional tests (default: 1)
- `--update-golden`: Update golden hashes JSON instead of just testing

## Test Markers

Tests can be marked with:
- `@pytest.mark.requires_test_data`: Marks tests that need test datasets

## Example Functional Test

```python
import pytest
import os
import subprocess

@pytest.mark.requires_test_data
def test_panorama_with_dataset(test_data_path, num_cpus):
    """Example functional test using test datasets."""
    if test_data_path is None:
        pytest.skip("Test data path not available")
    
    input_file = os.path.join(test_data_path, "sample.fasta")
    if not os.path.exists(input_file):
        pytest.skip(f"Required test file {input_file} not found")
    
    cmd = ["panorama", "command", "--input", input_file, "--cpu", num_cpus]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
```

## Continuous Integration

In CI environments, set the environment variable:
```yaml
env:
  PANORAMA_TEST_DATA_PATH: /path/to/test/data
```

Or pass as argument:
```bash
pytest --test-data-path=$TEST_DATA_PATH
```
