# Testing Guidelines for PANORAMA 🧪

Welcome to the PANORAMA testing guide! Whether you're adding your first test or refactoring an entire test suite, this
document will help you write tests that are reliable, maintainable, and actually useful. Good tests make everyone's life
easier - including yours!

```{note}
**New to testing or pytest?** No worries! We'll walk you through everything. Testing might seem daunting at first, but once you get the hang of it, you'll wonder how you ever coded without it. And if you get stuck, just ask - we're here to help!
```

## What We're Testing 🔭

PANORAMA uses [pytest](https://docs.pytest.org/en/stable/) as its testing framework. We organize our tests into two main
categories:

- **Unit tests** - Test individual functions, methods, and classes in isolation. Think of these as testing individual
  LEGO bricks.
- **Functional tests** - Test complete workflows and command-line interfaces. These test how the LEGO bricks work
  together to build something useful.

Both types are important! Unit tests catch bugs early and make debugging easier. Functional tests ensure that real-world
usage actually works.

## General Testing Philosophy 📝

### Framework and Structure

Here's how we organize tests in PANORAMA:

- **Use pytest exclusively** - It's powerful, flexible, and widely used
- **Prefer class-based organization** - Groups related tests together nicely
- **Write descriptive test names** - `test_merge_handles_empty_systems` beats `test_merge_1`
- **Use fixtures liberally** - Don't repeat yourself when setting up test data
- **Parametrize similar tests** - Test multiple cases without copy-pasting code

### Basic Test Structure

Here's what a well-organized test class looks like:

```python
class TestSystemMerging:
    """Test the System.merge() functionality."""

    @pytest.fixture
    def empty_system(self):
        """Create an empty system for testing."""
        return System(model=test_model, source="test")

    def test_merge_empty_systems(self, empty_system):
        """Merging two empty systems should work without errors."""
        other_system = System(model=test_model, source="test")
        empty_system.merge(other_system)
        assert len(empty_system) == 0

    @pytest.mark.parametrize("num_units", [1, 5, 10])
    def test_merge_multiple_units(self, empty_system, num_units):
        """Test merging systems with varying numbers of units."""
        # Test implementation
        pass
```

Clean, organized, and easy to understand!

## Unit Testing 🔬

Unit tests are your first line of defense against bugs. They test individual pieces of code in isolation.

### What Makes a Good Unit Test?

1. **Isolation** - Each test stands alone and doesn't depend on external systems or other tests
2. **Speed** - Unit tests should be fast (ideally under a second)
3. **Focused** - Test one thing at a time
4. **Reliable** - Same input always produces same output

### Best Practices for unit tests

**Test both success and failure:**

```python
def test_valid_gene_family(self):
    """Test that valid gene families are accepted."""
    gf = GeneFamily(name="valid_family", family_id=1)
    assert gf.name == "valid_family"


def test_invalid_gene_family_raises_error(self):
    """Test that invalid input raises appropriate error."""
    with pytest.raises(ValueError, match="Invalid family ID"):
        GeneFamily(name="test", family_id=-1)
```

**Test edge cases:**

- Empty inputs
- Very large inputs
- Boundary values
- None/null values
- Duplicate entries

**Use descriptive names:**

```python
# Good
def test_merge_fails_when_models_differ(self):
    pass


# Less helpful
def test_merge_2(self):
    pass
```

**Group related tests in classes:**

```python
class TestSystemUnit:
    """All tests related to SystemUnit functionality."""

    @pytest.fixture
    def basic_unit(self):
        """Shared fixture for the class."""
        return SystemUnit(...)

    def test_creation(self, basic_unit):
        pass

    def test_addition(self, basic_unit):
        pass
```

## Functional Testing 🚀

Functional tests verify that complete features work as users would actually use them.

### What Makes a Good Functional Test?

1. **Realistic** - Test real workflows with realistic data
2. **End-to-end** - Test the full pipeline, not just pieces
3. **User-focused** - Test what users actually do
4. **Thorough** - Verify outputs, not just that commands don't crash

### Best Practices for functional tests

**Use session-scoped fixtures for expensive operations:**

```python
@pytest.fixture(scope="session")
def test_pangenome():
    """Create test pangenome once for all tests."""
    # This might take a while, so we only do it once
    return create_test_pangenome()
```

**Mark tests that need external data:**

```python
@pytest.mark.requires_test_data
def test_annotation_pipeline(test_data_path):
    """Test the annotation pipeline with real data."""
    command = f"panorama annotate --pangenome {test_data_path}/test.h5"
    run_command(command)
```

**Test command-line interfaces exactly as users invoke them:**

```python
def test_systems_command():
    """Test the systems command with typical user parameters."""
    command = (
        f"panorama systems "
        f"--pangenomes {pangenome_list} "
        f"--models {model_file} "
        f"--source defensefinder"
    )
    result = run_command(command)
    assert result.returncode == 0
```

## Setting Up Your Test Environment 🛠️

### Install Test Dependencies

First, install the testing tools:

```shell
# From the PANORAMA root directory
pip install -e .[test]  # -e for editable mode (recommended)
```

This installs pytest, coverage tools, and everything else you need.

### Get the Test Data

We keep test data in a separate repository to keep the main repo lightweight:

```shell
# Clone the test data repository
git clone https://github.com/labgem/PANORAMA_test

# Set the environment variable (or use --test-data-path flag)
export PANORAMA_TEST_DATA_PATH=PANORAMA_test/
```

Now you're ready to run tests!

## Running Tests ▶️

### Local Testing

Here are the most common commands you'll use:

```bash
# Run all tests
pytest

# Run with coverage (shows which code is tested)
pytest --cov=panorama

# Run specific test file
pytest tests/test_systems.py

# Run specific test class
pytest tests/test_systems.py::TestSystemUnit

# Run specific test method
pytest tests/test_systems.py::TestSystemUnit::test_init_basic

# Run only tests marked as requiring test data
pytest -m "requires_test_data"

# Run tests in parallel (faster!)
pytest -n auto

# Verbose output (helpful for debugging)
pytest -v

# Show print statements (useful for debugging)
pytest -s
```

**If you haven't set the environment variable:**

```bash
pytest --test-data-path=/path/to/PANORAMA_test
```

### What the CI Does

When you push code, GitHub Actions automatically:

1. Downloads the test data from the separate repository
2. Installs all dependencies
3. Runs the full test suite with coverage
4. Reports results back to your pull request

The workflow looks like this:

```bash
# Download test data
git clone https://github.com/labgem/PANORAMA_test test_data

# Run tests with coverage
pytest --cov=panorama --test-data-path=test_data
```

If tests fail in CI but pass locally, it might be an environment issue. Don't hesitate to ask for help debugging!

## Creating Reusable Test Components 🔧

### Fixtures Are Your Friends

Fixtures let you reuse test setup code without repeating yourself:

```python
class TestFixture:
    """Base class for shared fixtures."""

    @pytest.fixture
    def model(self):
        """Create a test model."""
        return Model(
            name="test_model",
            min_mandatory=1,
            min_total=1,
            canonical=["canonical_1", "canonical_2"],
        )

    @pytest.fixture
    def functional_unit(self, model):
        """Create a test functional unit (depends on model fixture)."""
        fu = FuncUnit(name="test_unit", presence="mandatory", min_total=2)
        fu.model = model
        return fu
```

Now any test can use these fixtures by including them as parameters!

### Helper Methods

For complex setup that's specific to a test class, create helper methods:

```python
class TestGeneFamily:
    def create_gene_family(self, name, num_organisms=5):
        """Helper to create a gene family with organisms."""
        gf = GeneFamily(name=name, family_id=next_id())
        for i in range(num_organisms):
            org = Organism(name=f"org_{i}")
            gf.add_organism(org)
        return gf

    def test_family_with_many_organisms(self):
        """Test gene family with many organisms."""
        gf = self.create_gene_family("test", num_organisms=100)
        assert len(gf.organisms) == 100
```

## Testing Errors and Edge Cases ⚠️

### Always Test Error Conditions

Good code handles errors gracefully. Test that it does!

```python
def test_division_by_zero_raises_error(self):
    """Test that division by zero is handled properly."""
    with pytest.raises(ZeroDivisionError, match="cannot divide by zero"):
        calculator.divide(10, 0)
```

### Parametrize Error Tests

Test multiple invalid inputs efficiently:

```python
@pytest.mark.parametrize("invalid_input", [
    "not_a_number",
    None,
    [],
    {},
    -1,
    float('inf'),
])
def test_rejects_invalid_input(self, invalid_input):
    """Test that various invalid inputs are rejected."""
    with pytest.raises(TypeError):
        process_data(invalid_input)
```

This is much cleaner than writing six separate test methods!

## Test Quality and Coverage 📈

### Coverage Goals

We aim for high coverage, but smart coverage:

- **Target >90% line coverage** - This catches most bugs
- **Don't chase 100%** - Diminishing returns and sometimes impossible
- **Focus on critical paths** - Test important code thoroughly
- **Use coverage to find gaps** - Not as a grade to optimize

```bash
# Generate coverage report
pytest --cov=panorama --cov-report=html

# Open htmlcov/index.html in your browser to see details
```

Red lines in the coverage report? Those might need tests!

### What Makes Tests High Quality?

Good tests are:

- **Fast** - Unit tests run in milliseconds, functional tests in seconds
- **Reliable** - Same result every time, no flaky tests
- **Independent** - Can run in any order
- **Readable** - Someone else can understand what's being tested
- **Maintainable** - Easy to update when code changes

If a test is slow, flaky, or confusing, it needs improvement!

## Working with CI/CD 🔄

### How Tests Integrate with Development

Here's the typical flow:

1. You write code and tests locally
2. Run tests locally (`pytest`)
3. Push to your branch
4. GitHub Actions runs the full test suite
5. Tests pass → ready for review!
6. Tests fail → fix and push again

**Pro tip:** Run tests locally before pushing. It's much faster than waiting for CI!

### Understanding CI Failures

If tests fail in CI:

- Check the logs carefully
- Try to reproduce locally
- Could be a test data issue
- Could be an environment difference
- Ask for help if stuck!

## Troubleshooting Common Issues 🔍

### Fixtures Not Working

**Problem:** Fixture isn't available in my test

**Solution:** Check the fixture scope and location

```python
# Fixture in same class or conftest.py
@pytest.fixture
def my_fixture(self):  # or without 'self' in conftest.py
    return test_data
```

### Tests Pass Locally, Fail in CI

**Common causes:**

- File paths (use `Path` from `pathlib`, not hardcoded strings)
- Test data location (use `--test-data-path` or environment variable)
- Timing issues (avoid tests that depend on exact timing)
- Hidden dependencies between tests

### Slow Tests

**Solutions:**

- Use `session` scope for expensive fixtures
- Mock external calls
- Use smaller test datasets
- Run tests in parallel: `pytest -n auto`

### Test Data Issues

**Best practices:**

- Keep test data minimal but representative
- Use factories/builders to create test objects
- Don't hardcode paths - use fixtures or parameters
- Version control small test files, external repo for large ones

## Getting Help with Testing 🆘

Stuck on something? Here's what to do:

1. **Check the pytest documentation** - [docs.pytest.org](https://docs.pytest.org/)
2. **Look at existing tests** - Find similar tests in the codebase
3. **Ask in discussions** - Other developers have probably hit the same issue
4. **Open a draft PR** - Get feedback on your tests
5. **Tag a maintainer** - We're happy to help!

**Remember:** Writing good tests is a skill that takes practice. Your first tests might not be perfect, and that's
completely fine. We're here to help you improve!

## For New Contributors 🌍

We genuinely appreciate when contributors include tests with their code! Here's what you should know:

**Test expectations for PRs:**

- New features should include tests
- Bug fixes should include a regression test
- Don't worry about perfect coverage on your first try
- We'll help you improve tests during review

**Not sure how to test something?**

- Ask! Testing advice is part of code review
- Check similar existing tests for patterns
- Start simple - basic tests are better than no tests
- We can help you expand test coverage

**Common questions:**

- *"Do I need to test this tiny change?"* - If it's code, yes! Even simple tests catch bugs.
- *"How much should I test?"* - Test happy paths, error cases, and edge cases.
- *"My test is complicated, is that OK?"* - Maybe! We can help simplify it.

Testing makes PANORAMA better for everyone. Thanks for taking the time to test your code! 🎉

## Quick Reference 📋

```bash
# Setup (one time)
pip install -e .[test]
git clone https://github.com/labgem/PANORAMA_test

# Common test commands
pytest                                    # Run all tests
pytest --cov=panorama                    # Run with coverage
pytest tests/test_systems.py             # Run specific file
pytest -v                                # Verbose output
pytest -s                                # Show print statements
pytest -n auto                           # Parallel execution
pytest -k "test_merge"                   # Run tests matching pattern
pytest -m "requires_test_data"           # Run marked tests

# Before pushing
pytest --cov=panorama                    # Make sure tests pass!
```

Happy testing! 🧪