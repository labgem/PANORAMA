# Tips and Tricks for PANORAMA Developers

## Table of Contents
1. [Code Quality Fundamentals](#code-quality-fundamentals)
2. [Testing Strategy](#testing-strategy)
3. [Debugging and Profiling](#debugging-and-profiling)
4. [Git Workflow](#git-workflow)
5. [Common Pitfalls](#common-pitfalls)

---

## Code Quality Fundamentals

### Naming Conventions

Good names make code self-documenting:

```python
# Variables and functions: snake_case
gene_family_count = len(gene_families)
def calculate_coverage(genes):
    pass

# Classes: CamelCase
class GeneFamily:
    pass

# Constants: UPPER_SNAKE_CASE
MAX_ITERATIONS = 1000
DEFAULT_THRESHOLD = 0.95

# Private attributes: leading underscore
class System:
    def __init__(self):
        self._internal_cache = {}
```

**Be descriptive but concise:**
```python
# Good
def merge_systems(system1, system2):
    pass

# Too vague
def merge(s1, s2):
    pass

# Too verbose
def merge_two_system_objects_together(first_system, second_system):
    pass
```

### Error Handling

#### Handle Errors Gracefully

Always anticipate what could go wrong and handle it appropriately:

```python
def load_pangenome(filepath: str) -> Pangenome:
    """Load a pangenome from an HDF5 file."""
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Pangenome file not found: {filepath}")
    
    try:
        pangenome = Pangenome.from_file(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to load pangenome: {e}") from e
    
    return pangenome
```

#### Choose the Right Exception Type

Use appropriate exception types:

```python
# Invalid input from user
raise ValueError("threshold must be between 0 and 1")

# File doesn't exist
raise FileNotFoundError(f"Model file not found: {path}")

# Wrong type provided
raise TypeError(f"Expected GeneFamily, got {type(obj)}")

# Feature not implemented yet
raise NotImplementedError("Clustering method X not yet supported")

# Key doesn't exist in dict
raise KeyError(f"System '{system_id}' not found")
```

#### Validate Input Early

Check inputs at the start of functions:

```python
def calculate_similarity(family1: GeneFamily, family2: GeneFamily) -> float:
    """Calculate Jaccard similarity between two gene families."""
    # Validate inputs
    if not isinstance(family1, GeneFamily):
        raise TypeError(f"family1 must be GeneFamily, got {type(family1)}")
    
    if not isinstance(family2, GeneFamily):
        raise TypeError(f"family2 must be GeneFamily, got {type(family2)}")
    
    if len(family1) == 0 or len(family2) == 0:
        raise ValueError("Cannot calculate similarity for empty families")
    
    # Now we know inputs are valid
    intersection = len(family1.genes & family2.genes)
    union = len(family1.genes | family2.genes)
    return intersection / union
```

### Logging

Use Python's logging module for informational messages:

```python
import logging

logger = logging.getLogger(__name__)

def process_pangenome(pangenome):
    """Process a pangenome with logging."""
    logger.info(f"Processing pangenome with {len(pangenome.gene_families)} families")
    
    try:
        result = complex_operation(pangenome)
        logger.debug(f"Complex operation completed: {result}")
    except Exception as e:
        logger.error(f"Failed to process pangenome: {e}")
        raise
    
    logger.info("Processing completed successfully")
    return result
```

**Logging levels:**
- `DEBUG` - Detailed diagnostic information
- `INFO` - General informational messages
- `WARNING` - Something unexpected but not an error
- `ERROR` - Something failed
- `CRITICAL` - Serious failure

### Performance Best Practices

#### Use Generators for Large Datasets

```python
# Good: Memory efficient
def iter_gene_families(pangenome):
    for family in pangenome.families:
        yield family

# Less efficient: Loads everything into memory
def get_all_families(pangenome):
    return [family for family in pangenome.families]
```

#### Cache Expensive Computations

```python
from functools import lru_cache

@lru_cache(maxsize=1000)
def calculate_similarity(family_id1, family_id2):
    """Calculate similarity with caching."""
    # Expensive computation here
    pass
```

#### Use Built-in Functions

They're optimized in C and much faster:

```python
# Fast
total = sum(len(family) for family in families)

# Slower
total = 0
for family in families:
    total += len(family)
```

#### Use Sets for Membership Testing

```python
# Fast: O(1) lookup
gene_ids = set(family.gene_ids)
if gene_id in gene_ids:
    pass

# Slow: O(n) lookup
gene_ids = list(family.gene_ids)
if gene_id in gene_ids:
    pass
```

---

## Testing Strategy

### Unit Testing 🔬

Unit tests are your first line of defense against bugs. They test individual pieces of code in isolation.

#### What Makes a Good Unit Test?

1. **Isolation** - Each test stands alone and doesn't depend on external systems or other tests
2. **Speed** - Unit tests should be fast (ideally under a second)
3. **Focused** - Test one thing at a time
4. **Reliable** - The Same input always produces the same output

#### Test Both Success and Failure

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

#### Test Edge Cases

Always test:
- Empty inputs
- Very large inputs
- Boundary values
- None/null values
- Duplicate entries

#### Use Descriptive Names

```python
# Good
def test_merge_fails_when_models_differ(self):
    pass

# Less helpful
def test_merge_2(self):
    pass
```

#### Group Related Tests

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

### Functional Testing 🚀

Functional tests verify that complete features work as users would actually use them.

#### What Makes a Good Functional Test?

1. **Realistic** - Test real workflows with realistic data
2. **End-to-end** - Test the full pipeline, not just pieces
3. **User-focused** - Test what users actually do
4. **Thorough** - Verify outputs, not just that commands don't crash

#### Use Session-Scoped Fixtures

```python
@pytest.fixture(scope="session")
def test_pangenome():
    """Create test pangenome once for all tests."""
    # This might take a while, so we only do it once
    return create_test_pangenome()
```

#### Mark Tests That Need External Data

```python
@pytest.mark.requires_test_data
def test_annotation_pipeline(test_data_path):
    """Test the annotation pipeline with real data."""
    command = f"panorama annotate --pangenome {test_data_path}/test.h5"
    run_command(command)
```

#### Test Command-Line Interfaces

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

### Creating Reusable Test Components 🔧

#### Fixtures Are Your Friends

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

#### Helper Methods

For complex setup that's specific to a test class:

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

### Testing Errors and Edge Cases ⚠️

#### Always Test Error Conditions

```python
def test_division_by_zero_raises_error(self):
    """Test that division by zero is handled properly."""
    with pytest.raises(ZeroDivisionError, match="cannot divide by zero"):
        calculator.divide(10, 0)
```

#### Parametrize Error Tests

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

---

## Debugging and Profiling

### Using the Python Debugger 🔍

Don't just add print statements - use Python's debugger:

```python
# Add this line where you want to break
import pdb; pdb.set_trace()

# Or in Python 3.7+
breakpoint()
```

Common debugger commands:
- `n` - Next line
- `s` - Step into function
- `c` - Continue execution
- `p variable` - Print variable value
- `l` - List surrounding code
- `q` - Quit debugger

```{hint}
Your editor might integrate a debugger, that way you don't have to type commands manually.
```

### Better Print Debugging

If you must use print statements:

```python
# Basic print
print(f"DEBUG: family_count = {family_count}")

# Pretty print complex objects
from pprint import pprint
pprint(complex_dict)

# Print with context
import sys
print(f"DEBUG [{sys._getframe().f_code.co_name}]: value = {value}")
```

```{caution}
Remember to remove debug prints before committing!
```

### Performance Profiling with VizTracer

VizTracer creates visual timelines showing exactly where your code spends time:

```shell
# Profile a script
viztracer my_script.py --output profile.json

# Profile with specific arguments
viztracer panorama systems --pangenomes data.txt --models models.yml

# Open the visualization
vizviewer profile.json
# This opens a browser showing an interactive timeline
```

**What to look for in the timeline:**
- Functions that take a long time
- Functions are called very frequently
- Unexpected I/O operations
- Nested loops that could be optimized

---

## Git Workflow

### Writing Good Commits 📝

Good commit messages are like good lab notes - they help everyone (including future you) understand what happened and why.

#### The Basic Format

**One-line commits** for simple changes:

```bash
git commit -m "Fix off-by-one error in gene counting"
git commit -m "Add validation for empty gene families"
git commit -m "Update installation instructions for conda"
```

**Multi-line commits** when you need to explain more:

```bash
git commit -m "Optimize system clustering for large datasets

The previous implementation used nested loops that didn't scale well.
This commit introduces vectorized operations and caching that reduce
runtime from 2 hours to 25 minutes on 10k+ genome datasets.

Tested on: E. coli dataset (15k genomes), Staphylococcus (8k genomes)"
```

#### Commit Message Tips

- **Use imperative mood**: "Add feature" not "Added feature"
- **Be specific**: "Fix memory leak in clustering" beats "Fix bug"
- **Keep the first line under 50 characters** when possible
- **Explain the 'why' not the 'how'** - code shows how, commits explain why
- **Make atomic commits** - one logical change per commit

```{tip}
If you find yourself using "and" in a commit message, you might want to split it into multiple commits!
```

#### Small, Focused Commits

Break your work into digestible pieces:

```bash
# Good: Three clear, focused commits
git commit -m "Add merge method to System class"
git commit -m "Add unit tests for System.merge()"  
git commit -m "Document System.merge() in API reference"

# Less ideal: One big blob
git commit -m "Add merge feature with tests and docs"
```

This makes it easier to review, debug, and potentially revert changes if needed.

### Before You Push: The Checklist ✅

We've all pushed code and then immediately realized something was wrong. This checklist helps catch issues before they become embarrassing! 😅

#### 1. Run the Tests

```bash
# Quick check
pytest

# Full check with coverage (recommended)
pytest --cov=panorama

# Just test what you changed
pytest tests/test_my_feature.py
```

All tests should pass. If something fails, fix it before pushing. Your future self will thank you!

#### 2. Format with Black

We use Black to keep the code style consistent. No more debates about spaces and brackets!

```bash
# Format everything
black panorama/ tests/

# Check what would change (without modifying files)
black --check panorama/ tests/
```

Black makes code reviews smoother since we're focusing on logic, not style.

#### 3. Linting with Flake8

[flake8](https://flake8.pycqa.org/) catches potential bugs and style issues:

```shell
# Check the entire project
flake8 panorama/ tests/

# Check specific files
flake8 panorama/systems/system.py
```

Fix the issues flake8 reports before pushing. Most are quick fixes!

#### 4. Update Documentation

Documentation is code too! If you:

- **Added a new feature** → Update user documentation
- **Changed public APIs** → Update API reference
- **Added/modified functions** → Write/update docstrings

More information on how to write good documentation can be found in the ["how to build the documentation"](buildDoc.md).

#### 5. Review Your Own Changes

Before asking others to review, review yourself:

```bash
# What changed compared to dev?
git diff origin/dev

# Check your commit history
git log origin/dev..HEAD --oneline

# Make sure you didn't leave any debug code
grep -r "print(" panorama/  # Just an example!
```

#### 6. Update the VERSION File

Don't forget to bump the patch version! See the [Versioning section](#versioning-and-releases-🏷️) above.

### Handling Merge Conflicts

Conflicts happen - they're not a failure, just Git asking for your help to combine changes.

```bash
# Start the rebase
git rebase dev
# Git pauses on conflicts

# Open conflicting files and look for markers:
<<<<<<< HEAD
Your code
=======
Their code
>>>>>>> dev

# Edit to keep what you want, then:
git add resolved_file.py
git rebase --continue

# If things get messy, you can always abort and try again
git rebase --abort
```

**Stuck on conflicts?** Don't hesitate to ask for help! Ping a maintainer or open a draft PR and explain where you're stuck.

### Useful Git Commands 🛠️

Some handy commands to make your life easier:

```bash
# Beautiful commit history
git log --oneline --graph --all

# What changed in a specific commit?
git show abc123

# Temporarily save changes without committing
git stash
git stash pop  # Get them back

# Oops, need to change the last commit message?
git commit --amend -m "Better message"

# Interactive rebase to clean up commits before pushing
git rebase -i HEAD~3  # Last 3 commits

# Find which commit introduced a bug
git bisect start
```

---

## Common Pitfalls 🚧

### Mutable Default Arguments

```python
# Bad: Default list is shared between calls!
def add_family(families=[]):
    families.append(new_family)
    return families

# Good: Use None and create new list
def add_family(families=None):
    if families is None:
        families = []
    families.append(new_family)
    return families
```

### Catching All Exceptions

```python
# Bad: Hides all errors, even bugs!
try:
    result = risky_operation()
except:
    pass

# Good: Catch specific exceptions
try:
    result = risky_operation()
except ValueError as e:
    logger.error(f"Invalid value: {e}")
    raise
```

### Hardcoded Paths

```python
# Bad: Won't work on other systems
file_path = "/home/user/data/pangenome.h5"

# Good: Use Path and relative paths
from pathlib import Path
file_path = Path(__file__).parent / "data" / "pangenome.h5"
```

### String Formatting

```python
# Old style (avoid)
message = "Found %d genes in %s" % (count, family_name)

# Good: f-strings (Python 3.6+)
message = f"Found {count} genes in {family_name}"

# Also good: .format() for complex cases
message = "Found {count} genes in {name}".format(count=count, name=family_name)
```