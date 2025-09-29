# Development Methods and Best Practices 💻

Welcome to the PANORAMA development guide! This document covers the coding standards, tools, and practices we follow to keep the codebase clean, maintainable, and performant. Whether you're fixing a small bug or building a major feature, these guidelines will help you write code that fits right in.

```{note}
**First time working on PANORAMA?** Don't worry about memorizing everything here. Start coding, and refer back as needed. The tools we use will catch most issues automatically, and code review will help with the rest. We're here to support you!
```

## Setting Up Your Development Environment 🛠️

### Install Development Tools

We've bundled all the development tools you need in the `pyproject.toml`:

```shell
# From the PANORAMA root directory
pip install -e .[dev]  # -e for editable mode (recommended)
```

This installs:
- **Black** - Code formatter
- **VizTracer** - Performance profiling
- **Vitables** - Interface to open HDF5 files
- And more!

Now you're ready to develop !

## Code Style and Formatting 🎨

### PEP 8: Python's Style Guide

We follow [PEP 8](https://peps.python.org/pep-0008/), Python's official style guide. It covers things like:

- Indentation (4 spaces, not tabs)
- Line length (max 79 characters for code, 72 for comments)
- Naming conventions (snake_case for functions, CamelCase for classes)
- Import organization
- Whitespace usage

**The good news?** You don't need to memorize this! Our tools handle most of it automatically.

### Black: The Uncompromising Code Formatter

We use [Black](https://black.readthedocs.io/) to format all Python code. Black makes formatting decisions for you, so you can focus on logic instead of style debates.

**Format your code before committing:**

```shell
# Format all Python files in the project
black panorama/ tests/

# Check what would change without modifying files
black --check panorama/ tests/

# Format a specific file
black panorama/systems/system.py
```

Black is opinionated, which is actually great - no more discussions about where to put line breaks!

```{hint}
Set up your editor to run Black on save:
- **VSCode**: Install the Python extension, enable format on save
- **PyCharm**: Use the Black plugin
- **Vim**: Use the black.vim plugin
```

### Linting with flake8

While Black handles formatting, [flake8](https://flake8.pycqa.org/) catches potential bugs and style issues:

```shell
# Check the entire project
flake8 panorama/ tests/

# Check specific files
flake8 panorama/systems/system.py
```

Common issues flake8 catches:
- Unused imports
- Undefined variables
- Lines that are too long
- Overly complex functions
- Missing or extra whitespace

Fix the issues flake8 reports before pushing. Most are quick fixes!

## Writing Clean, Maintainable Code 📝

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

### Docstrings: Document Your Code

```{important}
We choose to use [Google style docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) 
for our code.
```

Every public function, class, and module should have a docstring:

```python
def detect_systems(pangenome, models, annotation_source):
    """
    Detect systems in a pangenome using provided models.

    This function scans the pangenome for gene families matching the provided
    system models and returns detected systems.

    Args:
        pangenome (Pangenome): The pangenome to analyze
        models (List[Model]): System models to search for
        annotation_source (str): Source of functional annotations to use

    Returns:
        List[System]: Detected systems ordered by completeness score

    Raises:
        ValueError: If annotation_source is not found in pangenome
        TypeError: If models is not a list
    """
    pass
```

**For simple functions, a one-line docstring is fine:**
```python
def get_gene_count(family):
    """Return the number of genes in the family."""
    return len(family.genes)
```


### Type Hints

Use type hints to make your code more readable and catch bugs early:

```python
from typing import List, Optional, Set, Dict

def merge_gene_families(
    families: List[GeneFamily],
    threshold: float = 0.95
) -> Set[GeneFamily]:
    """Merge similar gene families above the threshold."""
    pass

def find_system(
    system_id: str,
    pangenome: Pangenome
) -> Optional[System]:
    """Find a system by ID, or None if not found."""
    pass
```

Type hints are especially helpful for:
- Complex function signatures
- When None is a valid return value
- Collections (List, Dict, Set)
- When working with multiple types

## Error Handling and Validation ⚠️

### Handle Errors Gracefully

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

### Choose the Right Exception Type

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

### Validate Input Early

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

### Log Important Events

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

## Performance and Optimization 🚀

### Profile Before Optimizing

**"Premature optimization is the root of all evil"** - Donald Knuth

Don't guess where your code is slow - measure it! We use [VizTracer](https://github.com/gaogaotiantian/viztracer) for performance profiling.

### Using VizTracer

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
- Functions called very frequently
- Unexpected I/O operations
- Nested loops that could be optimized

### Common Performance Patterns

**Use generators for large datasets:**
```python
# Good: Memory efficient
def iter_gene_families(pangenome):
    for family in pangenome.families:
        yield family

# Less efficient: Loads everything into memory
def get_all_families(pangenome):
    return [family for family in pangenome.families]
```

**Cache expensive computations:**
```python
from functools import lru_cache

@lru_cache(maxsize=1000)
def calculate_similarity(family_id1, family_id2):
    """Calculate similarity with caching."""
    # Expensive computation here
    pass
```

**Use built-in functions (they're optimized in C):**
```python
# Fast
total = sum(len(family) for family in families)

# Slower
total = 0
for family in families:
    total += len(family)
```

**Use sets for membership testing:**
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

### When to Optimize

Optimize when:
- Profiling shows a clear bottleneck
- Users report performance issues
- You're working with very large datasets

Don't optimize when:
- Code is already fast enough
- It would make code much harder to read
- The slow part is I/O bound (disk/network)

Remember: **Readable code > Clever code**

## Code Organization 📁

### Module Structure

Organize code logically:

```
panorama/
├── systems/           # System detection and analysis
│   ├── __init__.py
│   ├── system.py     # Core System classes
│   ├── models.py     # Model definitions
│   └── detection.py  # Detection algorithms
├── annotation/        # Functional annotation
│   ├── __init__.py
│   └── annotate.py
└── utility/            # Shared utilities
    ├── __init__.py
    └── translate       # Module for translating models from set of sources
        └──macsymodel_translator.py
        └──padloc_translator.py
```

### Keep Functions Focused

Each function should do one thing well:

```python
# Good: Clear single purpose
def filter_low_coverage_families(families, min_coverage):
    """Filter gene families below minimum coverage."""
    return [f for f in families if f.coverage >= min_coverage]

def sort_families_by_size(families):
    """Sort families by number of genes."""
    return sorted(families, key=lambda f: len(f.genes), reverse=True)

# Less ideal: Doing too much
def process_families(families, min_coverage):
    """Filter and sort families."""
    filtered = [f for f in families if f.coverage >= min_coverage]
    sorted_families = sorted(filtered, key=lambda f: len(f.genes), reverse=True)
    return sorted_families
```

Small, focused functions are:
- Easier to test
- Easier to reuse
- Easier to understand
- Easier to debug

## Development Workflow 🔄

### Before Committing

Run this checklist before every commit:

```shell
# 1. Format code
black panorama/ tests/

# 2. Check linting
flake8 panorama/ tests/

# 3. Run tests
pytest

# 4. Check test coverage
pytest --cov=panorama

# Now you're ready to commit!
git add .
git commit -m "Your clear, descriptive message"
```

### Code Review Checklist

When reviewing code (yours or others'), check:

- **Functionality**: Does it work as intended?
- **Tests**: Are there adequate tests?
- **Readability**: Is the code clear and well-documented?
- **Performance**: Any obvious bottlenecks?
- **Error handling**: Are edge cases handled?
- **Style**: Does it follow our conventions?

## Common Pitfalls to Avoid 🚧

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

## Debugging Tips 🔍

### Use the Debugger

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

## Getting Better at Development 📚

### Learning Resources

- **PEP 8**: [peps.python.org/pep-0008](https://peps.python.org/pep-0008/)
- **Black docs**: [black.readthedocs.io](https://black.readthedocs.io/)
- **VizTracer**: [github.com/gaogaotiantian/viztracer](https://github.com/gaogaotiantian/viztracer)
- **Python docs**: [docs.python.org](https://docs.python.org/)

### Ask for Help

Development is a collaborative process! When you're stuck:

1. Check existing code for similar patterns
2. Search the issue tracker
3. Ask in discussions or chat
4. Tag a maintainer in your PR
5. Pair program with a colleague

No question is too basic. We all started somewhere, and we're here to help you grow as a developer!

## Quick Reference 📋

```shell
# Setup
pip install -e .[dev]

# Before committing
black panorama/ tests/
flake8 panorama/ tests/
pytest --cov=panorama

# Profiling
viztracer script.py
vizviewer result.json

# Debugging
python -m pdb script.py
```

**Key principles:**
- Follow PEP 8 (let Black handle formatting)
- Write docstrings for public APIs
- Handle errors explicitly
- Profile before optimizing
- Keep functions focused
- Write tests

Happy coding! 🚀