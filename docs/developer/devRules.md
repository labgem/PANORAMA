(dev-rules)=

# Development Methods and Best Practices 💻

Welcome to the PANORAMA development guide! This document covers the coding standards, tools, and practices we follow to
keep the codebase clean, maintainable, and performant. Whether you're fixing a small bug or building a major feature,
these guidelines will help you write code that fits right in.

```{note}
**First time working on PANORAMA?** Don't worry about memorizing everything here. 
Start coding, and refer back as needed. 
The tools we use will catch most issues automatically, and code review will help with the rest. 
```

## Setting Up Your Development Environment

### Install Development Tools

We've bundled all the development tools you need in the `pyproject.toml`:

```shell
# From the PANORAMA root directory
pip install -e .[dev]  # -e for editable mode (recommended)
```

This installs:

- **[Ruff](https://docs.astral.sh/ruff/)** - Code formatter
- **VizTracer** - Performance profiling
- **Vitables** - Interface to open HDF5 files
- And more!

Now you're ready to develop!

## Code Style and Formatting

```{tip}
Most IDE integrates plugins for code formatting and linting.
Look at your IDE's documentation for more information.
```

### PEP 8: Python's Style Guide

We follow as much as possible [PEP 8](https://peps.python.org/pep-0008/),
Python's official style guide. It covers things like:

- Indentation (4 spaces, not tabs)
- Line length (max 79 characters for code, 72 for comments)
- Naming conventions (snake_case for functions, CamelCase for classes)
- Import organization
- Whitespace usage

```{hint}
**The good news?** You don't need to memorize this! Our tools handle most of it automatically.
```

### Ruff: Fast, All-in-One Python Linter and Formatter

We use [Ruff](https://docs.astral.sh/ruff/) for both linting and formatting. Ruff is 10-100x faster than traditional
tools and replaces Black, flake8, isort, and more - all with a single tool and unified configuration.

#### Format and check your code before committing:

```shell
# Format all Python files in the project (replaces Black)
ruff format panorama/ tests/

# Check formatting without modifying files
ruff format --check panorama/ tests/

# Lint the code (replaces flake8)
ruff check panorama/ tests/

# Lint with auto-fix
ruff check --fix panorama/ tests/

# Format and lint a specific file
ruff format panorama/systems/system.py
ruff check panorama/systems/system.py
```

Ruff's formatter is compatible with Black (99% identical output), so you get the same uncompromising style without the debates!
Set up your editor to run Ruff on save:
- **VSCode**: Install the Ruff extension, enable format on save
- **PyCharm**: Use the Ruff plugin
- **Vim**: Use the ruff.vim plugin or ALE

#### Configuration

Ruff is configured in pyproject.toml:
```toml
[tool.ruff]
line-length = 120

[tool.ruff.lint]
select = [
    "E",   # pycodestyle errors
    "W",   # pycodestyle warnings
    "F",   # pyflakes
    "I",   # isort
    "B",   # flake8-bugbear
]
ignore = ["E203"]
```

Fix the issues Ruff reports before pushing. Most are quick fixes!

## Writing Clean, Maintainable Code

### Naming Conventions

We follow Python's naming conventions [PEP 8](https://peps.python.org/pep-0008/#naming-conventions)

### Docstrings: Document Your Code

We choose to use [Google style docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) 
for our code. 

Every public function, class, and module should have a docstring.

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
```

Type hints are especially helpful for:

- Complex function signatures
- When None is a valid return value
- Collections (List, Dict, Set)
- When working with multiple types

### Keep Functions Small

Each function should do one thing well.

## Error Handling and Validation

We try to anticipate and handle errors as early as possible:

- think about what could go wrong
- Choose the Right Exception Type
- Validate Input Early
- Log Important Events - use Python's logging module for informational messages

## Performance and Optimization

### Profile Before Optimizing

```{important}
**"Premature optimization is the root of all evil"** - Donald Knuth
```

Don't guess where your code is slow - measure it!
We use [VizTracer](https://github.com/gaogaotiantian/viztracer) for performance profiling.
But you're free to use any tool you like.

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

## Code Organization

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

### Ask for Help

Development is a collaborative process! When you're stuck:

1. Check existing code for similar patterns
2. Search the issue tracker
3. Ask in discussions or chat
4. Tag a maintainer in your PR
5. Pair program with a colleague
