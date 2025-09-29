# How to Contribute to PANORAMA 🤝

Welcome! We're thrilled that you're interested in contributing to PANORAMA. Whether you're fixing a typo, reporting a
bug, or building a major new feature, your contribution matters and we're here to support you through the process.

This guide will walk you through everything you need to know to contribute effectively. Don't worry if this seems like a
lot - you don't need to read everything at once! Start with what's relevant to your contribution and refer back as
needed.

## Quick Start for Different Contributions 🚀

### I Want to Fix Something or Add a Feature 💻

Ready to code? Here's your roadmap:

1. **Set up your development environment** → See [Development Setup](#development-setup-🛠️)
2. **Learn our git workflow** → See [Version Control Guide](link-to-git-guide)
3. **Follow coding standards** → See [Development Methods](link-to-dev-methods)
4. **Write tests** → See [Testing Guidelines](link-to-testing-guide)
5. **Update documentation** → See [Documentation Guide](link-to-doc-guide)
6. **Submit your Pull Request** → See [Pull Request Process](#pull-request-process-🔄)

### I Want to Improve Documentation 📚

Documentation improvements are always welcome!

1. **Find what needs improving** - Typos, unclear sections, missing examples
2. **Follow the documentation guide** → See [Building Documentation](link-to-doc-guide)
3. **Submit a Pull Request** with your changes

Even small documentation fixes are valuable - don't hesitate!

## Development Setup 🛠️

### Prerequisites

You'll need:

- **Python 3.8+** (we recommend 3.10 or newer)
- **Git** for version control
- A **GitHub account** for submitting contributions

### Clone the Repository

```shell
# Clone your fork (after forking on GitHub)
git clone https://github.com/YOUR_USERNAME/PANORAMA.git
cd PANORAMA

# Add the upstream repository
git remote add upstream https://github.com/labgem/PANORAMA.git
```

### Install in Development Mode

```shell
# Install with all development tools
pip install -e .[dev,test,doc]

# Or install specific extras
pip install -e .[dev]   # Just development tools
pip install -e .[test]  # Just testing tools
pip install -e .[doc]   # Just documentation tools
```

The `-e` flag installs in "editable" mode, so your code changes take effect immediately!

### Get Test Data

```shell
# Clone the test data repository
git clone https://github.com/labgem/PANORAMA_test

# Set the environment variable
export PANORAMA_TEST_DATA_PATH=PANORAMA_test/
```

## Our Development Workflow 🔄

We use a structured workflow to keep the codebase stable and organized:

### Branch Structure

- **`main`** - Stable, production-ready code
- **`dev`** - Latest development version (your starting point!)
- **Feature branches** - Your work happens here

### Creating a Feature Branch

```shell
# Start from dev
git checkout dev
git pull upstream dev

# Create your feature branch
git checkout -b feature/your-feature-name
```

See our [Git Version Control Guide](link-to-git-guide) for detailed workflow instructions.

## Coding Standards 📝

We follow consistent standards to keep the codebase maintainable:

### Code Style

- **Follow PEP 8** - Python's official style guide
- **Use Black** for formatting - Runs automatically, no debates!
- **Write docstrings** - Document your public functions and classes
- **Use type hints** - Helps catch bugs and improves readability

```shell
# Format your code before committing
black panorama/ tests/

# Check for style issues
flake8 panorama/ tests/
```

### Error Handling

- **Validate inputs** early in functions
- **Raise appropriate exceptions** (ValueError, TypeError, etc.)
- **Handle errors gracefully** - Don't let exceptions crash unexpectedly
- **Log important events** using Python's logging module

### Performance

- **Profile before optimizing** - Use VizTracer to find bottlenecks
- **Don't guess** where the slow parts are - measure!
- **Readable code first** - Only optimize when profiling shows a need

For detailed guidelines, see our [Development Methods Guide](link-to-dev-methods).

## Writing Tests 🧪

Tests are crucial - they catch bugs and give everyone confidence that changes work correctly.

### Test Requirements

- **New features must include tests** - Both unit and functional tests
- **Bug fixes should include regression tests** - Prevent the bug from coming back
- **Aim for high coverage** - But focus on meaningful tests, not just numbers

### Running Tests

```shell
# Run all tests
pytest

# Run with coverage
pytest --cov=panorama

# Run specific tests
pytest tests/test_systems.py
```

### Test Types

- **Unit tests** - Test individual functions and classes in isolation
- **Functional tests** - Test complete workflows and CLI commands

See our [Testing Guidelines](link-to-testing-guide) for comprehensive testing practices.

## Updating Documentation 📚

Documentation is just as important as code! When you make changes:

### What to Update

- **API documentation** - Auto-generated from docstrings, but may need regeneration
- **User guides** - If you add user-facing features
- **Developer docs** - If you change development processes or architecture

### Building Documentation Locally

```shell
cd docs
sphinx-autobuild source/ build/
# Open http://127.0.0.1:8000 in your browser
```

See our [Documentation Building Guide](link-to-doc-guide) for details.

## Pull Request Process 🔄

### Before Submitting

Make sure you've completed this checklist:

- ❓ Code is formatted with Black
- ❓ Tests pass locally (`pytest`)
- ❓ New features have tests
- ❓ Documentation is updated
- ❓ Commit messages are clear
- ❓ VERSION file is bumped (patch number)

### Creating the Pull Request

1. **Push your branch** to your fork
2. **Create a Pull Request** on GitHub targeting the `dev` branch
3. **Fill out the PR template** with details about your changes
4. **Link related issues** using "Fixes #123" or "Related to #456"

### PR Title Format

Use clear, descriptive titles:

- ✅ "Add support for CRISPR-Cas system detection"
- ✅ "Fix memory leak in pangenome clustering"
- ❌ "Update" or "Changes"

### Review Process

- **Be patient** - Reviews take time, especially for complex changes
- **Respond to feedback** - Address comments and questions
- **Be open to suggestions** - Reviewers want to help improve your code
- **Ask questions** if feedback is unclear

All PRs need at least one maintainer approval before merging.

## Versioning 🏷️

PANORAMA uses semantic versioning: `Major.Minor.Patch`

- **Patch** - Your PR increments this (bug fixes, small improvements)
- **Minor** - When `dev` merges into `main` (new features, accumulated changes)
- **Major** - Breaking changes (rare, discussed with maintainers)

Don't forget to update the VERSION file in your PR!

## Code Review Guidelines 👥

### For Authors

- **Respond promptly** to review comments
- **Don't take it personally** - Reviews are about code, not you
- **Mark conversations resolved** once you've addressed them
- **Ask for clarification** if feedback is confusing

### For Reviewers

When reviewing others' code:

- **Be constructive** and specific in feedback
- **Explain reasoning** behind suggestions
- **Acknowledge good work** when you see it
- **Distinguish** between required changes and suggestions

## Getting Help 🆘

Stuck on something? Here's what to do:

### Where to Ask

- **GitHub Discussions** - General questions and discussions
- **Issue tracker** - Bug reports and feature requests
- **Pull Request comments** - Questions about specific code
- **Tag maintainers** - We're here to help!

### What to Include

When asking for help, provide:

- What you're trying to do
- What you've already tried
- Error messages (complete, not snippets)
- Relevant code or commands
- Your environment (OS, Python version, PANORAMA version)

**Remember:** No question is too basic! We all started somewhere, and we want to support your contribution journey.

## Community Guidelines 🌍

### Be Respectful

- **Treat everyone with respect** - We're all here to make PANORAMA better
- **Be patient** - Everyone has different experience levels and other projects in progress
- **Be constructive** - Focus on ideas and code
- **Assume good intentions** - Most misunderstandings are just that

### Be Inclusive

We welcome contributors of all backgrounds and experience levels. 
Whether this is your first open source contribution or your thousandth, you're welcome here!

## Recognition 🌟

We value all contributions! Contributors are:

- Listed in the project's contributor list
- Mentioned in release notes for significant contributions
- Building their open source portfolio

Your work helps researchers worldwide - thank you! 🙏

## Quick Links 🔗

### Essential Guides

- [Git Version Control](link-to-git-guide) - Branch workflow, commits, and PR process
- [Development Methods](link-to-dev-methods) - Coding standards and best practices
- [Testing Guidelines](link-to-testing-guide) - Writing and running tests
- [Documentation Building](link-to-doc-guide) - Building and updating docs

### External Resources

- [GitHub Repository](https://github.com/labgem/PANORAMA)
- [Issue Tracker](https://github.com/labgem/PANORAMA/issues)
- [Test Data Repository](https://github.com/labgem/PANORAMA_test)
- [ReadTheDocs](link-to-readthedocs) - Official documentation

## Quick Reference 📋

```shell
# Setup
git clone https://github.com/YOUR_USERNAME/PANORAMA.git
cd PANORAMA
pip install -e .[dev,test,doc]

# Start working
git checkout dev
git pull upstream dev
git checkout -b feature/my-feature

# Before committing
black panorama/ tests/
flake8 panorama/ tests/
pytest --cov=panorama

# Update version
echo "1.2.6" > VERSION  # Increment patch number

# Push and create PR
git push -u origin feature/my-feature
# Create PR on GitHub
```

## Ready to Contribute? 🎉

Pick an issue labeled `good first issue` or `help wanted`, or tackle something you've noticed that needs fixing. Don't
hesitate to ask questions - we're here to help you succeed!

Thank you for contributing to PANORAMA. Your work helps advance bioinformatics research worldwide! 🚀