# Documentation Build Guide 📚

This guide helps developers build and maintain the PANORAMA documentation before merging changes to the main branch. The
documentation system uses Sphinx with the PyData theme and MyST parser for enhanced Markdown support.

```{danger}
When merging or creating pull requests to main, ReadTheDocs will automatically detect changes and update the online 
documentation. Ensure your documentation builds without errors and follows the established formatting standards before 
submitting.
```

```{note}
The documentation framework is built on [Sphinx](https://www.sphinx-doc.org/en/master/) with the community-supported 
[PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/index.html). 
MyST parser enables enhanced Markdown features including cross-references, admonitions, and code execution.
```

## Environment Setup 🛠️

### Prerequisites

- Python environment with PANORAMA installed
- Sphinx and documentation dependencies
- Access to PANORAMA source code

### Installation

Documentation requirements are specified in the `sphinx_requirements.txt` file and included in the project's
`pyproject.toml`. Install using pip:

```shell
# From PANORAMA root directory
pip install .[doc]  # Install with documentation dependencies
```

**Alternative installation:**

```shell
# Direct installation from requirements file
pip install -r docs/sphinx_requirements.txt
```

## Building Documentation 🏗️

### Standard HTML Build
```{note}
For now this method does not give the actual pretty results, with CSS and JS.
We encourage to use the [Live Development Server](#live-development-server) instead.
```
Generate static HTML documentation using Sphinx's build command:

```shell
# Navigate to documentation directory
cd $PANORAMA_ROOT/docs/

# Build HTML documentation
sphinx-build -b html . build/
```

**Using Makefile (recommended):**

```shell
cd $PANORAMA_ROOT/docs/
make html
```

**Clean build (recommended for troubleshooting):**

```shell
make clean && make html
```

### Live Development Server

Use `sphinx-autobuild` for real-time documentation preview during development:

```shell
cd $PANORAMA_ROOT/docs/
sphinx-autobuild . build/ 
```

**Features:**

- **Automatic Refresh**: Changes trigger immediate page reloads
- **Live Preview**: View modifications without manual rebuilds
- **Error Detection**: Real-time build error notifications
- **Cross-platform**: Works on all major operating systems

```{note}
The `readthedocs-sphinx-search` package only functions on ReadTheDocs hosting. 
Local builds will show "[INFO] Docs are not being served on Read the Docs" - this is expected behavior.
```

### Build Verification

After building, verify the documentation:

1. **Check Build Logs**: Review console output for warnings or errors
2. **Open HTML Files**: Navigate to `build/index.html` in a web browser
3. **Test Navigation**: Verify all links and cross-references work correctly
4. **Validate Formatting**: Ensure code blocks, tables, and admonitions render properly

## Modifying Existing Documentation ✏️

### User Documentation

User-facing documentation files are located in `source/user/`:

**File Types:**

- **Command Documentation**: One file per PANORAMA command
- **Installation Guides**: Setup instructions for different platforms

**Editing Guidelines:**

- Use clear, action-oriented headings
- Include practical command-line examples
- Add troubleshooting sections for complex procedures
- Reference other documentation sections using MyST cross-references

### Modeler Documentation

Model creation documentation is in `source/modeler/`:
**Content Areas:**

- Model Definition: Guidelines for defining biological system models
- Model Format Specifications: File formats and data structures for models
- Validation Procedures: Model testing and quality assurance protocols
- Integration Workflows: How to incorporate new models into PANORAMA
- Best Practices: Recommended approaches for model development and maintenance

### Developer Documentation

Development-focused content is in `source/developer/`:

**Content Areas:**

- **Code Guidelines**: PEP compliance and style standards
- **Git Workflows**: Branch management and contribution processes
- **Testing**: Unit test creation and CI/CD integration
- **Architecture**: Core system design and extension points

### API Documentation

API documentation is automatically generated from docstrings but requires manual updates for new modules:

```shell
# Regenerate API documentation
sphinx-apidoc -o source/api $PANORAMA_ROOT/panorama -f

# Convert RST to Markdown (if using rst2myst)
rst2myst convert source/api/*.rst
rm source/api/*.rst
```

```{hint}
[rst2myst](https://rst-to-myst.readthedocs.io/en/latest/) is a tool for converting RST to Markdown. 
It can be incompatible with the environment. 
If so, we suggest to create another environment to install it.
```
## Adding New Documentation 📄

### Command Documentation Workflow

1. **Create Command File**:
   ```shell
   touch source/user/command_name.md  # or other name that seems relevant
   ```

2. **Add Content Structure**:
   ```markdown
   # Command Name 🧬
   
   Brief description of command functionality.
   
   ## Command Line Usage 🚀
   ## Arguments ⚙️  
   ## Output Files 📁
   ## Examples 💡
   ## Troubleshooting 🛠️
   ```

3. **Update Table of Contents**:
   Add to `source/user/index.md` toctree:
   ```markdown
        ```{toctree}
        :maxdepth: 2
        ...  
        user/command_name
        ```
   ```

### Developer Guidelines

For new development documentation:

1. **Create Guideline File**: Use descriptive filename reflecting content focus
2. **Follow Template Structure**: Match the existing developer documentation format
3. **Add Cross-references**: Link to related API documentation and user guides
4. **Update Index**: Add to developer documentation toctree

## Documentation Standards 📋

### Content Guidelines

**User Documentation:**

- One file per PANORAMA command with descriptive titles
- Installation instructions for all supported platforms
- Issue reporting and enhancement request procedures
- Avoid referencing internal code functions

**Developer Documentation:**

- PEP compliance and coding standards
- Git workflow and GitHub contribution guidelines
- Unit testing procedures and CI/CD integration
- Core architecture documentation for central components
- Extension and plugin development guides

**API Documentation:**

- Comprehensive docstring coverage for all public functions
- Type hints and parameter descriptions
- Usage examples for complex functions
- Cross-references between related components

### Formatting Standards

**Headings:**

- Use descriptive, action-oriented headings
- Follow consistent emoji usage for visual hierarchy
- Maintain proper heading levels (H1 > H2 > H3)

**Code Examples:**

- Include complete, runnable command examples
- Add explanatory comments for complex operations
- Show expected output where helpful
- Use proper syntax highlighting

**Cross-references:**

- Link related sections using MyST references
- Reference API documentation from user guides
- Maintain a consistent link text and formatting
