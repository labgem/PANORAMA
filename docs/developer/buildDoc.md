(build-doc)=

# Building the Documentation 📚

This guide will help you build and preview the PANORAMA documentation locally before merging your changes. Whether
you're fixing a typo or adding a whole new section, testing your docs locally saves time and catches issues early!

```{danger}
When you merge your branch into `main`, ReadTheDocs will automatically rebuild and deploy the documentation online.
 Make sure everything looks good locally first - broken docs are visible to everyone! 
```

## Installing Documentation Dependencies 📦

The documentation requires some specific packages to build. We've made this easy for you!

All required packages are listed in the [sphinx_requirements.txt](../sphinx_requirements.txt) file. But here's an even
simpler way - the [pyproject.toml](../../pyproject.toml) includes everything you need:

```shell
# From the PANORAMA root directory
pip install -e .[doc]  # -e for editable mode (recommended for development)

# Or without editable mode
pip install .[doc]
```

This installs PANORAMA plus all the documentation tools: Sphinx, MyST-Parser, the PyData Sphinx theme, and more.

## Building and Previewing Documentation 🔨

Want to see your changes in real-time as you write? Use `sphinx-autobuild` (installed with the doc dependencies):

```shell
cd docs  # Navigate to the docs folder
sphinx-autobuild source/ build/

# The server will start and give you a URL like:
# Serving on http://127.0.0.1:8000
# Open this URL in your browser
```

The documentation will automatically rebuild whenever you save a file. Refresh your browser to see changes! 🎉

```{tip}
Sometimes changes don't show up immediately. If this happens, try to clean your build with `make clean`.
```

```{note}
You might see a message about [readthedocs-sphinx-search](https://readthedocs-sphinx-search.readthedocs.io/en/latest/) 
not working locally - that's expected! This package only works on the live ReadTheDocs site. 
The message `[INFO] Docs are not being served on Read the Docs, readthedocs-sphinx-search will not work.` is harmless.
```

## Modifying Existing Documentation ✏️

### User and Developer Guides

To update existing documentation, open the relevant Markdown file and edit it! The structure is straightforward:

- **User documentation**: `docs/user/`
- **Modeler documentation**: `docs/modeler/`
- **Developer documentation**: `docs/developer/`
- **API reference**: `docs/api/`

Save your changes and check them in your browser with `sphinx-autobuild` running.

### API Documentation

**Good news** - the API documentation updates automatically when you change docstrings in the code! Modify your
docstrings and rebuild.

**However**, if you add a new package or module, you'll need to regenerate the API docs. See
the [Updating API documentation](#updating-api) section below.

(heading-adding)=

## Adding to Existing Documentation ➕

### Adding New Command Documentation

When you add a new command to PANORAMA, document it for users:

1. **Create a new Markdown file** in `docs/source/user/` named after your command (e.g., `detect_systems.md`)
2. **Write clear documentation** with examples, parameters, and use cases
3. **Add it to the table of contents** by editing the `toctree` in the User Guide index:

```markdown
    ```{toctree}
    user/installation
    user/quickstart
    user/your_new_command  # Add this line (without .md extension)
    ```
```

### Adding New Developer Guidelines

Got ideas for improving the developer docs? Awesome! We welcome contributions.

- **If it fits an existing guide**, add your content there
- **If it needs its own file**, create a new Markdown file with a descriptive name
- **Add it to the Developer Guide toctree** so people can find it

We're pretty flexible about what goes into developer docs - if you think it'll help someone, it probably belongs!

(updating-api)=

### Updating API Documentation

When you add new packages or modules to PANORAMA, regenerate the API reference:

```shell
# From the PANORAMA root directory
sphinx-apidoc -o docs/source/api panorama/ -f
```

```{attention}
The `sphinx-apidoc` command generates ReStructuredText (.rst) files, but we use Markdown. 
Convert them using the instructions in the [RST to Markdown section](#rst2md) below.
```

## Creating Documentation from Scratch 🏗️

```{warning}
This section is for major documentation overhauls. If you're considering this, 
please discuss with other maintainers first! Usually, you just need to modify existing docs.
```

If you really need to start fresh, here's how:

### Initialize with Sphinx

```shell
DOCS=path/to/documentation/folder
sphinx-quickstart $DOCS
```

Follow the prompts:

- **Separate source and build directories?** → `y`
- **Project name** → `PANORAMA`
- **Author name** → Your name
- **Project release** → Current version number
- **Project language** → `en` (or press Enter)

This creates the basic structure: `source/`, `build/`, `Makefile`, etc.

### Configuration File

Replace the contents of `source/conf.py` with this configuration:

```python
# Configuration file for the Sphinx documentation builder.
#
# For the full list of options, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from pathlib import Path

# -- Project information -----------------------------------------------------

project = 'PANORAMA'
copyright = '2023, LABGeM'
author = '{your name}'

# Get version from VERSION file
release = open(Path(__file__).resolve().parents[2] / "VERSION").read().rstrip()

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    'sphinx_search.extension',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
```

(rst2md)=

### Converting ReStructuredText to Markdown

Sphinx uses ReStructuredText (.rst) by default, but we prefer Markdown because it's simpler and more widely used. Here's
how to convert:

```{note}
The `rst-to-myst` package isn't compatible with our Sphinx version, so use a separate environment for the conversion.
```

```shell
# In a separate environment
pip install rst-to-myst[sphinx]

# Convert the files
rst2myst convert source/index.rst
# or for multiple files
rst2myst convert api/*.rst

# Back in your PANORAMA environment, clean up
rm source/index.rst  # or api/*.rst
```

### Including the README

Want to include the project README in your documentation? Smart move - less duplication!

Add this to your `index.md`:

```markdown
    ```{include} ../../README.md
    :relative-images:
    ```
```

### Documentation Structure

Here's how we organize PANORAMA docs:

#### User Documentation

User docs should be practical and example-driven:

1. **One file per command** - Each command gets its own guide
2. **Installation guide** - Help users get started
3. **Contributing guide** - How to report issues or request features
4. **No code references** - Focus on usage, not implementation

Write for bioinformaticians who want to use PANORAMA, not necessarily code it.

#### Developer Documentation

Developer docs are for people working on PANORAMA's code:

1. **PEP standards** - Code style and Python conventions
2. **Git workflow** - How we use version control (you're reading one now!)
3. **Testing guide** - Writing and running tests
4. **Documentation guide** - How to improve these docs (meta!)
5. **Architecture deep-dives** - Explain complex parts of the codebase. Feel free to reference code, classes, and
   implementation details here.

#### API Documentation

Generate API docs automatically from your docstrings and reference API elements in your docs: `{ref}\package panorama\`

```{tip}
The `sphinx.ext.autosectionlabel` extension can create duplicate label warnings when rebuilding. 
If this happens, just rename or remove duplicate section titles in the affected files.
```

```{tip}
`sphinx-apidoc` creates a `modules.md` file that's not actually used. Feel free to delete it to avoid warnings.
```

## Getting Help 🆘

Documentation giving you trouble? Here's what to do:

1. **Check Sphinx documentation** - [sphinx-doc.org](https://www.sphinx-doc.org/)
2. **Check MyST documentation** - [mystmd.org/guide](https://mystmd.org/guide)
3. **Check MyST parser documentation** - [myst-parser.readthedocs.io](https://myst-parser.readthedocs.io/en/latest/index.html)
4. **Check PyData Sphinx theme documentation** - [pydata-sphinx-theme.readthedocs.io](https://pydata-sphinx-theme.readthedocs.io/en/stable/index.html)
5. **Ask in discussions** - Other contributors can help
6. **Open a draft PR** - Get feedback on your documentation changes

Remember: good documentation is just as important as good code. Thanks for taking the time to document your work! 🙏