# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from pathlib import Path

# -- Project information -----------------------------------------------------

project = "PANORAMA"
copyright = "2025, LABGeM"
html_show_copyright = True
author = "Jérôme Arnoux"

# The full version, including alpha/beta/rc tags
release = (
    open(Path(__file__).resolve().parents[1] / "VERSION").read().rstrip()
)  # Get release number in the VERSION file


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    # "sphinxcontrib.jquery",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx_search.extension",
    "sphinxcontrib.mermaid",
    "sphinx_wagtail_theme",
]

source_suffix = {".md": "markdown"}

# Prefix document path to section labels, to use:
# `path/to/file:heading` instead of just `heading`
autosectionlabel_prefix_document = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

suppress_warning = ["myst.header", "autosectionlabel.*"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'
html_theme = "sphinx_wagtail_theme"
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# These are options specifically for the Wagtail Theme.
html_theme_options = {
    "project_name": "PANORAMA",
    "logo": "img/wagtail-logo-circle.svg",
    "logo_alt": "Wagtail",
    "logo_height": 59,
    "logo_url": "/",
    "logo_width": 45,
    "github_url": "https://github.com/jpjarnoux/PANORAMA",
    # "header_links": "Top 1|http://example.com/one, Top 2|http://example.com/two",
    # "footer_links": ",".join(
    #     [
    #         "About Us|http://example.com/",
    #         "Contact|http://example.com/contact",
    #         "Legal|http://example.com/dev/null",
    #     ]
    # ),
}
# html_show_sphinx = False
