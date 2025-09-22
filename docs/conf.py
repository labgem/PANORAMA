# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Imports -----------------------------------------------------------------

import datetime
from pathlib import Path

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

# -- Project information -----------------------------------------------------

project = "PANORAMA"
organization = "LABGeM"
author = "Jérôme Arnoux"
year = datetime.date.today().year
copyright = f"{'2025' if year == 2025 else f'2025-{year}'}, {organization}, {author}"

html_show_copyright = True

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
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.duration",
    "sphinx_search.extension",
    "sphinx.ext.napoleon",  # Extension for NumPy and Google style docstrings
    "sphinx.ext.extlinks",
    "sphinx.ext.todo",  # Remove warning todo
    "sphinx_design",
    "sphinxcontrib.jquery",
    "sphinxcontrib.mermaid",
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
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    "requirements.txt",
    "developper/figure_script.md",
    "developper/draw_spot_script.md",
    "developper/draw.md",
    "developper/system_asso_script.md",
    "developper/write_flat_script.md",
    "developper/conserved_spot_script.md"
]

suppress_warnings = [
    "myst.header",
    # "autosectionlabel.*",
    "toc.not_included"
]

# The name of the default role for inline references
default_role = "py:obj"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'
html_theme = "pydata_sphinx_theme"
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Source Buttons
html_context = {
    "github_user": "labgem",
    "github_repo": "PANORAMA",
    "github_version": "main",
    "doc_path": "docs",
}

html_logo = (
    "https://labgem.genoscope.cns.fr/wp-content/uploads/2021/06/GENOSCOPE-LABGeM.jpg"
)

# These are options specifically for the Wagtail Theme.
html_theme_options = {
    # Navigation depth and collapsing sidebars
    "show_nav_level": 2,
    "navigation_depth": 10,
    # Page Table of Contents
    "show_toc_level": 10,
    "secondary_sidebar_items": ["page-toc", "edit-this-page", "sourcelink"],
    "use_edit_page_button": True,
    # Header links
    "external_links": [
        # {
        #     "url": "https://doi.org/10.1093/bioinformatics/btad214",
        #     "name": "Paper",
        # },
        {
            "url": "https://github.com/labgem/PPanGGOLiN",
            "name": "PPanGGOLiN",
        }
    ],
    "header_links_before_dropdown": 5,  # default value
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/labgem/PANORAMA",
            "icon": "fa-brands fa-github",
        },
        # {
        #     "name": "PyPI",
        #     "url": "https://pypi.org/project/panorama",
        #     "icon": "fa-custom fa-pypi",
        # },
    ],
    "icon_links_label": "Quick Links",
    # Sphinx indices
    "primary_sidebar_end": ["indices.html"],
    # Annoucement banners
    "announcement": "PANORAMA just released!",
    # Back to Top button
    "back_to_top_button": True,
    # Branding and logo
    "logo": {
        "text": "PANORAMA ",
        "alt_text": "PANORAMA documentation - Home",
        "image_light": "https://labgem.genoscope.cns.fr/wp-content/uploads/2021/06/GENOSCOPE-LABGeM.jpg",
        "image_dark": "https://labgem.genoscope.cns.fr/wp-content/uploads/2021/06/GENOSCOPE-LABGeM.jpg",
    },
    # Configure pygments theme
    "pygments_light_style": "tango",
    "pygments_dark_style": "monokai",
    "navbar_align": "left",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
    # "footer_links": ",".join(
    #     [
    #         "About Us|http://example.com/",
    #         "Contact|http://example.com/contact",
    #         "Legal|http://example.com/dev/null",
    #     ]
    # ),
    "navigation_with_keys": "True",
}
# html_show_sphinx = False

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "PANORAMA"

# -- Options for napoleon extension ------------------------------------------

napoleon_include_init_with_doc = True
napoleon_include_special_with_doc = True
napoleon_include_private_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_rtype = False
napoleon_custom_sections = [('Returns', 'params_style'), ('Raises', 'params_style')]
