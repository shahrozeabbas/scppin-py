# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'scPPIN-py'
copyright = '2025, scPPIN-py Contributors'
author = 'scPPIN-py Contributors'
release = '0.3.0'
version = '0.3.0'

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # Support for Google and NumPy style docstrings
    'sphinx.ext.viewcode',  # Add links to source code
    'sphinx.ext.intersphinx',  # Link to other projects' documentation
    'sphinx.ext.autosummary',  # Generate summary tables
    'sphinx.ext.autodoc.typehints',  # Type hints support
    'sphinx_copybutton',  # Copy button for code blocks
    'myst_parser',  # Markdown support
]

# Napoleon settings for parsing docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Autodoc settings
autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented_params'
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

# Autosummary settings
autosummary_generate = True

# Intersphinx mapping - links to external documentation
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'igraph': ('https://python-igraph.readthedocs.io/en/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'scanpy': ('https://scanpy.readthedocs.io/en/stable/', None),
    'anndata': ('https://anndata.readthedocs.io/en/stable/', None),
}

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Furo theme (dark, modern, minimal)
html_theme = 'furo'
html_static_path = ['_static']

# Furo theme options
html_theme_options = {
    "dark_css_variables": {
        "color-brand-primary": "#9575cd",
        "color-brand-content": "#9575cd",
        "color-background-primary": "#181818",
        "color-background-secondary": "#1e1e1e",
        "color-background-border": "#2e2e2e",
        "color-foreground-primary": "#ffffff",
        "color-foreground-secondary": "#d0d0d0",
        "color-foreground-muted": "#a0a0a0",
        "color-foreground-border": "#404040",
        "color-inline-code-background": "#2e2e2e",
        "color-sidebar-background": "#161616",
        "color-sidebar-background-border": "#2e2e2e",
        "color-sidebar-search-background": "#1e1e1e",
        "color-sidebar-item-background": "#2e2e2e",
    },
    "light_css_variables": {
        "color-brand-primary": "#7b68ee",
        "color-brand-content": "#7b68ee",
    },
    "top_of_page_button": None,
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "navigation_style": "sidebar",
}

html_title = 'scPPIN-py Documentation'
html_logo = None  # Can add logo later
html_favicon = None  # Can add favicon later

# -- Options for MyST parser -------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "substitution",
]

