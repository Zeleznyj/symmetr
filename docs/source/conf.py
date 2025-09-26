# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Symmetr'
copyright = '2025, Jakub Zelezny'
author = 'Jakub Zelezny'
release = '0.9.10'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",   # optional, for Google/NumPy style docstrings
    "sphinx.ext.viewcode",   # optional, adds links to highlighted source code
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    "collapse_navigation": False,  # for sphinx_rtd_theme
    "sticky_navigation": True,
    "navigation_depth": 3,         # only top-level sections in sidebar
    "titles_only": True
}

# Exclude imported members and include members from __all__
autodoc_default_options = {
    'members': True,
    'imported-members': False,  # Ignore imported items
    'special-members': None,
    'exclude-members': '__weakref__',  # Exclude specific members if needed
}

# Generate autosummary pages
#autosummary_generate = True
#autosummary_imported_members = False  # Don't include imported members
