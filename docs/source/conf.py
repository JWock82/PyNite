# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))
# sys.path.insert(0, os.path.abspath('../../PyNite'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Pynite'
copyright = '2023, D. Craig Brinck, SE'
author = 'D. Craig Brinck, SE'
release = '0.0.94'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc']

# Control how class names are displayed
autoclass_content = 'class'

# autodoc settings
autodoc_mock_imports = ['numpy', 'IPython', 'vtk', 'pdfkit']  # Mock import dependencies

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
html_theme = 'bizstyle'
html_theme_options = {
    "sidebarwidth": "300px",
    "navigation_depth": 1
}
html_static_path = ['_static']
