# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

#This is needed to allow the highlighting when the docs are built online
from pygments.lexers import PythonLexer
from pygments.lexers.special import TextLexer
from sphinx.highlighting import lexers

# Map ipython3 (not recognized) to Python3
lexers['ipython3'] = PythonLexer()
lexers['ipy'] = PythonLexer()  # optional alias
###########################

import sys
from pathlib import Path

sys.path.insert(0, str(Path('..', 'src').resolve()))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Imodulator'
copyright = '2025, Duarte Silva, Kaan Sünnetçioğlu'
author = 'Duarte Silva, Kaan Sünnetçioğlu'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
    'sphinxcontrib.bibtex',
	'sphinx.ext.todo',
	'nbsphinx',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

collapse_navigation = False

# -- Configurations for todo extension -----------------------------------
todo_include_todos = True

# -- Configurations for napoleon extension -----------------------------------
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True

# -- Configurations for autodoc-typehints extension ---------------------------
always_use_bars_union = True
typehints_defaults = 'braces'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 3,
}

html_sidebars = {
    '**': [
        'globaltoc.html',   # the main sidebar (global navigation)
        'localtoc.html',    # add a new sidebar with the current page’s headings
        'searchbox.html',   # keep the search box
    ],
}

# -- Options for the bibtex extension ----------------------------------------
bibtex_bibfiles = ['refs.bib']

# -- Options for nbsphinx extension ------------------------------------------
# Optional: prevent notebooks from executing during build
nbsphinx_execute = 'never'

# Recommended: ignore notebook checkpoints
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# -- Mock imports for autodoc -----------------------------------------------
autodoc_mock_imports = [
    "gmsh",
	"pygmsh"
]