# Git branching model

I'm using the Git branching model described
[here](http://nvie.com/posts/a-successful-git-branching-model/).

# Documentation

Documentation is handled by Sphinx. To force sphinx to rebuild all the html files run `python -m sphinx-build -M html . ./_build -a` in the cmd when you're in the docs folder.

Make sure you also run the following in your environment:

`pip install sphinxcontrib-bibtex`
`pip install sphinx_autodoc_typehints`
`pip install sphinx-rtd-theme`
`pip install sphinx`
`pip install nbsphinx`

# How to release a new version

When you want to release a new version, follow these directions.

1. Make sure you merge proper commits to from the development to the main branch
2. Bump the version number in the `pyproject.toml`