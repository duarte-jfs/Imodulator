# Git branching model

I'm using the Git branching model described
[here](http://nvie.com/posts/a-successful-git-branching-model/).

# Code formatting

This project uses [Ruff](https://docs.astral.sh/ruff/) as its formatter. The rules live
in `[tool.ruff]` / `[tool.ruff.format]` in `pyproject.toml` and are auto-discovered — no
editor setup is required.

- Install the pinned Ruff: `pip install -e ".[dev]"`
- Format: `ruff format .` (or a single file, e.g. `ruff format src/imodulator/PhotonicPolygon.py`)
- Check without writing: `ruff format --check .`

Ruff discovers the project config by walking up from the file being formatted, so any
editor that runs `ruff format` picks up the same rules. Personal CLI `--config`/flag
overrides win over the project file — drop those when working on this repo so your
formatting matches everyone else's.

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

(https://packaging.python.org/en/latest/tutorials/packaging-projects/#generating-distribution-archives)

1. Make sure you merge proper commits to from the development to the main branch
2. Bump the version number in the `pyproject.toml`
3. Bump the version number in `docs/conf.py`
4. Build the docs locally
5. Build the docs in ReadTheDocs
6. Generate the distribution archives by running

> python -m build

7. Upload the .tar.gz and .whl files to the github release section
>pip install --upgrade twine

>twine upload --repository testpypi dist/*
8. Test uploading the package to TestPyPI to make sure it all goes well
9. If it all went ok, then go to PyPI and hope for the best