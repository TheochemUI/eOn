# Working with the documentation

`eON` is a relatively complex project, with both `C++` and `Python` sources,
along with a large set of options.

```{versionchanged} 2.1_TBA
The input file format transitioned to TOML from previously being an enhanced INI format.
```

## Setup

Although we use `micromamba` for handling system dependencies, the documenation
is handled via the `python` ecosystem. Namely:

- [PDM](https://pdm-project.org/en/latest/) is used to track versions, and provide [development groups](https://pdm-project.org/latest/usage/dependency/#add-development-only-dependencies)

To facilitate interactions with `pdm`,
[pipx](https://pipx.pypa.io/latest/installation/) is recommended.

## Building locally

````{margin}
```{note}
`pipx` simplifies running Python commands, and `pdm` handles version updates better, syncing nicely with the `pyproject.toml` so no Python dependencies not needed by the client should be in the `environment.yml`
```
````

```{code-block} bash
# Setup dependencies
pipx run pdm sync
# Need to install for autodoc-pydantic
pipx run pdm run pip install . -vvv
pipx run pdm run sphinx-build -b html docs/source docs/build/html
```

This can be viewed locally with an HTTP server.

```{code-block} bash
python -m http.server docs/build/html
```

## Writing documentation

We use `myst` markdown via the `myst-parser` extension for almost everything,
however, the `pydantic` schema is handled by `autodoc-pydantic` which requires
`rst` directives only, so:
- The docstrings for the configuration are formatted with `rst`.
- `eval-rst` is required to wrap the configuration stanzas in the `myst`
  markdown files

## Additions

The following sections detail methods to add functionality to the documentation.

## Adding extensions

Additions to the build process are handled by the `pdm` development group `docs`, so additions are done via:

```{code-block} bash
pipx run pdm add -dG docs "sphinxcontrib-bibtex"
```

## Adding citations

Citations are handled in a `.bib` file which is exported via `better-bibtex`
with Zotero. Kindly do not modify these by hand.

Note that because we need local bibliographies, as noted in [the
documentation](https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#local-bibliographies)
we need to use key prefixes.
