# Working with the documentation

`eON` is a relatively complex project, with both `C++` and `Python` sources,
along with a large set of options.

```{versionchanged} 2.1
The input file format transitioned to TOML from previously being an enhanced INI format.
```

## Setup

Although we use `micromamba` for handling system dependencies, the documenation
is handled via the `python` ecosystem. Namely:

- [PDM](https://pdm-project.org/en/latest/) is used to track versions, and provide [development groups](https://pdm-project.org/latest/usage/dependency/#add-development-only-dependencies)

To facilitate interactions with `pdm`,
[pipx](https://pipx.pypa.io/latest/installation/) is recommended.

## Building locally

```{code-block} bash
pipx run pdm run sphinx-build -b html docs/source docs/build/html
```

This can be viewed locally with an HTTP server.

```{code-block} bash
python -m http.server docs/build/html
```

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