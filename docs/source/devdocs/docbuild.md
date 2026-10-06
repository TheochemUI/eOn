---
myst:
  html_meta:
    "description": "Instructions for building the eOn documentation with pixi, after installing eonclient into the docs-mta prefix."
    "keywords": "eOn docs, build documentation, Sphinx, MyST, pixi, autodoc-pydantic"
---

# Working with the documentation

`eOn` is a relatively complex project, with both `C++` and `Python` sources,
along with a large set of options.

## Setup

The documentation environment is `docs-mta` in `pixi.toml`. It adds the
`docs`, `develop`, and `metatomic` features. Install eonclient into that
prefix before the book build. The documentation job runs
`meson setup bbdir --prefix=$CONDA_PREFIX --libdir=lib --buildtype release -Dwith_metatomic=True -Dtorch_version=2.10 -Dpip_metatomic=True`
and `meson install --skip-subprojects -C bbdir` under
`pixi run -e docs-mta`, then regenerates the tutorial figures. A missing
figure fails that job.

## Building locally

```{code-block} bash
pixi run -e docs-mta makedocs
```

`makedocs` builds the Sphinx book and writes the sibling doxyYoda C++
API tree to `docs/build/html/api-cpp/`. The book nav and landing page
point at that tree; the Doxygen mainpage points back at the book.

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

Additions to the build go in the `docs` feature of `pixi.toml`.

## Adding citations

Citations live in a `.bib` file exported from Zotero via `better-bibtex`.

Local bibliographies, as noted in [the
documentation](https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#local-bibliographies),
need key prefixes.
