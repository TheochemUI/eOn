---
myst:
  html_meta:
    "description": "Instructions for building and working with the eOn documentation, including setup with uv, local building, and adding citations and extensions."
    "keywords": "eOn docs, build documentation, Sphinx, MyST, uv, autodoc-pydantic"
---

# Working with the documentation

`eOn` is a relatively complex project, with both `C++` and `Python` sources,
along with a large set of options.

```{versionchanged} 3.1_TBA
The input file format transitioned to TOML from previously being an enhanced INI format.
```

## Setup

Although we use `pixi` for handling system dependencies, the documentation
is handled via the `python` ecosystem. Namely:

- [uv](https://docs.astral.sh/uv/) is used to manage the lockfile and dependency groups defined in `pyproject.toml`

## Building locally

````{margin}
```{note}
`uv` manages dependency groups declared in `pyproject.toml` under
`[dependency-groups]`.  System-level dependencies (compilers, etc.)
are handled separately by `pixi`.
```
````

```{code-block} bash
# Setup dependencies
uv sync --group docs --no-install-project
# Need to install for autodoc-pydantic
uv pip install . -vvv --no-build-isolation
uv run sphinx-build -b html docs/source docs/build/html
```

This can be viewed locally with an HTTP server.

```{code-block} bash
python -m http.server docs/build/html
```

### Extra-index pitfall

The `pixi.toml` workspace sets `extra-index-urls` pointing at the PyTorch CPU
wheel index and uses `index-strategy = "unsafe-best-match"`.  When pixi
activates an environment, these propagate as `UV_EXTRA_INDEX_URL` and
`PIP_EXTRA_INDEX_URL` environment variables.

If you regenerate `uv.lock` from inside a pixi shell, `uv lock` will resolve
packages against the PyTorch index as well as PyPI.  The PyTorch index carries
only a subset of wheels (often only the latest CPython pre-release), so
packages like `markupsafe` can resolve to versions with no compatible wheel for
the CI Python version.

To avoid this, always regenerate the lockfile with the extra-index variables
unset:

```{code-block} bash
UV_EXTRA_INDEX_URL="" PIP_EXTRA_INDEX_URL="" uv lock --upgrade
```

Alternatively, run `uv lock` outside of any pixi environment.

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

Additions to the build process are handled by the `docs` dependency group, so additions are done via:

```{code-block} bash
uv add --group docs "sphinxcontrib-bibtex"
```

## Adding citations

Citations are handled in a `.bib` file which is exported via `better-bibtex`
with Zotero. Kindly do not modify these by hand.

Note that because we need local bibliographies, as noted in [the
documentation](https://sphinxcontrib-bibtex.readthedocs.io/en/latest/usage.html#local-bibliographies)
we need to use key prefixes.
