---
myst:
  html_meta:
    "description": "A release checklist for eOn developers, outlining the steps for preparing a new software release, including updating changelogs and versions."
    "keywords": "eOn release checklist, software release, versioning, changelog, towncrier"
---

# Checklist

- [ ] Bump version in `pyproject.toml` -- this is the single source of truth.
  All other consumers (`docs/source/conf.py`, `client/get_version.py`,
  and the meson/CMake builds) read from it automatically.
  + [ ] Also bump `version` in `pixi.toml` to match.

- [ ] Lockfiles updated

``` bash
uvx pixi-to-conda-lock pixi.lock --output condaEnvs
```

- [ ] Changelog build and edited if necessary

``` bash
pixi r -e dev towncrier build --draft --version=X.Y.Z --date "$(date +%Y-%m-%d)"
```

- [ ] Generate release notes, and also the publication page, e.g. see `releases/vX.Y.Z/`

To generate the archive since auto-generated archives of tags have [security
issues](https://github.blog/open-source/git/update-on-the-future-stability-of-source-code-archives-and-hashes/):

```bash
git archive --format=tar vX.Y.Z | xz -9 > eon-vX.Y.Z.tar.xz
```

Paper related tags do not need full release machinery, just a lightweight tag:

```bash
git tag -a arXiv_2510.06030v3 45c077cb
# feat(arxiv): descriptive
#
# more details including the arxiv link
```
