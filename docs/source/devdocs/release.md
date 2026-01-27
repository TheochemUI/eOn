---
myst:
  html_meta:
    "description": "A release checklist for eOn developers, outlining the steps for preparing a new software release, including updating changelogs and versions."
    "keywords": "eOn release checklist, software release, versioning, changelog, towncrier"
---

# Checklist

- [ ] Lockfiles updated

``` bash
uvx pixi-to-conda-lock pixi.lock --output condaEnvs
```

- [ ] Changelog build and edited if necessary

``` bash
pixi r -e dev towncrier build --draft --version=2.8.0 --date "$(date +%Y-%m-%d)"
```

- [ ] Ensure the versions are set consistently
  + [ ] `eon` server
  + [ ] `eonclient` C++ build

- [ ] Generate release notes, and also the publication page, e.g. see `releases/v2.8.0`

To generate the archive since auto-generated archives of tags have [security
issues](https://github.blog/open-source/git/update-on-the-future-stability-of-source-code-archives-and-hashes/)..

```bash
git archive --format=tar v2.9.0 | xz -9 > eon-v2.9.0.tar.xz
```

Paper related tags do not need full release machinery, just a lightweight tag:

```bash
git tag -a arXiv_2510.06030v3 45c077cb
# feat(arxiv): descriptive
#
# more details including the arxiv link
```
