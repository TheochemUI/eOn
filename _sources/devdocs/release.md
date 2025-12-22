---
myst:
  html_meta:
    "description": "A release checklist for EON developers, outlining the steps for preparing a new software release, including updating changelogs and versions."
    "keywords": "EON release checklist, software release, versioning, changelog, towncrier"
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
