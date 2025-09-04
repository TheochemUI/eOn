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
