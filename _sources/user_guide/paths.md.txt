---
myst:
  html_meta:
    "description": "Configuration for specifying the directory paths for calculations and output in eOn."
    "keywords": "eOn paths, configuration, directories, file paths, output"
---

# Paths

```{admonition} Legacy configuration
:class: note
These path settings are primarily used by the server-client architecture.
For standalone client usage (the common case), the default paths are sufficient
and this section can be left unconfigured.
```

The akmc, basin hopping, parallel replica, and escape rate drivers
create the `results` directory before they open the log or `lockfile`.

## Configuration

```{code-block} ini
[Paths]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.PathsConfig
```
