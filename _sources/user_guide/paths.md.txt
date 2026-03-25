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

## Configuration

```{code-block} ini
[Paths]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.PathsConfig
```
