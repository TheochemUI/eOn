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

`info.txt` under the results directory is replaced in one step. A write that
fails leaves the previous file.

## Configuration

```{code-block} ini
[Paths]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.PathsConfig
```
