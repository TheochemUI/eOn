---
myst:
  html_meta:
    "description": "Documentation for the structure comparison settings in EON, used to configure similarity measures and equivalence thresholds."
    "keywords": "EON structure comparison, point group, bond order, configuration comparison"
---

# Structure Comparison

## Configuration

```{code-block} ini
[Structure Comparison]
```

```{versionchanged} 3.1_TBA
In TOML, this will be `[Structure_Comparison]`
```


```{eval-rst}
.. autopydantic_model:: eon.schema.StructureComparisonConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: SC_
keyprefix: sc-
---
```
