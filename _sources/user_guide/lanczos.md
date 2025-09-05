---
myst:
  html_meta:
    "description": "Configuration guide for the Lanczos method in EON, used for determining the lowest eigenmode of a system."
    "keywords": "EON Lanczos, lowest eigenmode, saddle point search, dynamics"
---

# Lanczos

The Lanczos method for determining the lowest eigenmode along the lines of {cite:t}`lcz-malekDynamicsLennardJonesClusters2000`.

## Configuration

```{code-block} ini
[Lanczos]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.LanczosConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: LCZ_
keyprefix: lcz-
---
```
