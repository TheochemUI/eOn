---
myst:
  html_meta:
    "description": "Configuration options for calculating the Hessian matrix in EON, used for evaluating prefactors in harmonic transition state theory."
    "keywords": "EON Hessian, vibrational analysis, prefactor, harmonic transition state theory, hTST"
---

# Hessian

The Hessian matrix, in mass-weighted coordinates, is used to evaluate prefactors
in hamonic transition state theory rate calculations.

## Configuration

```{code-block} ini
[Hessian]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.HessianConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: HESS_
keyprefix: hess-
---
```
