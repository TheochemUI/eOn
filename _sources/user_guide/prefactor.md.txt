---
myst:
  html_meta:
    "description": "Guide to harmonic transition state theory (hTST) prefactor calculation in eOn."
    "keywords": "eOn prefactor, hTST, harmonic transition state theory, rate constant, Vineyard"
---

# Prefactor

```{versionchanged} 2.12
Fixed Hessian size validation bug where `min2Freqs` was never checked.
Fixed float precision in QQ-HTST physical constants.
```

The prefactor job computes the harmonic transition state theory (hTST) rate
prefactor using the Vineyard formula. This relates vibrational frequencies at
a minimum and saddle point to determine the attempt frequency for a transition.

The prefactor job requires pre-computed minimum and saddle point structures.
It uses the [Hessian](project:hessian.md) to obtain vibrational frequencies
at each stationary point.

## Usage

```{code-block} ini
[Main]
job = prefactor

[Prefactor]
rate_estimation = htst
filter_scheme = fraction
filter_fraction = 0.9
```

The `filter_scheme` controls which frequencies contribute to the prefactor.
Setting `filter_fraction = 0.9` includes frequencies within 90% of the full
range, filtering out very high-frequency modes that may be numerical artifacts.

## Configuration

```{code-block} ini
[Prefactor]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.PrefactorConfig
```
