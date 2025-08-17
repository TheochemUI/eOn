# Dimer method

The dimer method of {cite:t}`dm-henkelmanDimerMethodFinding1999` with
improvements by {cite:t}`dm-heydenEfficientMethodsFinding2005` and
{cite:t}`dm-kastnerSuperlinearlyConvergingDimer2008` for estimating the lowest
Eigenmode using only first derivatives.

An overview may be found in {cite:t}`dm-olsenComparisonMethodsFinding2004`.

The dimer separation is set in the **[Main]** section with the
`finiteDifference` parameter.

The method of {cite:t}`dm-melanderRemovingExternalDegrees2015` is also
implemented for use with gas phase systems.

```{note}
There is no point removing rotations for an extended system. Rotation removal may be more detrimental as noted in {cite:t}`goswamiBayesianHierarchicalModels2025a`.
```

```{versionadded} 2.5_TBA
The Gaussian Process Regression accelerated dimer in C++ from {cite:t}`goswamiEfficientImplementationGaussian2025a`.
```

## Configuration

```{code-block} toml
[Dimer]
```

```{eval-rst}
.. autopydantic_model:: eon.schema.DimerConfig
```

## References

```{bibliography}
---
style: alpha
filter: docname in docnames
labelprefix: DM_
keyprefix: dm-
---
```
