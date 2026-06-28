---
myst:
  html_meta:
    "description": "Documentation for the Dimer method in eOn, used for finding transition states and lowest eigenmodes using only first derivatives."
    "keywords": "eOn Dimer method, transition state search, saddle point, eigenmode"
---

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
There is no point removing rotations for an extended system. Rotation removal may be more detrimental as noted in {cite:t}`dm-goswamiBayesianHierarchicalModels2025a`.
```

```{versionadded} 2.5
The Gaussian Process Regression accelerated dimer in C++ from {cite:t}`dm-goswamiEfficientImplementationGaussian2025a`.
```

## Rotation backends

The softest-mode estimate used in the dimer rotation step can be chosen with
`rotation_backend` under **[Dimer]** (default `classical`). Min-mode-following
translation along that mode is unchanged; only how τ / the lowest curvature is
obtained differs.

| Value | Mode estimation | Notes |
| --- | --- | --- |
| `classical` | Constrained dimer rotation (Heyden / Kästner–Sherwood style) | Default; uses `opt_method`, torque limits, and rotation budgets as before |
| `lanczos` | Finite-difference Lanczos min-mode | Shares the client Lanczos implementation; see {doc}`lanczos` |
| `davidson` | Finite-difference Davidson min-mode | Same role as Lanczos with a different iterative subspace |
| `lor` | Locally optimal rotation (LOR), Algorithm I | {cite:t}`dm-lengEfficientSoftestMode2013`; at most one new FD force per rotation iteration via Hessian–vector products and force translation of prior H·N / H·F⊥ products |

```{versionadded} 2.15
`rotation_backend` (`classical` \| `lanczos` \| `davidson` \| `lor`). The LOR path
implements Leng et al. {cite:p}`dm-lengEfficientSoftestMode2013`: 2×2 then 3×3 Ritz
problems in the rotation subspace, force translation for H·P₃ (linear action of H
on the unit Gram–Schmidt residual of the trial direction, not Gram–Schmidt on
H·P in ambient space), residual / `rotations_max` / stall stops, and a best-mode
restore when the Ritz sequence is non-monotonic under FD noise.
```

Example:

```{code-block} ini
[Dimer]
improved = True
rotation_backend = lor
```

`improved = True` (default) uses `ImprovedDimer`; non-classical backends skip the
classical `IDimerRot` loop and call the selected min-mode backend through a
shared dispatch helper. Rotation iteration budgets for LOR follow
`rotations_max` (not geometry `max_iterations`).

## Configuration

```{code-block} ini
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
