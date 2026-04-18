---
myst:
  html_meta:
    "description": "Guide to the Lanczos eigenmode method in eOn for saddle point searches."
    "keywords": "eOn Lanczos, lowest eigenmode, saddle point search, min-mode following"
---

# Lanczos

```{versionchanged} 2.12
Force call tracking fixed to use the job's potential instance.
```

The Lanczos method determines the lowest curvature mode of the potential energy
surface using the iterative Lanczos algorithm
{cite:t}`lcz-malekDynamicsLennardJonesClusters2000`. It is an alternative to
the [dimer method](project:dimer.md) for min-mode following saddle searches.

## Force call costs per eigenmode computation

Each min-mode method has different force evaluation costs per eigenmode
computation step:

| Method | Force calls per step | How it works |
|--------|---------------------|--------------|
| Lanczos | 1 + N (N = Lanczos iters, typically 3--7) | Builds Krylov subspace; 1 initial + 1 per iteration |
| Improved Dimer | 2 + N (N = rotations, up to 20) | Forward/backward FD + 1 per trial rotation |
| ARTn | 1 per `artn_step` call | External force eval; internal Lanczos is implicit |

The Lanczos method converges the lowest eigenmode in fewer force calls because
the Krylov basis extracts more curvature information per evaluation than
repeated dimer rotations.

### Benchmark: LJ38 TS optimization (100 structures)

```{versionadded} 2.13
```

| Method | Avg force calls | Median | Min | Max | Failures |
|---------|-----------------|--------|-----|------|----------|
| Lanczos | 238 | 178 | 84 | 1625 | 0/100 |
| ARTn | 426 | 269 | 100 | 2462 | 10/100 |
| Dimer | 523 | 359 | 151 | 2996 | 0/100 |

Force calls are counted via `Potential::forceCallCounter` (incremented on every
`Matter::computePotential` call). Cached evaluations (positions unchanged) do
not count.

## When to Use Lanczos vs Dimer

Both methods find the same lowest eigenmode. The Lanczos method is generally
more efficient (fewer force calls) while the Dimer is more robust and
integrates with GP acceleration:

- **Lanczos**: 1 gradient evaluation per iteration. Fewer total force calls.
  Natural for large systems. Integrates with OCINEB climbing image refinement.
- **Improved Dimer**: 2 evaluations per iteration (finite difference). More
  force calls but integrates with the [AtomicGPDimer](project:dimer.md) for
  GP-accelerated searches.
- **ARTn**: Push-based exploration. Finds different saddles than
  gradient-following methods. See [ARTn](project:artn.md).

## Usage

The Lanczos method is selected in the saddle search configuration:

```{code-block} ini
[Saddle Search]
min_mode_method = lanczos

[Lanczos]
tolerance = 0.01
max_iterations = 20
```

`max_iterations` controls how many Lanczos iterations are performed per
eigenmode computation. Unlike the dimer's `rotations_max`, each Lanczos
iteration costs exactly 1 gradient evaluation.

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
