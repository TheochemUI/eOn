---
myst:
  html_meta:
    "description": "The eOn finite_difference job, a curvature scan at the least-coordinated atom."
    "keywords": "eOn finite difference, curvature.dat, epicenter"
---

# Finite-difference curvature

`job = finite_difference` reads `pos.con`. It picks the least-coordinated
atom and displaces that atom together with its neighbors inside
`[Structure Comparison] neighbor_cutoff`. A frozen coordinate is left
out of the displacement and out of its norm. The client writes
`curvature.dat`, with columns `dR` and `curvature`, and `results.dat`.
The displacements use step sizes `1e-7`, `1e-6`, `1e-5`, `1e-4`,
`1e-3`, `5e-3`, `0.01`, `0.05`, and `0.1`. This scan is not the Hessian
`fd_scheme`, and `[Main] finite_difference` does not set these steps.
The schema spelling `finite_differences` does not select this job.

```{code-block} ini
[Main]
job = finite_difference
```
