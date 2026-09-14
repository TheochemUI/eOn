---
myst:
  html_meta:
    "description": "eOn v3.2.1: IDPP honors frozen atoms, pyeonclient optimizer methods, rgpot 3.1 D3/D4/Expr, ASE PBC wrap order."
    "keywords": "eOn v3.2.1, IDPP constraints, opt_method, DFT-D3, ExprPot, ASE PBC"
---

## [v3.2.1] - 2026-09-13

Patch on `v3.2.0`. IDPP path init keeps frozen atoms fixed
(TheochemUI/eOn#410). pyeonclient can set `opt_method`,
`neb_opt_method`, and refine (#406). rgpot wrap is **v3.2.0** with
D3/D4, ExprPot, and MOPACPot. ASE `from_ase` applies PBC before wrap.

```{toctree}
:maxdepth: 2
:caption: Release notes

release-notes
```
