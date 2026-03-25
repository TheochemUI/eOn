---
myst:
  html_meta:
    "description": "Installing LAMMPS for use with eOn."
    "keywords": "eOn LAMMPS, installation"
---

# LAMMPS

```{versionchanged} 2.12
LAMMPS is no longer a build-time dependency. eOn loads `liblammps` at
runtime via `dlopen`. Simply install LAMMPS in the same environment.
```

## Quick Install

```{code-block} bash
# conda-forge (recommended)
conda install -c conda-forge lammps

# Or pixi
pixi add lammps
```

No eOn rebuild required. See the [LAMMPS potential guide](project:../user_guide/lammps_pot.md)
for usage details.
