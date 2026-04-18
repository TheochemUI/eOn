---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  display_name: Python 3
  language: python
  name: python3
myst:
  html_meta:
    "description": "Tutorial on visualizing eOn optimization trajectories -- minimization, saddle search, and NEB -- using rgpycrumbs and chemparseplot."
    "keywords": "eOn visualization, trajectory plotting, NEB landscape, saddle search, minimization, rgpycrumbs, chemparseplot"
---

# Visualizing Optimization Trajectories

eOn can output structured per-iteration data alongside trajectory movies when
`write_movies = true` is set in the `[Debug]` section of `config.ini`. These
files enable rich 2D reaction-valley projections and convergence plots via
[rgpycrumbs](https://rgpycrumbs.rgoswami.me) and
[chemparseplot](https://chemparseplot.rgoswami.me).

This tutorial walks through the full workflow: running eOn with PET-MAD on the
vinyl alcohol -> acetaldehyde keto-enol tautomerization, then generating energy
profiles, convergence panels, and 2D optimization landscapes.

## Setup

```{code-cell} python
:tags: [remove-cell]

import os
import subprocess
from contextlib import chdir
from pathlib import Path
import tempfile

os.environ["RGPYCRUMBS_AUTO_DEPS"] = "1"

# Ensure eonclient can find metatensor/torch shared libs at runtime
import sys
import torch, metatensor.torch, metatomic.torch
_lib_paths = {os.path.join(sys.prefix, "lib")}
for _mod in [torch, metatensor.torch, metatomic.torch]:
    _d = os.path.dirname(_mod.__file__)
    _lib_paths.add(_d)
    _lib = os.path.join(_d, "lib")
    if os.path.isdir(_lib):
        _lib_paths.add(_lib)
_extra = ":".join(_lib_paths)
os.environ["LD_LIBRARY_PATH"] = f"{_extra}:{os.environ.get('LD_LIBRARY_PATH', '')}"

# Download PET-MAD xs v1.5.0 model
repo_id = "lab-cosmo/upet"
tag = "v1.5.0"
url = f"https://huggingface.co/{repo_id}/resolve/main/models/pet-mad-xs-{tag}.ckpt"
model_dir = Path("models")
model_dir.mkdir(exist_ok=True)
model_path = model_dir / f"pet-mad-xs-{tag}.pt"

if not model_path.exists():
    subprocess.run(["mtt", "export", url, "-o", str(model_path)], check=True)
print(f"Model: {model_path}")
```

```{code-cell} python
:tags: [remove-cell]

from rgpycrumbs.eon.helpers import write_eon_config

# Tutorial working directory
tutorial_dir = Path(tempfile.mkdtemp(prefix="eon_tutorial_"))
print(f"Working in: {tutorial_dir}")
```

## Minimization

We start by minimizing a slightly distorted vinyl alcohol molecule to
demonstrate the minimization trajectory output.

```{code-cell} python
:tags: [remove-cell]

import shutil
import numpy as np
import readcon

# Read the vinyl alcohol reactant and perturb it
src = Path("data/reactant.con")
atoms = readcon.read_con_as_ase(str(src))[0]
np.random.seed(42)
atoms.positions += np.random.randn(*atoms.positions.shape) * 0.1

min_dir = tutorial_dir / "minimization"
min_dir.mkdir()
readcon.write_con(str(min_dir / "pos.con"), [readcon.ConFrame.from_ase(atoms)])

min_config = {
    "Main": {"job": "minimization"},
    "Potential": {"potential": "metatomic"},
    "Metatomic": {"model_path": str(model_path.resolve())},
    "Optimizer": {"opt_method": "lbfgs", "converged_force": 0.0514221, "max_iterations": 200},
    "Debug": {"write_movies": True},
}
write_eon_config(min_dir, min_config)

with chdir(min_dir):
    result = subprocess.run(["eonclient"], capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"eonclient failed (exit {result.returncode})")
print(f"Minimization done: {list(min_dir.glob('minimization.*'))}")
```

### Energy profile

```{code-cell} python
from rgpycrumbs.eon.plt_min import main as plt_min_main
from click.testing import CliRunner

runner = CliRunner()
runner.invoke(plt_min_main, [
    "--job-dir", str(min_dir), "--prefix", "minimization",
    "--plot-type", "profile", "--dpi", "150",
])
```

### Convergence panel

```{code-cell} python
runner.invoke(plt_min_main, [
    "--job-dir", str(min_dir), "--prefix", "minimization",
    "--plot-type", "convergence", "--dpi", "150",
])
```

## Nudged Elastic Band

We run an OCI-NEB calculation on the vinyl alcohol -> acetaldehyde
keto-enol tautomerization using PET-MAD.

```{code-cell} python
:tags: [remove-cell]

neb_dir = tutorial_dir / "neb"
neb_dir.mkdir()
shutil.copy("data/reactant.con", neb_dir / "reactant.con")
shutil.copy("data/product.con", neb_dir / "product.con")

neb_config = {
    # Mirror eon_orchestrator's vinyl-alcohol NEB config
    # (https://github.com/HaoZeke/eon_orchestrator examples/vinyl_alcohol).
    # 17 images + SIDPP path + energy-weighted springs + climbing image +
    # OCI-NEB MMF refinement converges the keto-enol tautomerization
    # barrier (~6 eV with PET-MAD-XS) cleanly.
    "Main": {"job": "nudged_elastic_band", "random_seed": 706253457},
    "Potential": {"potential": "metatomic"},
    "Metatomic": {"model_path": str(model_path.resolve()), "device": "cpu"},
    "Nudged Elastic Band": {
        "images": 17,
        # initialization
        "initializer": "sidpp",
        "sidpp_growth_alpha": 0.33,
        "minimize_endpoints": "false",
        # energy-weighted springs
        "energy_weighted": "true",
        "ew_ksp_min": 0.972,
        "ew_ksp_max": 9.72,
        "ew_trigger": 0.5,
        # climbing image once roughly converged
        "climbing_image_method": "true",
        "climbing_image_converged_only": "true",
        "ci_after": 0.5,
        "ci_after_rel": 0.8,
        # OCI-NEB MMF refinement at the saddle (recommended,
        # see Goswami et al. 2026)
        "ci_mmf": "true",
        "ci_mmf_after": 0.1,
        "ci_mmf_after_rel": 0.8,
        "ci_mmf_penalty_strength": 1.5,
        "ci_mmf_penalty_base": 0.4,
        "ci_mmf_angle": 0.9,
        "ci_mmf_nsteps": 1000,
    },
    "Optimizer": {
        "opt_method": "lbfgs",
        "max_move": 0.1,
        # 10^-3 Ha/Bohr; standard chemistry precision
        "converged_force": 0.0514221,
        "max_iterations": 1000,
    },
    "Debug": {"write_movies": True},
}
write_eon_config(neb_dir, neb_config)
with chdir(neb_dir):
    result = subprocess.run(["eonclient"], capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
        raise RuntimeError(f"eonclient failed (exit {result.returncode})")
print(f"NEB done: {len(list(neb_dir.glob('neb_*.dat')))} steps")
```

### Energy profile

```{code-cell} python
from rgpycrumbs.eon.plt_neb import main as plt_neb_main

runner.invoke(plt_neb_main, [
    "--input-dat-pattern", str(neb_dir / "neb_*.dat"),
    "--con-file", str(neb_dir / "neb.con"),
    "--ira-kmax", "14",
    "--force-recompute",
    "--plot-type", "profile",
    "--plot-structures", "crit_points",
    "--strip-renderer", "xyzrender",
    "--perspective-tilt", "8",
    "--zoom-ratio", "0.3",
    "--dpi", "150",
])
```

### 2D reaction landscape

```{code-cell} python
runner.invoke(plt_neb_main, [
    "--input-dat-pattern", str(neb_dir / "neb_*.dat"),
    "--input-path-pattern", str(neb_dir / "neb_path_*.con"),
    "--con-file", str(neb_dir / "neb.con"),
    "--ira-kmax", "14",
    "--force-recompute",
    "--plot-type", "landscape",
    "--surface-type", "grad_imq",
    "--show-pts",
    "--landscape-path", "all",
    "--project-path",
    "--plot-structures", "crit_points",
    "--strip-renderer", "xyzrender",
    "--strip-dividers",
    "--perspective-tilt", "8",
    "--show-legend",
    "--dpi", "150",
])
```

## See also

For producing your own NEB trajectories beyond this tutorial:

- The [`eon-pet-neb`](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
  example in the [lab-cosmo/atomistic-cookbook](https://github.com/lab-cosmo/atomistic-cookbook):
  a step-by-step PET-MAD NEB walkthrough using ASE for path setup. This is
  the canonical metatomic-consumer integration test for eOn.
- [HaoZeke/eon_orchestrator](https://github.com/HaoZeke/eon_orchestrator):
  Snakemake-orchestrated workflow for batches of NEB calculations with
  PET-MAD, including IRA pre-alignment, endpoint minimization, energy
  profiles, and 2D RMSD landscapes. The vinyl alcohol config used in this
  tutorial is taken from `examples/vinyl_alcohol/`.

## Further reading

- [rgpycrumbs documentation](https://rgpycrumbs.rgoswami.me) for advanced
  options (structure rendering, multiple trajectory overlays, custom themes)
- {doc}`/user_guide/minimization` and {doc}`/user_guide/saddle_search` for
  output file details
- {doc}`dict_config` for programmatic configuration
