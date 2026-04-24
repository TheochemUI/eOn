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

eOn can embed structured per-iteration metadata directly into trajectory movie
frames when `write_movies = true` is set in the `[Debug]` section of
`config.ini`. These metadata-rich `.con` movies enable rich 2D reaction-valley
projections and convergence plots via
[rgpycrumbs](https://rgpycrumbs.rgoswami.me) and
[chemparseplot](https://chemparseplot.rgoswami.me).

For one compatibility window, `write_deprecated_outs = true` can also be used
to emit the older minimization and saddle-search `.dat` sidecars. NEB continues
to write its existing `neb.dat` / `neb_*.dat` outputs alongside the `.con`
movies.

This tutorial walks through two complete workflows with PET-MAD:

- a single-ended minimization of a perturbed vinyl alcohol structure
- an OCI-NEB calculation for vinyl alcohol -> acetaldehyde

For each workflow we run `eonclient`, then visualize the result through the
public `rgpycrumbs` dispatcher so the tutorial matches current command-line
usage.

## Setup

```{code-cell} python
:tags: [remove-cell]

import os
import subprocess
import sys
from pathlib import Path
import tempfile

os.environ["RGPYCRUMBS_AUTO_DEPS"] = "1"

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

```{code-cell} python
:tags: [remove-cell]

from IPython.display import Image, display

plot_dir = tutorial_dir / "plots"
plot_dir.mkdir(exist_ok=True)

def run_cli_or_raise(args, *, cwd=None, timeout=300):
    print("$", " ".join(args))
    result = subprocess.run(
        args,
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"command failed (exit {result.returncode}): {' '.join(args)}")
    return result


def run_rgpycrumbs(group, command, *args, cwd=None, timeout=600):
    # 600s caps first-call rgpycrumbs overhead on the weakest GHA runner.
    # Landscape plots (grad_matern / grad_imq surface fits) dominate wall
    # time; profile and convergence panels land well under a minute locally,
    # but cold matplotlib + CLI import on ubuntu-latest has been observed
    # past the old 180s ceiling.
    return run_cli_or_raise(
        [sys.executable, "-m", "rgpycrumbs.cli", group, command, *args],
        cwd=cwd,
        timeout=timeout,
    )


def run_eon(job_dir, *, timeout=300):
    eonclient = str(Path(sys.executable).with_name("eonclient"))
    result = subprocess.run(
        [eonclient],
        cwd=job_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if result.returncode != 0:
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        raise RuntimeError(f"eonclient failed (exit {result.returncode})")
    return result


def show_plot(path):
    display(Image(filename=path))
```

## Minimization

We start by minimizing a slightly distorted vinyl alcohol molecule to
demonstrate the single-ended optimization outputs consumed by
`rgpycrumbs eon plt-min`.

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
    # Current plt-min profile/convergence commands still consume the legacy
    # minimization.dat sidecar, while the 2D landscape uses the metadata-rich
    # minimization.con movie.
    "Debug": {"write_movies": True, "write_deprecated_outs": True},
}
write_eon_config(min_dir, min_config)

run_eon(min_dir)
print(f"Minimization done: {list(min_dir.glob('minimization.*'))}")
```

### Energy profile

```{code-cell} python
min_profile = plot_dir / "min_profile.png"
run_rgpycrumbs(
    "eon", "plt-min",
    "--job-dir", str(min_dir),
    "--prefix", "minimization",
    "--plot-type", "profile",
    "--dpi", "150",
    "--output", str(min_profile),
)
show_plot(min_profile)
```

### 2D optimization landscape

```{code-cell} python
min_landscape = plot_dir / "min_landscape.png"
run_rgpycrumbs(
    "eon", "plt-min",
    "--job-dir", str(min_dir),
    "--prefix", "minimization",
    "--plot-type", "landscape",
    "--project-path",
    "--surface-type", "grad_matern",
    "--plot-structures", "endpoints",
    "--strip-renderer", "xyzrender",
    "--perspective-tilt", "8",
    "--dpi", "150",
    "--output", str(min_landscape),
)
show_plot(min_landscape)
```

### Convergence panel

```{code-cell} python
min_convergence = plot_dir / "min_convergence.png"
run_rgpycrumbs(
    "eon", "plt-min",
    "--job-dir", str(min_dir),
    "--prefix", "minimization",
    "--plot-type", "convergence",
    "--dpi", "150",
    "--output", str(min_convergence),
)
show_plot(min_convergence)
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
run_eon(neb_dir)
print(f"NEB done: {len(list(neb_dir.glob('neb_*.dat')))} steps")
```

### Energy profile

```{code-cell} python
neb_profile = plot_dir / "neb_profile.png"
run_rgpycrumbs(
    "eon", "plt-neb",
    "--input-dat-pattern", str(neb_dir / "neb_*.dat"),
    "--con-file", str(neb_dir / "neb.con"),
    "--ira-kmax", "14",
    "--force-recompute",
    "--rc-mode", "rmsd",
    "--plot-type", "profile",
    "--plot-structures", "crit_points",
    "--strip-renderer", "xyzrender",
    "--perspective-tilt", "8",
    "--zoom-ratio", "0.15",
    "--dpi", "150",
    "--output-file", str(neb_profile),
)
show_plot(neb_profile)
```

### 2D reaction landscape

```{code-cell} python
neb_landscape = plot_dir / "neb_landscape.png"
run_rgpycrumbs(
    "eon", "plt-neb",
    "--input-dat-pattern", str(neb_dir / "neb_*.dat"),
    "--input-path-pattern", str(neb_dir / "neb_path_*.con"),
    "--con-file", str(neb_dir / "neb.con"),
    "--ira-kmax", "14",
    "--force-recompute",
    "--rc-mode", "rmsd",
    "--plot-type", "landscape",
    "--surface-type", "grad_imq",
    "--show-pts",
    "--landscape-path", "all",
    "--project-path",
    "--plot-structures", "crit_points",
    "--strip-renderer", "xyzrender",
    "--strip-dividers",
    "--perspective-tilt", "8",
    "--zoom-ratio", "0.2",
    "--show-legend",
    "--dpi", "150",
    "--output-file", str(neb_landscape),
)
show_plot(neb_landscape)
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
