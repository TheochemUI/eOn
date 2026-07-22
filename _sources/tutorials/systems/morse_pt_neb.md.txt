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
    "description": "Runnable Morse Pt NEB tutorial with eonclient and modern rgpycrumbs plt-neb 1D/2D plots (built-in potential)."
    "keywords": "eOn Morse Pt, NEB, plt-neb, reaction landscape, documentation system"
---

# Morse Pt NEB (built-in potential)

A compact **Pt surface / adatom** NEB using the built-in `morse_pt` potential.
Geometries match `benchmarks/data/neb_morse_pt/`. The system is small, uses only
a built-in potential, and is suitable for docs CI.

Workflow:

1. write `config.ini` and copy endpoints
2. run `eonclient`
3. plot **1D profiles** and a **2D reaction-valley landscape** with the current
   `rgpycrumbs eon plt-neb` conventions (full history, 1:1 Å panel, structure strip)

## Prerequisites

```{code-block} bash
which eonclient
python -c "import rgpycrumbs, chemparseplot, readcon"
```

## Setup and run

```{code-cell} python
:tags: [remove-cell]

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from IPython.display import Image, display
from rgpycrumbs.eon.helpers import write_eon_config

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("RGPYCRUMBS_AUTO_DEPS", "1")

# myst-nb runs with CWD = this file's directory (tutorials/systems).
DATA = Path("data/morse_pt_neb")
assert (DATA / "reactant.con").is_file(), DATA.resolve()
assert (DATA / "product.con").is_file(), DATA.resolve()

work = Path(tempfile.mkdtemp(prefix="eon_morse_pt_neb_"))
plot_dir = work / "plots"
plot_dir.mkdir()
shutil.copy(DATA / "reactant.con", work / "reactant.con")
shutil.copy(DATA / "product.con", work / "product.con")

cfg = {
    "Main": {"job": "nudged_elastic_band", "random_seed": 42},
    "Potential": {"potential": "morse_pt"},
    "Nudged Elastic Band": {
        "images": 7,
        "initializer": "linear",
        "spring": 5.0,
        "climbing_image_method": "true",
        "climbing_image_converged_only": "true",
    },
    "Optimizer": {
        "opt_method": "lbfgs",
        "converged_force": 0.01,
        "max_iterations": 200,
        "max_move": 0.2,
    },
    "Debug": {"write_movies": True},
}
write_eon_config(work, cfg)


def run_eon(job_dir, *, timeout=300):
    eonclient = str(Path(sys.executable).with_name("eonclient"))
    if not Path(eonclient).exists():
        eonclient = "eonclient"
    r = subprocess.run(
        [eonclient],
        cwd=job_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"eonclient failed ({r.returncode})")
    return r


def run_rgpycrumbs(*args, cwd=None, timeout=600):
    cmd = [sys.executable, "-m", "rgpycrumbs.cli", *args]
    print("$", " ".join(cmd))
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    if r.stdout:
        print(r.stdout)
    if r.stderr:
        print(r.stderr, file=sys.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"rgpycrumbs failed ({r.returncode})")
    return r


def show_plot(path):
    display(Image(filename=str(path)))


print("Workdir:", work)
run_eon(work)
print("neb steps:", len(list(work.glob("neb_*.dat"))))
print("paths:", len(list(work.glob("neb_path_*.con"))))
```

## 1D energy profiles (full history)

Plot every optimization step, highlight the final band, and attach a structure
strip for **all** images on the path (`--plot-structures all`). The saddle is
marked with a gold star on the final profile.

```{code-cell} python
neb_profile = plot_dir / "neb_1d.png"
run_rgpycrumbs(
    "eon",
    "plt-neb",
    "--con-file",
    str(work / "neb.con"),
    "--input-dat-pattern",
    str(work / "neb_*.dat"),
    "--plot-type",
    "profile",
    "--highlight-last",
    "--show-pts",
    "--plot-structures",
    "all",
    "--strip-renderer",
    "xyzrender",
    "--strip-dividers",
    "--xyzrender-config",
    "paton",
    "--rotation",
    "90x,0y,0z",
    "--show-legend",
    "--title",
    "Morse Pt NEB path optimization",
    "--dpi",
    "150",
    "--output-file",
    str(neb_profile),
    cwd=work,
)
show_plot(neb_profile)
```

## 2D reaction-valley landscape

Use the full optimization history for the surface (`--landscape-path all`),
project into progress / orthogonal deviation (`--project-path`), and keep a
**true 1:1 Å** panel (`Δs = Δd`). The strip under the map shows every band image
in two rows when there are many structures.

```{code-cell} python
neb_land = plot_dir / "neb_2d.png"
run_rgpycrumbs(
    "eon",
    "plt-neb",
    "--con-file",
    str(work / "neb.con"),
    "--input-dat-pattern",
    str(work / "neb_*.dat"),
    "--input-path-pattern",
    str(work / "neb_path_*.con"),
    "--plot-type",
    "landscape",
    "--rc-mode",
    "path",
    "--landscape-mode",
    "surface",
    "--landscape-path",
    "all",
    "--surface-type",
    "grad_imq",
    "--project-path",
    "--show-pts",
    "--highlight-last",
    "--plot-structures",
    "all",
    "--strip-renderer",
    "xyzrender",
    "--strip-dividers",
    "--xyzrender-config",
    "paton",
    "--rotation",
    "90x,0y,0z",
    "--show-legend",
    "--title",
    "Morse Pt NEB-RMSD surface",
    "--dpi",
    "150",
    "--output-file",
    str(neb_land),
    cwd=work,
)
show_plot(neb_land)
```

## See also

- {doc}`../visualization` — PET-MAD / vinyl alcohol visualization walkthrough
- {doc}`/user_guide/neb` — NEB parameters
- [atomistic-cookbook eon-pet-neb](https://atomistic-cookbook.org/examples/eon-pet-neb/eon-pet-neb.html)
  — metatomic oxadiazole NEB with the same plot stack
