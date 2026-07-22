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
    "description": "Runnable LJ cluster minimization with eonclient and modern rgpycrumbs plt-min profile, landscape, and convergence plots."
    "keywords": "eOn LJ cluster, minimization, plt-min, optimization landscape"
---

# LJ cluster minimization (built-in potential)

Minimize a small **Lennard-Jones cluster** (`ljcluster`) from
`benchmarks/data/min_lj_cluster/`. Movies are required for the modern plot
stack.

Plotting conventions:

- energy **profile** and **convergence** from the minimization trajectory
- a **single** 2D landscape for this job (one endpoint → one RMSD frame), with
  title from `--label` and **initial** / **minimized** captions
- relative energy on the landscape colorbar for readable ticks

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

DATA = Path("data/lj_min")
assert (DATA / "pos.con").is_file(), DATA.resolve()

work = Path(tempfile.mkdtemp(prefix="eon_lj_min_"))
plot_dir = work / "plots"
plot_dir.mkdir()
shutil.copy(DATA / "pos.con", work / "pos.con")

cfg = {
    "Main": {"job": "minimization", "random_seed": 42},
    "Potential": {"potential": "ljcluster"},
    "Optimizer": {
        "opt_method": "lbfgs",
        "converged_force": 0.01,
        "max_iterations": 500,
        "max_move": 0.2,
    },
    # Movies for plt-min landscape; optional legacy .dat for older tools.
    "Debug": {"write_movies": True, "write_deprecated_outs": True},
}
write_eon_config(work, cfg)


def run_eon(job_dir, *, timeout=300):
    eonclient = str(Path(sys.executable).with_name("eonclient"))
    if not Path(eonclient).exists():
        eonclient = "eonclient"
    r = subprocess.run(
        [eonclient], cwd=job_dir, capture_output=True, text=True, timeout=timeout
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
print("outputs:", sorted(p.name for p in work.glob("minimization.*")))
```

## Energy profile

```{code-cell} python
out = plot_dir / "min_1d.png"
run_rgpycrumbs(
    "eon",
    "plt-min",
    "--job-dir",
    str(work),
    "--label",
    "lj cluster",
    "--plot-type",
    "profile",
    "--dpi",
    "150",
    "-o",
    str(out),
)
show_plot(out)
```

## 2D optimization landscape

One landscape for this trajectory only. The title becomes **Lj cluster
minimization**; endpoints are labeled **initial** / **minimized**.

```{code-cell} python
out = plot_dir / "min_2d.png"
run_rgpycrumbs(
    "eon",
    "plt-min",
    "--job-dir",
    str(work),
    "--label",
    "lj cluster",
    "--plot-type",
    "landscape",
    "--project-path",
    "--surface-type",
    "grad_imq",
    "--plot-structures",
    "endpoints",
    "--strip-renderer",
    "xyzrender",
    "--strip-dividers",
    "--xyzrender-config",
    "paton",
    "--rotation",
    "90x,0y,0z",
    "--dpi",
    "150",
    "-o",
    str(out),
)
show_plot(out)
```

## Convergence

```{code-cell} python
out = plot_dir / "min_conv.png"
run_rgpycrumbs(
    "eon",
    "plt-min",
    "--job-dir",
    str(work),
    "--label",
    "lj cluster",
    "--plot-type",
    "convergence",
    "--dpi",
    "150",
    "-o",
    str(out),
)
show_plot(out)
```

## Multiple endpoints (reactant *and* product)

When both ends of a path are minimized, plot **two** landscapes (one job dir
each). Overlay only the 1D profile / convergence panels:

```{code-block} bash
# separate 2D frames
python -m rgpycrumbs.cli eon plt-min --job-dir min_reactant --label reactant \
  --plot-type landscape --project-path --plot-structures endpoints -o min_2D_reactant.png
python -m rgpycrumbs.cli eon plt-min --job-dir min_product --label product \
  --plot-type landscape --project-path --plot-structures endpoints -o min_2D_product.png

# shared 1D overlays
python -m rgpycrumbs.cli eon plt-min \
  --job-dir min_reactant --label reactant \
  --job-dir min_product --label product \
  --plot-type profile -o min_1D.png
```

## See also

- {doc}`morse_pt_neb` — built-in NEB + `plt-neb`
- {doc}`/user_guide/minimization`
- {doc}`../visualization` — PET-MAD vinyl alcohol visualization
