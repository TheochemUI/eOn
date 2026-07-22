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
    "description": "Runnable Morse Pt dimer saddle search with eonclient and rgpycrumbs plt-saddle visualization."
    "keywords": "eOn saddle search, dimer, Morse Pt, plt-saddle"
---

# Morse Pt saddle search (built-in potential)

A small **dimer** saddle search on `morse_pt` using the geometries from
`benchmarks/data/one_pt_saddle_search/`. After the client finishes, plot the
trajectory with `rgpycrumbs eon plt-saddle` (same single-ended landscape stack
as minimization).

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

DATA = Path("data/pt_saddle")
assert (DATA / "pos.con").is_file(), DATA.resolve()

work = Path(tempfile.mkdtemp(prefix="eon_pt_saddle_"))
plot_dir = work / "plots"
plot_dir.mkdir()
shutil.copy(DATA / "pos.con", work / "pos.con")
if (DATA / "displacement.con").is_file():
    shutil.copy(DATA / "displacement.con", work / "displacement.con")

cfg = {
    "Main": {"job": "saddle_search", "temperature": 300, "random_seed": 42},
    "Potential": {"potential": "morse_pt"},
    "Optimizer": {
        "converged_force": 0.01,
        "max_iterations": 500,
        "opt_method": "lbfgs",
        "max_move": 0.2,
    },
    "Saddle Search": {
        "min_mode_method": "dimer",
        "displace_magnitude": 0.1,
        "displace_radius": 3.3,
        "max_energy": 10.0,
    },
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
print("outputs:", sorted(p.name for p in work.iterdir() if p.suffix in {".con", ".dat"}))
```

## Profile and landscape

```{code-cell} python
# Profile (energy / eigenvalue vs iteration when available)
out_p = plot_dir / "saddle_profile.png"
try:
    run_rgpycrumbs(
        "eon",
        "plt-saddle",
        "--job-dir",
        str(work),
        "--label",
        "pt saddle",
        "--plot-type",
        "profile",
        "--dpi",
        "150",
        "-o",
        str(out_p),
    )
    show_plot(out_p)
except Exception as exc:
    print("plt-saddle profile skipped:", exc)
```

```{code-cell} python
out_l = plot_dir / "saddle_2d.png"
try:
    run_rgpycrumbs(
        "eon",
        "plt-saddle",
        "--job-dir",
        str(work),
        "--label",
        "pt saddle",
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
        "--dpi",
        "150",
        "-o",
        str(out_l),
    )
    show_plot(out_l)
except Exception as exc:
    print("plt-saddle landscape skipped:", exc)
```

## See also

- {doc}`/user_guide/saddle_search` and {doc}`/user_guide/dimer`
- {doc}`lj_minimization` for the shared single-ended plot conventions
