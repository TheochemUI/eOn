#!/usr/bin/env python3
"""Regenerate the LJ tutorial figures on a docs build.

Minimizes ``docs/lj13.con`` with the built-in Lennard-Jones potential
(``write_movies`` and the legacy minimization sidecar), then asks
``rgpycrumbs eon plt-min`` for the profile and landscape PNGs under
``docs/source/fig/generated/``. The directory is gitignored; docs CI
treats a missing or empty PNG as a failed pipeline.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_POS = _REPO / "docs" / "lj13.con"
_OUT = _REPO / "docs" / "source" / "fig" / "generated"
_MIN_PNG_BYTES = 1000

_CONFIG = """\
[Main]
job = minimization
random_seed = 42

[Potential]
potential = lj

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 200
max_move = 0.2

[Debug]
write_movies = true
write_deprecated_outs = true
"""


def _eonclient() -> str:
    sibling = Path(sys.executable).with_name("eonclient")
    if sibling.is_file():
        return str(sibling)
    found = shutil.which("eonclient")
    if not found:
        raise SystemExit("eonclient is not on PATH")
    return found


def _run(cmd: list[str], cwd: Path, timeout: int) -> None:
    print("$", " ".join(cmd), flush=True)
    result = subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        timeout=timeout,
        check=False,
        capture_output=True,
    )
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.returncode != 0:
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        raise SystemExit(f"command failed ({result.returncode}): {cmd[0]}")


def _require_png(path: Path) -> None:
    if not path.is_file() or path.stat().st_size < _MIN_PNG_BYTES:
        raise SystemExit(f"tutorial figure missing or empty: {path}")
    print(f"wrote {path} ({path.stat().st_size} bytes)")


def main() -> None:
    if not _POS.is_file():
        raise SystemExit(f"missing geometry: {_POS}")
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("RGPYCRUMBS_AUTO_DEPS", "1")
    _OUT.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="eon_lj13_min_") as tmp:
        work = Path(tmp)
        shutil.copy(_POS, work / "pos.con")
        (work / "config.ini").write_text(_CONFIG)
        _run([_eonclient()], work, timeout=300)
        for name in ("minimization.con", "minimization.dat"):
            if not (work / name).is_file():
                raise SystemExit(f"eonclient did not write {name}")
        plots = (
            ("profile", _OUT / "lj13_min_profile.png", []),
            (
                "landscape",
                _OUT / "lj13_min_landscape.png",
                [
                    "--project-path",
                    "--surface-type",
                    "grad_imq",
                    "--plot-structures",
                    "endpoints",
                ],
            ),
        )
        for plot_type, dest, extra in plots:
            _run(
                [
                    sys.executable,
                    "-m",
                    "rgpycrumbs.cli",
                    "eon",
                    "plt-min",
                    "--job-dir",
                    str(work),
                    "--label",
                    "LJ13",
                    "--prefix",
                    "minimization",
                    "--plot-type",
                    plot_type,
                    "--dpi",
                    "150",
                    *extra,
                    "-o",
                    str(dest),
                ],
                work,
                timeout=600,
            )
            _require_png(dest)


if __name__ == "__main__":
    main()
