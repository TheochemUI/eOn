#!/usr/bin/env python
"""CPython AKMC / server dispatcher.

Entry points:
  python -m eon.server
  python -m eon
  eon-server   (console script → eon.server:main)
"""

from __future__ import annotations

import shutil
from io import StringIO
from pathlib import Path
from typing import Callable, Optional

from eon.config import ConfigClass


def select_job_runner(job: str) -> Optional[Callable[[ConfigClass], None]]:
    """Return the registered server runner for *job*, or None for fallback.

    Job modules are imported lazily so ``import eon.server`` stays light and
    does not require the full AKMC stack until a job is actually selected.
    """
    key = job.strip().lower()
    if key == "akmc":
        from eon import akmc

        return akmc.main
    if key in ("parallel_replica", "unbiased_parallel_replica"):
        from eon import parallelreplica

        return parallelreplica.main
    if key == "basin_hopping":
        from eon import basinhopping

        return basinhopping.main
    if key == "escape_rate":
        from eon import escaperate

        return escaperate.main
    return None


def _warn_pos_con_in_potfiles(config: ConfigClass) -> None:
    pot = Path(config.path_pot)
    fnames = [p.name for p in pot.iterdir()] if pot.is_dir() else []
    if "pos.con" in fnames:
        print(
            "WARNING: pos.con found in potfiles path. Are you sure you want this? "
            "It will overwrite the pos.con in the calculation directory when your "
            "jobs are being run."
        )


def _fallback_single_job(config: ConfigClass) -> None:
    """Submit the current working directory as one client job via communicator."""
    from eon import communicator
    from eon import fileio as io

    config.path_scratch = config.path_root
    comm = communicator.get_communicator(config)

    job: dict = {}
    for f in Path(".").iterdir():
        if not f.is_file():
            continue
        parts = f.name.split(".")
        if len(parts) > 1:
            f_passed = parts[0] + "." + parts[1]
            job[f_passed] = StringIO(f.read_text())
    job["id"] = "output"
    if Path("output_old").is_dir():
        shutil.rmtree("output_old")
    if Path("output").is_dir():
        shutil.move("output", "output_old")
    comm.submit_jobs([job], {})


def server(config: ConfigClass | None = None) -> None:
    """Initialize config and dispatch to the configured main job.

    *config* defaults to the process-edge module singleton for CLI
    compatibility; callers should prefer constructing ``ConfigClass()`` and
    passing it explicitly.
    """
    if config is None:
        config = ConfigClass()
    config.init()
    _warn_pos_con_in_potfiles(config)

    job = config.main_job.lower()
    runner = select_job_runner(job)
    if runner is not None:
        runner(config)
    else:
        _fallback_single_job(config)


def main() -> None:
    """Console-script entry point for ``eon-server``."""
    server()


if __name__ == "__main__":
    main()
