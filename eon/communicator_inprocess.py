"""In-process communicator: jobs run as pyeonclient.Matter (no eonclient binary).

This is the Matter-through path:
  Structure/Atoms  →  pyeonclient.Matter  →  relax / (future jobs)  →  Structure

Requires the pyeonclient extension (``-Dwith_pyeonclient=true``). Cluster/MPI
still use the file-based communicators.
"""

from __future__ import annotations

import logging
import os
from io import StringIO
from typing import Any

import numpy as np

from eon.communicator import Communicator, CommunicatorError

logger = logging.getLogger("communicator")


def _require_pyeonclient():
    try:
        import pyeonclient as pc
    except ImportError as e:
        raise CommunicatorError(
            "inprocess communicator needs pyeonclient "
            "(build with -Dwith_pyeonclient=true)"
        ) from e
    return pc


def _params_from_invariants(pc, invariants: dict) -> Any:
    """Load Parameters from config.ini bytes if present, else defaults."""
    params = pc.Parameters()
    # invariants values are (StringIO, mode) or StringIO
    for name, val in invariants.items():
        base = os.path.basename(name)
        if base not in ("config.ini", "config"):
            continue
        content = val[0] if isinstance(val, tuple) else val
        if hasattr(content, "getvalue"):
            text = content.getvalue()
        elif hasattr(content, "read"):
            text = content.read()
        else:
            text = str(content)
        # write temp for Parameters.load
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".ini", delete=False) as fh:
            fh.write(text)
            path = fh.name
        try:
            params.load(path)
        finally:
            os.unlink(path)
        return params
    params.potential = pc.PotType.LJ
    params.quiet = True
    params.write_log = False
    return params


def _structure_from_job_con(job: dict, key: str = "pos.con"):
    from eon import fileio as io

    blob = job.get(key)
    if blob is None:
        raise CommunicatorError(f"job missing {key}")
    if hasattr(blob, "getvalue"):
        text = blob.getvalue()
    elif hasattr(blob, "read"):
        if hasattr(blob, "seek"):
            blob.seek(0)
        text = blob.read()
    else:
        text = str(blob)
    return io.loadcon(StringIO(text))


def _results_dat(status: int, energy: float, force_calls: int, job_type: str) -> str:
    return (
        f"{status} termination_reason\n"
        f"{'GOOD' if status == 0 else 'FAIL'} termination_reason_text\n"
        f"{job_type} job_type\n"
        f"{energy:.12e} potential_energy\n"
        f"{force_calls} total_force_calls\n"
    )


class LocalInProcess(Communicator):
    """Run client work in-process via Matter (nanobind), not a subprocess."""

    def __init__(self, scratchpath, bundle_size=1, config=None):
        if config is None:
            raise TypeError("LocalInProcess requires a ConfigClass instance")
        Communicator.__init__(self, scratchpath, bundle_size, config=config)
        self._pc = _require_pyeonclient()
        self._finished: list[dict] = []

    def get_queue_size(self):
        return 0

    def get_number_in_progress(self):
        return 0

    def cancel_state(self, state):
        return 0

    def submit_jobs(self, data, invariants):
        """Run each job dict in-process. Supports minimization-like jobs.

        Job dict keys (legacy file names kept for explorer compatibility):
          * pos.con — reactant structure (StringIO of .con text)
          * id — job id string
        """
        pc = self._pc
        params = _params_from_invariants(pc, invariants)
        pot = pc.make_potential(params)

        for job in data:
            jid = job.get("id", "job")
            try:
                structure = _structure_from_job_con(job, "pos.con")
            except Exception as e:
                logger.exception("inprocess: failed to parse pos.con for %s", jid)
                raise CommunicatorError(str(e)) from e

            from pyeonclient.bridge import structure_to_matter, matter_to_structure

            matter = structure_to_matter(structure, pot, params)
            # Default job for in-process path: minimize (Matter.relax).
            # Full JobType dispatch lands as more C++ entry points are bound.
            # relax returns (Matter, converged: bool); default is non-inplace.
            matter, converged = matter.relax(
                inplace=True, quiet=True, write_movie=False, checkpoint=False
            )
            out = matter_to_structure(matter)

            import eon.fileio as fio

            min_io = StringIO()
            fio.savecon(min_io, out)
            min_io.seek(0)

            energy = float(matter.potential_energy)
            fcalls = int(matter.force_calls)
            status = 0 if converged else 1
            results = StringIO(
                _results_dat(status, energy, fcalls, "minimization")
            )

            self._finished.append(
                {
                    "id": jid,
                    # Bundle index within the job, as the file-based
                    # communicators set it; one task per job here.
                    "number": 0,
                    "name": str(jid),
                    "min.con": min_io,
                    "results.dat": results,
                    # also expose Matter for callers that want zero re-parse
                    "_matter": matter,
                    "_structure": out,
                    "_energy": energy,
                    "_converged": converged,
                }
            )
            logger.info(
                "inprocess job %s: converged=%s E=%.6f fcalls=%s",
                jid,
                converged,
                energy,
                fcalls,
            )

    def get_results(self, resultspath=None, keep_result=None):
        """Return finished job dicts (legacy StringIO + Matter fields)."""
        out = self._finished
        self._finished = []
        # Compatibility: some callers expect a list of lists from unbundle
        if not out:
            return []
        return out
