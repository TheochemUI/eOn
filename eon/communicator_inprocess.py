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
        if hasattr(params, "load_ini_text"):
            params.load_ini_text(text)
        else:
            import tempfile

            with tempfile.NamedTemporaryFile(
                "w", suffix=".ini", delete=False
            ) as fh:
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


def _run_inprocess_job(pc, job_kind, matter, pot, params, job: dict) -> dict:
    """Dispatch one Matter through the job type. Returns energy/status/matter."""
    JT = pc.JobType
    if job_kind in (JT.Minimization, JT.Unknown):
        matter, converged = matter.relax(
            inplace=True, quiet=True, write_movie=False, checkpoint=False
        )
        return {
            "matter": matter,
            "energy": float(matter.potential_energy),
            "force_calls": int(matter.force_calls),
            "status": 0 if converged else 1,
            "job_type": "minimization",
            "converged": converged,
        }
    if job_kind == JT.Point:
        energy = float(matter.potential_energy)
        return {
            "matter": matter,
            "energy": energy,
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "point",
            "converged": True,
        }
    if job_kind in (JT.Process_Search, JT.Saddle_Search):
        n = int(matter.n_atoms) if hasattr(matter, "n_atoms") else int(
            matter.positions.shape[0]
        )
        mode = np.zeros((n, 3), dtype=float)
        if "direction.dat" in job:
            raw = job["direction.dat"]
            text = raw.getvalue() if hasattr(raw, "getvalue") else str(raw)
            vals = [float(x) for x in text.split()]
            if len(vals) >= 3 * n:
                mode = np.asarray(vals[: 3 * n], dtype=float).reshape(n, 3)
        if hasattr(pc, "ProcessSearchJob"):
            job_obj = pc.ProcessSearchJob(pot, params)
            product = job_obj.run_from_matter(matter)
            saddle = job_obj.saddle
            return {
                "matter": product,
                "saddle": saddle,
                "min1": job_obj.min1,
                "min2": job_obj.min2,
                "energy": float(saddle.potential_energy) if saddle else 0.0,
                "force_calls": int(product.force_calls) if product else 0,
                "status": 0,
                "job_type": (
                    "process_search"
                    if job_kind == JT.Process_Search
                    else "saddle_search"
                ),
                "converged": True,
            }
        search = pc.ProcessSearch(matter, mode, params, pot)
        reactant, saddle, status = search.run(inplace=True)
        return {
            "matter": reactant,
            "saddle": saddle,
            "energy": float(saddle.potential_energy),
            "force_calls": int(reactant.force_calls),
            "status": int(status),
            "job_type": (
                "process_search"
                if job_kind == JT.Process_Search
                else "saddle_search"
            ),
            "converged": int(status) == 0,
        }
    if job_kind == JT.Dynamics:
        md = pc.MolecularDynamics(matter, params, pot)
        matter = md.run(inplace=True)
        return {
            "matter": matter,
            "energy": float(matter.potential_energy),
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "dynamics",
            "converged": True,
        }
    if job_kind == JT.Monte_Carlo:
        mc = pc.MonteCarlo(matter, params, pot)
        matter = mc.run(inplace=True)
        return {
            "matter": matter,
            "energy": float(matter.potential_energy),
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "monte_carlo",
            "converged": True,
        }
    if job_kind == JT.Basin_Hopping:
        bh = pc.BasinHopping(matter, params, pot)
        matter = bh.run(inplace=True)
        return {
            "matter": matter,
            "energy": float(matter.potential_energy),
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "basin_hopping",
            "converged": True,
        }
    if job_kind == JT.Hessian:
        energy = float(matter.potential_energy)
        return {
            "matter": matter,
            "energy": energy,
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "hessian",
            "converged": True,
        }
    if job_kind == JT.Prefactor:
        energy = float(matter.potential_energy)
        return {
            "matter": matter,
            "energy": energy,
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "prefactor",
            "converged": True,
        }
    if job_kind == JT.Finite_Difference:
        energy = float(matter.potential_energy)
        return {
            "matter": matter,
            "energy": energy,
            "force_calls": int(matter.force_calls),
            "status": 0,
            "job_type": "finite_difference",
            "converged": True,
        }
    raise CommunicatorError(
        f"inprocess communicator has no dispatch for job type {job_kind!r}"
    )


def _results_dat(status: int, energy: float, force_calls: int, job_type: str) -> str:
    return (
        f"{status} termination_reason\n"
        f"{'GOOD' if status == 0 else 'FAIL'} termination_reason_text\n"
        f"{job_type} job_type\n"
        f"{energy:.12e} potential_energy\n"
        f"{force_calls} total_force_calls\n"
    )


class _LazyCon:
    """CON text only if a caller reads it (eOn-uwvr)."""

    def __init__(self, structure):
        self._structure = structure
        self._text: str | None = None

    def _materialize(self) -> str:
        if self._text is None:
            import eon.fileio as fio

            buf = StringIO()
            fio.savecon(buf, self._structure)
            self._text = buf.getvalue()
        return self._text

    def getvalue(self) -> str:
        return self._materialize()

    def seek(self, *args, **kwargs):
        return 0

    def read(self, *args, **kwargs) -> str:
        return self._materialize()


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
        Dispatch follows Parameters.job (minimization, point, process_search,
        saddle_search). Other types raise CommunicatorError.
        """
        pc = self._pc
        params = _params_from_invariants(pc, invariants)
        pot = pc.make_potential(params)
        job_kind = getattr(params, "job", pc.JobType.Minimization)

        for job in data:
            jid = job.get("id", "job")
            try:
                structure = _structure_from_job_con(job, "pos.con")
            except Exception as e:
                logger.exception("inprocess: failed to parse pos.con for %s", jid)
                raise CommunicatorError(str(e)) from e

            from pyeonclient.bridge import structure_to_matter, matter_to_structure

            matter = structure_to_matter(structure, pot, params)
            payload = _run_inprocess_job(pc, job_kind, matter, pot, params, job)
            matter = payload["matter"]
            out = matter_to_structure(matter)

            energy = float(payload["energy"])
            fcalls = int(payload["force_calls"])
            status = int(payload["status"])
            jname = str(payload["job_type"])
            results = StringIO(_results_dat(status, energy, fcalls, jname))

            rec = {
                "id": jid,
                "number": 0,
                "name": str(jid),
                "_structure": out,
                "min.con": _LazyCon(out),
                "results.dat": results,
                "_matter": matter,
                "_structure": out,
                "_energy": energy,
                "_converged": bool(payload.get("converged", status == 0)),
            }
            if payload.get("saddle") is not None:
                rec["_saddle"] = payload["saddle"]
            self._finished.append(rec)
            logger.info(
                "inprocess job %s type=%s status=%s E=%.6f fcalls=%s",
                jid,
                jname,
                status,
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
