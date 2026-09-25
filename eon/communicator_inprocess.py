"""In-process communicator: jobs run as pyeonclient.Matter (no eonclient binary).

Geometry crosses the job dict as a Structure (numpy working set) or a
readcon.ConFrame. Matter is the potential-bearing client object. Result
saddle and product values are ConFrames. This path does not read or write
.con text.
"""

from __future__ import annotations

import logging
import os
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np

from eon.communicator import Communicator, CommunicatorError

logger = logging.getLogger("communicator")

_GEOMETRY_KEYS = ("structure", "conframe", "reactant", "pos", "pos.con")


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
        base = Path(name).name
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


def _is_structure(obj) -> bool:
    from eon.structure import Structure

    return isinstance(obj, Structure)


def _is_conframe(obj) -> bool:
    cls = type(obj)
    return cls.__name__ == "ConFrame" and "readcon" in cls.__module__


def _is_con_text(obj) -> bool:
    if _is_structure(obj) or _is_conframe(obj):
        return False
    return (
        hasattr(obj, "getvalue")
        or hasattr(obj, "read")
        or isinstance(obj, (str, bytes))
    )


def _structure_from_geometry(blob, where: str):
    """Structure or ConFrame to the numpy working set. Rejects .con text."""
    if _is_con_text(blob):
        raise CommunicatorError(
            f"{where} must be a Structure or ConFrame, not .con text"
        )
    if _is_structure(blob):
        return blob
    if _is_conframe(blob):
        from eon.structure import Structure

        return Structure.from_conframe(blob)
    raise CommunicatorError(
        f"{where} must be a Structure or ConFrame, got {type(blob).__name__}"
    )


def _structure_from_job(job: dict, invariants: dict | None):
    for key in _GEOMETRY_KEYS:
        if key in job and job[key] is not None:
            return _structure_from_geometry(job[key], f"job[{key!r}]")
    if invariants:
        for key in _GEOMETRY_KEYS:
            if key not in invariants or invariants[key] is None:
                continue
            val = invariants[key]
            blob = val[0] if isinstance(val, tuple) else val
            if _is_structure(blob) or _is_conframe(blob):
                return _structure_from_geometry(blob, f"invariants[{key!r}]")
    raise CommunicatorError(
        "inprocess job needs a Structure or ConFrame "
        "(structure, conframe, reactant, pos, or pos.con)"
    )


def _conframe_of(structure):
    return structure.to_conframe()


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
        """Run each job dict in-process.

        Geometry is a Structure or a readcon.ConFrame on one of
        ``structure``, ``conframe``, ``reactant``, ``pos``, or ``pos.con``,
        or the same objects on the invariant dict. ``.con`` text is refused.
        The result carries ``product`` and, for a saddle search, ``saddle``
        as ConFrames. Nothing is written to disk.
        """
        pc = self._pc
        params = _params_from_invariants(pc, invariants)
        pot = pc.make_potential(params)
        job_kind = getattr(params, "job", pc.JobType.Minimization)

        from pyeonclient.bridge import structure_to_matter, matter_to_structure

        for job in data:
            jid = job.get("id", "job")
            try:
                structure = _structure_from_job(job, invariants)
            except CommunicatorError:
                logger.exception("inprocess: missing geometry for %s", jid)
                raise
            except Exception as e:
                logger.exception("inprocess: failed to read geometry for %s", jid)
                raise CommunicatorError(str(e)) from e

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
                "product": _conframe_of(out),
                "results.dat": results,
                "_matter": matter,
                "_structure": out,
                "_energy": energy,
                "_converged": bool(payload.get("converged", status == 0)),
            }
            saddle = payload.get("saddle")
            if saddle is not None:
                rec["saddle"] = _conframe_of(matter_to_structure(saddle))
                rec["_saddle"] = saddle
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
        """Return finished job dicts. Geometry is ConFrame, not a file."""
        out = self._finished
        self._finished = []
        if not out:
            return []
        return out
