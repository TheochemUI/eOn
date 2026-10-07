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

from eon.cancel import CancelToken
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
    # readcon's pyo3 class reports its module as "builtins", so match the
    # class itself rather than a module name.
    try:
        from readcon import ConFrame
    except ImportError:
        return False
    return isinstance(obj, ConFrame)


def _is_con_text(obj) -> bool:
    if _is_structure(obj) or _is_conframe(obj):
        return False
    return (
        hasattr(obj, "getvalue")
        or hasattr(obj, "read")
        or isinstance(obj, (str, bytes))
    )


def _con_text_of(blob) -> str:
    if hasattr(blob, "getvalue"):
        return blob.getvalue()
    if hasattr(blob, "read"):
        if hasattr(blob, "seek"):
            blob.seek(0)
        return blob.read()
    if isinstance(blob, bytes):
        return blob.decode()
    return str(blob)


def _structure_from_geometry(blob, where: str):
    """Structure, ConFrame or .con text to the numpy working set.

    The AKMC drivers (explorer, basin hopping, escape rate, parallel replica)
    still send ``pos.con`` as a StringIO of .con text, so text is parsed here
    rather than refused.
    """
    if _is_con_text(blob):
        from eon import fileio as io

        try:
            return io.loadcon(StringIO(_con_text_of(blob)))
        except Exception as e:
            raise CommunicatorError(f"{where}: could not read .con text ({e})") from e
    if _is_structure(blob):
        return blob
    if _is_conframe(blob):
        from eon.structure import Structure

        return Structure.from_conframe(blob)
    raise CommunicatorError(
        f"{where} must be a Structure, ConFrame or .con text, "
        f"got {type(blob).__name__}"
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
            if blob is not None:
                return _structure_from_geometry(blob, f"invariants[{key!r}]")
    raise CommunicatorError(
        "inprocess job needs a Structure, ConFrame or .con text "
        "(structure, conframe, reactant, pos, or pos.con)"
    )


def _conframe_of(structure):
    return structure.to_conframe()


class _LazyCon:
    """.con text for callers that read min.con or saddle.con, built on first read."""

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


def _run_inprocess_job(pc, job_kind, matter, pot, params, job: dict, token=None) -> dict:
    """Dispatch one Matter through the job type. Returns energy/status/matter."""
    if token is not None:
        token.raise_if_cancelled()
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


def _job_result(
    status: int,
    energy: float,
    force_calls: int,
    job_type: str,
    *,
    cancelled: bool = False,
) -> dict:
    """Typed in-process result. Cluster adapters format text from this dict."""
    if cancelled:
        reason = "cancelled"
    elif status == 0:
        reason = "GOOD"
    else:
        reason = "FAIL"
    return {
        "termination_reason": status,
        "termination_reason_text": reason,
        "job_type": job_type,
        "potential_energy": energy,
        "total_force_calls": force_calls,
    }


def _results_dat(
    status: int,
    energy: float,
    force_calls: int,
    job_type: str,
    *,
    cancelled: bool = False,
) -> str:
    data = _job_result(
        status, energy, force_calls, job_type, cancelled=cancelled
    )
    lines = []
    for key, val in data.items():
        if isinstance(val, float):
            lines.append(f"{val:.12e} {key}")
        else:
            lines.append(f"{val} {key}")
    return "\n".join(lines) + "\n"


class LocalInProcess(Communicator):
    """Run client work in-process via Matter (nanobind), not a subprocess."""

    def __init__(self, scratchpath, bundle_size=1, config=None):
        if config is None:
            raise TypeError("LocalInProcess requires a ConfigClass instance")
        Communicator.__init__(self, scratchpath, bundle_size, config=config)
        self._pc = _require_pyeonclient()
        self._finished: list[dict] = []
        self.token = CancelToken()
        self._in_submit = False
        self._cpp_token = None

    def get_queue_size(self):
        return 0

    def get_number_in_progress(self):
        return 0

    def cancel_state(self, state):
        """Ask the running batch to stop. Idle calls return 0.

        The Python token stops the next dispatch. The C++ token is what
        Matter polls inside a compiled relax, NEB, or saddle search.
        """
        del state
        if not self._in_submit:
            return 0
        self.token.cancel()
        cpp = self._cpp_token
        if cpp is not None and not cpp.requested():
            cpp.request()
        return 1

    def submit_jobs(self, data, invariants):
        """Run each job dict in-process.

        Geometry is a Structure, a readcon.ConFrame or .con text on one of
        ``structure``, ``conframe``, ``reactant``, ``pos``, or ``pos.con``,
        or the same objects on the invariant dict. The result carries
        ``product`` and, for a saddle search, ``saddle`` as ConFrames, and
        the same geometry as ``min.con`` and ``saddle.con`` for the drivers
        that read .con text. Nothing is written to disk.
        """
        pc = self._pc
        params = _params_from_invariants(pc, invariants)
        pot = pc.make_potential(params)
        job_kind = getattr(params, "job", pc.JobType.Minimization)

        self._in_submit = True
        try:
            self._submit_loop(data, invariants, pc, pot, params, job_kind)
        finally:
            self._in_submit = False
            self.token.reset()

    def _submit_loop(self, data, invariants, pc, pot, params, job_kind):
        from pyeonclient.bridge import structure_to_matter, matter_to_structure

        for job in data:
            self.token.raise_if_cancelled()
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
            if hasattr(pc, "CancelToken") and hasattr(matter, "set_cancel_token"):
                self._cpp_token = pc.CancelToken()
                matter.set_cancel_token(self._cpp_token)
                if self.token.cancelled:
                    self._cpp_token.request()
            cancelled_exc = getattr(pc, "JobCancelled", None)
            try:
                payload = _run_inprocess_job(
                    pc, job_kind, matter, pot, params, job, token=self.token
                )
            except Exception as exc:
                if cancelled_exc is None or not isinstance(exc, cancelled_exc):
                    raise
                payload = {
                    "matter": matter,
                    "energy": 0.0,
                    "force_calls": int(getattr(matter, "force_calls", 0) or 0),
                    "status": 1,
                    "job_type": str(getattr(job_kind, "name", job_kind)),
                    "converged": False,
                    "cancelled": True,
                }
            finally:
                self._cpp_token = None
            matter = payload["matter"]
            out = matter_to_structure(matter)

            energy = float(payload["energy"])
            fcalls = int(payload["force_calls"])
            status = int(payload["status"])
            jname = str(payload["job_type"])
            cancelled = bool(payload.get("cancelled"))
            job_result = _job_result(
                status, energy, fcalls, jname, cancelled=cancelled
            )

            rec = {
                "id": jid,
                "number": 0,
                "name": str(jid),
                "product": _conframe_of(out),
                "min.con": _LazyCon(out),
                "job_result": job_result,
                "_matter": matter,
                "_structure": out,
                "_energy": energy,
                "_converged": bool(payload.get("converged", status == 0)),
            }
            saddle = payload.get("saddle")
            if saddle is not None:
                saddle_structure = matter_to_structure(saddle)
                rec["saddle"] = _conframe_of(saddle_structure)
                rec["saddle.con"] = _LazyCon(saddle_structure)
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
