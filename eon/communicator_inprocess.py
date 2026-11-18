"""In-process communicator: jobs run as pyeonclient.Matter (no eonclient binary).

Geometry crosses the job dict as a Structure (numpy working set) or a
readcon.ConFrame. Matter is the potential-bearing client object. Result
saddle and product values are ConFrames. This path does not read or write
.con text.
"""

from __future__ import annotations

import logging
import os
import hashlib
import re
import threading
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np

from eon.cancel import CancelToken
from eon.communicator import Communicator, CommunicatorError

logger = logging.getLogger("communicator")

_GEOMETRY_KEYS = ("structure", "conframe", "reactant", "pos", "pos.con")
_PROCESS_CWD_LOCK = threading.Lock()


def _require_pyeonclient():
    try:
        import pyeonclient as pc
    except ImportError as e:
        raise CommunicatorError(
            "inprocess communicator needs pyeonclient "
            "(build with -Dwith_pyeonclient=true)"
        ) from e
    return pc


def _config_text(invariants: dict) -> str | None:
    for name, val in invariants.items():
        if Path(name).name not in ("config.ini", "config"):
            continue
        content = val[0] if isinstance(val, tuple) else val
        if hasattr(content, "getvalue"):
            return content.getvalue()
        if hasattr(content, "read"):
            if hasattr(content, "seek"):
                content.seek(0)
            return content.read()
        return str(content)
    return None


def _session_config_fingerprint(text: str | None) -> str:
    stable = re.sub(
        r"(?mi)^\s*random_seed\s*=.*$",
        "",
        text or "",
    )
    return hashlib.sha256(stable.encode("utf-8")).hexdigest()


def _params_from_invariants(pc, invariants: dict) -> Any:
    """Load Parameters from config.ini bytes if present, else defaults."""
    params = pc.Parameters()
    text = _config_text(invariants)
    if text is not None:
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


def _geometry_fingerprint(structure) -> str:
    digest = hashlib.sha256()
    for value in (structure.r, structure.box, structure.mass, structure.free):
        array = np.ascontiguousarray(value)
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.tobytes())
    digest.update("\0".join(structure.names).encode("utf-8"))
    return digest.hexdigest()


def _state_key(job: dict) -> str:
    if "state" in job:
        return str(job["state"])
    return str(job.get("id", "state")).split("_", 1)[0]


def _con_io(structure) -> StringIO:
    from eon import fileio as fio

    stream = StringIO()
    fio.savecon(stream, structure)
    stream.seek(0)
    return stream


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
            if "displacement.con" in job and "direction.dat" in job:
                from pyeonclient.bridge import structure_to_matter

                displacement_structure = _structure_from_geometry(
                    job["displacement.con"], "job['displacement.con']"
                )
                displacement = structure_to_matter(
                    displacement_structure, pot, params
                )
                product = job_obj.run_from_matter_displaced(
                    matter, displacement, mode
                )
            else:
                product = job_obj.run_from_matter(matter)
            saddle = job_obj.saddle
            if product is None or saddle is None:
                raise CommunicatorError(
                    "inprocess: process search returned no product or saddle"
                )
            return {
                "matter": product,
                "saddle": saddle,
                "min1": job_obj.min1,
                "min2": job_obj.min2,
                "energy": float(saddle.potential_energy),
                "potential_energy_reactant": float(
                    job_obj.min1.potential_energy
                ),
                "potential_energy_product": float(
                    job_obj.min2.potential_energy
                ),
                "potential_energy_saddle": float(
                    saddle.potential_energy
                ),
                "barrier_reactant_to_product": float(
                    job_obj.barrier_forward
                ),
                "barrier_product_to_reactant": float(
                    job_obj.barrier_reverse
                ),
                "displacement_saddle_distance": float(
                    pc.structure_distance(matter, saddle)
                ),
                "prefactor_reactant_to_product": float(
                    job_obj.prefactor_forward
                ),
                "prefactor_product_to_reactant": float(
                    job_obj.prefactor_reverse
                ),
                "force_calls_breakdown": {
                    "minimization": int(job_obj.force_calls_minimization),
                    "saddle": int(job_obj.force_calls_saddle),
                    "prefactors": int(job_obj.force_calls_prefactors),
                },
                "force_calls": int(product.force_calls),
                "status": int(job_obj.status),
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
    """Typed in-process result. ``results.dat`` is derived from this dict."""
    if cancelled:
        reason = "cancelled"
    elif status == 0:
        reason = "GOOD"
    else:
        reason = "FAIL"
    return {
        "status_code": status,
        "status_text": reason,
        "job_type": job_type,
        "potential_energy": energy,
        "force_calls": {"total": force_calls},
    }


def _results_dat(
    status: int,
    energy: float,
    force_calls: int,
    job_type: str,
    *,
    cancelled: bool = False,
) -> str:
    return _results_dat_from_job_result(
        _job_result(
            status, energy, force_calls, job_type, cancelled=cancelled
        )
    )


def _results_dat_from_job_result(job_result: dict) -> str:
    from eon_schema.jobs import job_result_legacy_dict

    data = job_result_legacy_dict(job_result)
    lines = []
    for key, val in data.items():
        if isinstance(val, float):
            lines.append(f"{val:.12e} {key}")
        else:
            lines.append(f"{val} {key}")
    return "\n".join(lines) + "\n"


class LocalInProcess(Communicator):
    """Run client work in-process via Matter (nanobind), not a subprocess."""

    results_are_synchronous = True

    def __init__(self, scratchpath, bundle_size=1, config=None):
        if config is None:
            raise TypeError("LocalInProcess requires a ConfigClass instance")
        Communicator.__init__(self, scratchpath, bundle_size, config=config)
        self._pc = _require_pyeonclient()
        self._finished: list[dict] = []
        self.token = CancelToken()
        self._in_submit = False
        self._cpp_token = None
        self._state_sessions: dict[str, dict[str, Any]] = {}
        self._state_minimizations: dict[str, int] = {}
        self._reuse_state_sessions = bool(
            getattr(config, "comm_reuse_state_sessions", True)
        )

    @property
    def state_session_count(self):
        return len(self._state_sessions)

    def state_minimization_count(self, state):
        return self._state_minimizations.get(str(state), 0)

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

        self._in_submit = True
        try:
            with _PROCESS_CWD_LOCK:
                previous = Path.cwd()
                try:
                    workdir = self._stage_potential_files()
                    if workdir is not None:
                        os.chdir(workdir)
                    self._submit_loop(data, invariants, pc)
                finally:
                    os.chdir(previous)
        finally:
            self._in_submit = False
            self.token.reset()

    def _stage_potential_files(self):
        from eon import fileio

        files = fileio.load_potfiles(
            getattr(self.config, "path_pot", "") or ""
        )
        if not files:
            return None
        workdir = Path(self.scratchpath) / ".inprocess-potential"
        workdir.mkdir(parents=True, exist_ok=True)
        root = workdir.resolve()
        for filename, (contents, permissions) in files.items():
            destination = (workdir / filename).resolve()
            if root not in destination.parents:
                raise CommunicatorError(
                    f"potential support file escapes workspace: {filename}"
                )
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(contents.getvalue())
            destination.chmod(permissions)
        return workdir

    def _submit_loop(self, data, invariants, pc):
        from pyeonclient.bridge import structure_to_matter, matter_to_structure

        for job in data:
            self.token.raise_if_cancelled()
            setup_force_calls = 0
            jid = job.get("id", "job")
            config_source = (
                job if "config.ini" in job or "config" in job else invariants
            )
            params = _params_from_invariants(pc, config_source)
            config_text = _config_text(config_source)
            config_fingerprint = _session_config_fingerprint(config_text)
            job_kind = getattr(params, "job", pc.JobType.Minimization)
            try:
                structure = _structure_from_job(job, invariants)
            except CommunicatorError:
                logger.exception("inprocess: missing geometry for %s", jid)
                raise
            except Exception as e:
                logger.exception("inprocess: failed to read geometry for %s", jid)
                raise CommunicatorError(str(e)) from e

            stateful = job_kind in (
                pc.JobType.Process_Search,
                pc.JobType.Saddle_Search,
            )
            state_id = _state_key(job)
            fingerprint = _geometry_fingerprint(structure)
            session = (
                self._state_sessions.get(state_id)
                if stateful and self._reuse_state_sessions
                else None
            )
            if (
                session is None
                or session["fingerprint"] != fingerprint
                or session["config_fingerprint"] != config_fingerprint
            ):
                pot = pc.make_potential(params)
                matter = structure_to_matter(structure, pot, params)
                if stateful and bool(
                    getattr(params, "process_search_minimize_first", False)
                ):
                    setup_calls_before = int(pot.force_call_counter)
                    matter, converged = matter.relax(
                        inplace=True,
                        quiet=True,
                        write_movie=False,
                        checkpoint=False,
                    )
                    if not converged:
                        raise CommunicatorError(
                            f"inprocess: reactant minimization failed for state {state_id}"
                        )
                    setup_force_calls = (
                        int(pot.force_call_counter) - setup_calls_before
                    )
                    self._state_minimizations[state_id] = (
                        self._state_minimizations.get(state_id, 0) + 1
                    )
                if stateful and self._reuse_state_sessions:
                    self._state_sessions[state_id] = {
                        "potential": pot,
                        "matter": matter,
                        "parameters": params,
                        "fingerprint": fingerprint,
                        "config_fingerprint": config_fingerprint,
                    }
            else:
                pot = session["potential"]
                matter = session["matter"]
            if hasattr(pc, "CancelToken") and hasattr(matter, "set_cancel_token"):
                self._cpp_token = pc.CancelToken()
                matter.set_cancel_token(self._cpp_token)
                if self.token.cancelled:
                    self._cpp_token.request()
            cancelled_exc = getattr(pc, "JobCancelled", None)
            calls_before = int(pot.force_call_counter)
            try:
                payload = _run_inprocess_job(
                    pc, job_kind, matter, pot, params, job, token=self.token
                )
            except Exception as exc:
                if cancelled_exc is None or not isinstance(exc, cancelled_exc):
                    raise
                payload = {
                    "matter": matter,
                    "energy": float(matter.potential_energy),
                    "force_calls": int(getattr(matter, "force_calls", 0) or 0),
                    "status": 1,
                    "job_type": str(getattr(job_kind, "name", job_kind)),
                    "converged": False,
                    "cancelled": True,
                }
            finally:
                self._cpp_token = None
            matter = payload["matter"]
            if "force_calls_breakdown" in payload:
                payload["force_calls_breakdown"]["minimization"] += (
                    setup_force_calls
                )
                payload["force_calls"] = sum(
                    payload["force_calls_breakdown"].values()
                )
            else:
                payload["force_calls"] = (
                    int(pot.force_call_counter)
                    - calls_before
                    + setup_force_calls
                )
            out = matter_to_structure(matter)

            energy = float(payload["energy"])
            fcalls = int(payload["force_calls"])
            status = int(payload["status"])
            jname = str(payload["job_type"])
            cancelled = bool(payload.get("cancelled"))
            job_result = _job_result(
                status, energy, fcalls, jname, cancelled=cancelled
            )
            if "force_calls_breakdown" in payload:
                job_result["force_calls"].update(
                    payload["force_calls_breakdown"]
                )
            for key in (
                "potential_energy_reactant",
                "potential_energy_product",
                "potential_energy_saddle",
                "barrier_reactant_to_product",
                "barrier_product_to_reactant",
                "displacement_saddle_distance",
                "prefactor_reactant_to_product",
                "prefactor_product_to_reactant",
            ):
                if key in payload:
                    job_result[key] = payload[key]
            results = StringIO(
                _results_dat_from_job_result(job_result)
            )

            rec = {
                "id": jid,
                "number": 0,
                "name": str(jid),
                "product": _conframe_of(out),
                "min.con": _con_io(out),
                "job_result": job_result,
                "results.dat": results,
                "_matter": matter,
                "_structure": out,
                "_energy": energy,
                "_converged": bool(payload.get("converged", status == 0)),
            }
            if "min1" in payload:
                reactant = matter_to_structure(payload["min1"])
                rec["reactant"] = _conframe_of(reactant)
                rec["reactant.con"] = _con_io(reactant)
            rec["product.con"] = _con_io(out)
            if "direction.dat" in job:
                rec["mode.dat"] = StringIO(
                    _con_text_of(job["direction.dat"])
                )
            saddle = payload.get("saddle")
            if saddle is not None:
                saddle_structure = matter_to_structure(saddle)
                rec["saddle"] = _conframe_of(saddle_structure)
                rec["saddle.con"] = _con_io(saddle_structure)
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
