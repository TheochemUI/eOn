"""amsel discover_decide for a superbasin graph and for one state's table.

``discover_decide_for_superbasin`` is the only gate. It calls
``amsel.discover_decide_status`` and returns one status from
``DISCOVER_DECIDE_STATUSES``.

* ``accepted`` / ``retightened`` — the candidate is one basin.
* ``split_required`` — use the primary transient set anchored at the entry.
* ``rejected_no_metastable_basin`` — ordinary single-state KMC.
* ``fallback_single`` — a failed amsel call (``available`` is true).
* ``unavailable`` — amsel is not importable, or ``on_error`` is
  ``unavailable_mcamc`` (``available`` is false).

``[amsel] discover_decide`` does not require ``use_mcamc``.
``amsel_discover_exit`` reads the current state's process table. A barrier
strictly below ``e_min_init`` is in-basin. A higher barrier is an exit.
A table with no in-basin edge still exits: the transient set is the entry
state and the products are absorbing. A product column of -1 is that exit
too: the absorbing label is a 32-bit id, and the hop creates the product
state from the process id. The exit time and channel come from
MRM (mean time) or FPTA (one sample). The kernels are the ``amsel``
package. A missing package logs ``unavailable`` and returns no exit.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
from pathlib import Path
from typing import Any, Callable, Mapping

logger = logging.getLogger("superbasin.amsel_gate")

# One status vocabulary. ``unavailable`` pairs with available false.
# ``fallback_single`` pairs with available true. The pair
# (unavailable, true) is not a result of this module.
DISCOVER_DECIDE_STATUSES = (
    "unavailable",
    "fallback_single",
    "accepted",
    "retightened",
    "split_required",
    "rejected_no_metastable_basin",
    "rejected",
    "split_cached",
)
_EXIT_STATUSES = ("accepted", "retightened", "split_required")
_EXIT_KERNELS = ("mrm", "fpta")
# amsel StateId is u32. State directories are numbered from 0. An unlinked
# product stays out of that range until the hop creates the directory.
_UNLINKED_STATE_BASE = 1 << 31


class AmselSuperbasinReject(RuntimeError):
    """Raised when discover_decide rejects the candidate superbasin."""

    def __init__(self, status: str, detail: str = ""):
        self.status = status
        super().__init__(detail or status)


def _process_barrier_eV(proc: Mapping[str, Any]) -> float:
    for key in ("barrier", "saddle_energy", "barrier_eV"):
        if key in proc and proc[key] is not None:
            try:
                return float(proc[key])
            except (TypeError, ValueError):
                pass
    # Fallback: infer from rate is unreliable; use large barrier so edge is
    # "high TS" unless rate is huge — prefer explicit barrier fields.
    return 1.0


def build_graph_from_superbasin(superbasin: Any, entry_number: int) -> tuple[
    list[int],
    list[tuple[int, int, float]],
    list[float],
]:
    """Return (candidate_states, rate_triples, ts_energies) for amsel."""
    candidates = [int(n) for n in superbasin.state_numbers]
    if int(entry_number) not in candidates:
        candidates = [int(entry_number)] + candidates
    rates: list[tuple[int, int, float]] = []
    barriers: list[float] = []
    for number in superbasin.state_numbers:
        procs = superbasin.state_dict[number].get_process_table()
        for process_key, proc in list(procs.items()):
            rate = float(proc.get("rate", 0.0) or 0.0)
            if rate <= 0.0:
                continue
            product = _product_id(proc, process_key)
            if product is None:
                continue
            rates.append((int(number), int(product), rate))
            barriers.append(_process_barrier_eV(proc))
    return candidates, rates, barriers


def discover_decide_for_superbasin(
    superbasin: Any,
    entry_state: Any,
    *,
    e_min_init: float = 0.5,
    e_min_step: float = 0.05,
    e_min_floor: float = 0.05,
    cv_threshold: float = 10.0,
    on_error: str = "fallback_single",
) -> dict[str, Any]:
    """Run amsel.discover_decide_status; return a structured result dict.

    Keys: status (str), primary_transient (list[int]|None), raw (tuple|None),
    available (bool). If amsel is not importable, available=False and status
    is ``unavailable`` (caller should proceed with legacy MCAMC).

    ``on_error`` decides what a failing amsel call returns: ``raise``
    re-raises, ``unavailable_mcamc`` reports amsel as unavailable so MCAMC
    runs as before, and ``fallback_single`` (the default) reports
    ``fallback_single`` so the caller takes an ordinary single-state step.
    """
    try:
        from amsel import discover_decide_status
    except ImportError:
        return {
            "available": False,
            "status": "unavailable",
            "primary_transient": None,
            "raw": None,
            "reason": "amsel not importable",
        }

    entry = int(entry_state.number)
    candidates, rate_triples, barriers = build_graph_from_superbasin(
        superbasin, entry
    )
    if not rate_triples:
        return {
            "available": True,
            "status": "rejected_no_metastable_basin",
            "primary_transient": None,
            "raw": None,
            "reason": "no positive-rate processes in superbasin tables",
        }

    # Align barrier vector with rate triples (required arity)
    while len(barriers) < len(rate_triples):
        barriers.append(0.5)
    barriers = barriers[: len(rate_triples)]
    # The cutoff separates fast in-basin edges from exits. It is not lifted
    # to the largest barrier: that would make every edge transient and leave
    # the basin without an absorbing border.
    e_init = float(e_min_init)

    try:
        raw = discover_decide_status(
            entry,
            [int(c) for c in candidates],
            [(int(a), int(b), float(r)) for a, b, r in rate_triples],
            [float(x) for x in barriers],
            e_init,
            float(e_min_step),
            float(e_min_floor),
            float(cv_threshold),
        )
    except Exception as exc:
        logger.warning("amsel discover_decide_status failed: %s", exc)
        policy = str(on_error or "fallback_single")
        if policy == "raise":
            raise
        if policy == "unavailable_mcamc":
            return {
                "available": False,
                "status": "unavailable",
                "primary_transient": None,
                "raw": None,
                "reason": f"{type(exc).__name__}: {exc}",
            }
        # fallback_single: do not look like a successful AMSEl skip
        return {
            "available": True,
            "status": "fallback_single",
            "primary_transient": None,
            "raw": None,
            "reason": f"{type(exc).__name__}: {exc}",
        }

    status = str(raw[0]) if isinstance(raw, (list, tuple)) and raw else str(raw)
    primary = None
    if isinstance(raw, (list, tuple)) and len(raw) > 1:
        # DiscoverDecisionStatus: (status, transient, absorbing, rates, ...)
        try:
            primary = [int(x) for x in raw[1]]
        except Exception:
            primary = None
    return {
        "available": True,
        "status": status,
        "primary_transient": primary,
        "raw": raw,
        "reason": "",
    }


def _split_cache_path(superbasin: Any) -> str | None:
    path = getattr(superbasin, "path", None)
    if not path:
        return None
    return str(path) + ".amsel_split.json"


def load_persisted_split(superbasin: Any) -> list[int] | None:
    cache = _split_cache_path(superbasin)
    if not cache or not Path(cache).is_file():
        return None
    try:
        with Path(cache).open(encoding="utf-8") as fh:
            data = json.load(fh)
        keep = [int(x) for x in data.get("keep", [])]
        return keep or None
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None


def persist_split(superbasin: Any, keep: list[int]) -> None:
    cache = _split_cache_path(superbasin)
    if not cache:
        return
    try:
        with Path(cache).open("w", encoding="utf-8") as fh:
            json.dump({"keep": [int(x) for x in keep]}, fh)
    except OSError as exc:
        logger.warning("amsel persist_split failed: %s", exc)


def apply_gate_to_superbasin(
    superbasin: Any,
    entry_state: Any,
    decision: Mapping[str, Any],
) -> str:
    """Mutate superbasin in-place for split_required; raise on reject.

    Returns the status string for logging.
    """
    status = str(decision.get("status", "unavailable"))
    if not decision.get("available") or status in (
        "unavailable",
        "accepted",
        "retightened",
        "fallback_single",
    ):
        return status
    cached = load_persisted_split(superbasin)
    if cached:
        keep = [n for n in superbasin.state_numbers if int(n) in set(cached)]
        if keep and keep != list(superbasin.state_numbers):
            superbasin.state_numbers = keep
            superbasin.states = [superbasin.state_dict[n] for n in keep]
        return "split_cached"
    if status == "split_required":
        primary = decision.get("primary_transient") or []
        if not primary:
            logger.warning(
                "amsel split_required but no primary_transient; keeping full basin"
            )
            return status
        primary_set = set(int(x) for x in primary)
        # Entry must remain for MCAMC indexing (st2i[entry_state.number]).
        entry_n = int(entry_state.number)
        primary_set.add(entry_n)
        # Restrict member list to primary transient states that exist
        keep = [n for n in superbasin.state_numbers if int(n) in primary_set]
        if not keep:
            logger.warning(
                "amsel primary_transient %s disjoint from superbasin; keeping full",
                primary,
            )
            return status
        logger.info(
            "amsel discover_decide=split_required: restricting superbasin %s -> %s",
            superbasin.state_numbers,
            keep,
        )
        superbasin.state_numbers = keep
        superbasin.states = [superbasin.state_dict[n] for n in keep]
        # Drop dict entries not in keep (optional consistency)
        for n in list(superbasin.state_dict.keys()):
            if n not in keep:
                del superbasin.state_dict[n]
        # Persist so in-memory membership matches on-disk after resume.
        if hasattr(superbasin, "write_data"):
            try:
                superbasin.write_data()
            except Exception as exc:
                logger.warning("amsel split: write_data failed: %s", exc)
        persist_split(superbasin, [int(n) for n in keep])
        return status
    if status in ("rejected_no_metastable_basin", "rejected"):
        raise AmselSuperbasinReject(
            status,
            "amsel discover_decide rejected superbasin %s (entry %s): %s"
            % (superbasin.state_numbers, entry_state.number, decision.get("reason", "")),
        )
    return status


def status_pair(status: str, available: bool) -> bool:
    """Return true when ``(status, available)`` is in the one vocabulary."""
    text = str(status)
    if text not in DISCOVER_DECIDE_STATUSES:
        return False
    if text == "unavailable":
        return available is False
    if text == "fallback_single":
        return available is True
    return True


class _ProcessView:
    """State id with a fixed process table. Border states carry an empty table."""

    def __init__(self, number: int, table: Mapping[Any, Any] | None = None):
        self.number = int(number)
        self._table = {} if table is None else table

    def get_process_table(self):
        return self._table


class AmselExit:
    """One basin exit chosen by MRM or FPTA."""

    __slots__ = (
        "exit_number",
        "proc_id",
        "product_number",
        "mean_time",
        "time_is_sample",
        "barrier",
        "rate",
        "kernel",
    )

    def __init__(
        self,
        exit_number: int,
        proc_id: int,
        product_number: int,
        mean_time: float,
        time_is_sample: bool,
        barrier: float,
        rate: float,
        kernel: str,
    ):
        self.exit_number = int(exit_number)
        self.proc_id = proc_id
        self.product_number = int(product_number)
        self.mean_time = float(mean_time)
        self.time_is_sample = bool(time_is_sample)
        self.barrier = float(barrier)
        self.rate = float(rate)
        self.kernel = str(kernel)


def _resolve_state(resolve_state: Callable[[int], Any] | None, number: int) -> Any:
    if resolve_state is None:
        return None
    try:
        return resolve_state(int(number))
    except Exception:
        return None


def _unlinked_state_id(process_key: Any, proc: Mapping[str, Any]) -> int:
    """32-bit absorbing label for a product column that is still -1.

    The process id does not fit in ``amsel``'s ``u32`` state id. The label
    is only the graph id. The hop creates the product state from the
    process key.
    """
    text = "\0".join(
        (
            str(process_key),
            str(proc.get("barrier")),
            str(proc.get("product_energy")),
            str(proc.get("saddle_energy")),
        )
    )
    digest = hashlib.blake2s(text.encode("utf-8"), digest_size=4).digest()
    offset = int.from_bytes(digest, "little") & 0x7FFFFFFF
    return _UNLINKED_STATE_BASE + offset


def _product_id(proc: Mapping[str, Any], process_key: Any = None) -> int | None:
    product = proc.get("product")
    if product is None:
        return None
    try:
        pid = int(product)
    except (TypeError, ValueError):
        return None
    if pid < 0:
        return _unlinked_state_id(process_key, proc)
    return pid


def basin_view_from_state(
    entry_state: Any,
    resolve_state: Callable[[int], Any] | None,
    e_min_init: float,
) -> Any:
    """States reached by barriers below ``e_min_init``, plus exit products.

    Exit products are candidate ids with an empty process table, so their
    own edges are not fed to discover_decide. A barrier strictly below
    ``e_min_init`` enqueues the product when the state list can resolve it.
    """
    from types import SimpleNamespace

    tables: dict[int, Any] = {int(entry_state.number): entry_state}
    pending = [entry_state]
    border: list[int] = []
    cutoff = float(e_min_init)
    while pending:
        state = pending.pop()
        for process_key, proc in list(state.get_process_table().items()):
            pid = _product_id(proc, process_key)
            if pid is None:
                continue
            if _process_barrier_eV(proc) < cutoff:
                if pid in tables:
                    continue
                other = None
                if pid < _UNLINKED_STATE_BASE:
                    other = _resolve_state(resolve_state, pid)
                if other is None:
                    other = _ProcessView(pid, {})
                tables[pid] = other
                pending.append(other)
            elif pid not in tables and pid not in border:
                border.append(pid)
    numbers = list(tables)
    for pid in border:
        if pid not in tables:
            tables[pid] = _ProcessView(pid, {})
            numbers.append(pid)
    return SimpleNamespace(
        state_numbers=numbers,
        state_dict=tables,
        states=[tables[n] for n in numbers],
        id=None,
        path=None,
    )


def _has_in_basin_edge(view: Any, e_min_init: float) -> bool:
    cutoff = float(e_min_init)
    for state in view.state_dict.values():
        for process_key, proc in state.get_process_table().items():
            if _product_id(proc, process_key) is None:
                continue
            if _process_barrier_eV(proc) < cutoff:
                return True
    return False


def _log_status(decision: Mapping[str, Any]) -> str:
    status = str(decision.get("status", "unavailable"))
    logger.info(
        "amsel discover_decide status=%s available=%s primary=%s",
        status,
        decision.get("available"),
        decision.get("primary_transient"),
    )
    return status


def _partition(decision: Mapping[str, Any]) -> tuple[list[int], list[int], list] | None:
    raw = decision.get("raw")
    if not isinstance(raw, (list, tuple)) or len(raw) < 4:
        return None
    try:
        transient = [int(x) for x in raw[1]]
        absorbing = [int(x) for x in raw[2]]
        rates = [(int(a), int(b), float(r)) for a, b, r in raw[3]]
    except (TypeError, ValueError):
        return None
    if not transient or not absorbing or not rates:
        return None
    return transient, absorbing, rates


def _sample_index(weights: list[float], draw: float) -> int | None:
    total = 0.0
    for weight in weights:
        if weight > 0.0 and math.isfinite(weight):
            total += weight
    if total <= 0.0:
        return None
    mark = min(max(float(draw), 0.0), 1.0) * total
    walked = 0.0
    last = None
    for index, weight in enumerate(weights):
        if weight > 0.0 and math.isfinite(weight):
            walked += weight
            last = index
            if walked >= mark:
                return index
    return last


def _invoke_kernel(
    name: str,
    transient: list[int],
    absorbing: list[int],
    rates: list,
    entry: int,
    draw: float,
) -> tuple[float, list[float], bool] | None:
    import amsel

    if name == "mrm" and hasattr(amsel, "mrm"):
        tau, rate_to_absorbing, _x = amsel.mrm(transient, absorbing, rates, int(entry))
        return float(tau), [float(x) for x in rate_to_absorbing], False
    if name == "fpta" and hasattr(amsel, "fpta"):
        t_exit, weights = amsel.fpta(
            transient, absorbing, rates, int(entry), float(draw)
        )
        return float(t_exit), [float(x) for x in weights], True
    if not hasattr(amsel, "AmcProblem"):
        return None
    problem = amsel.AmcProblem(
        transient=transient, absorbing=absorbing, rates=rates
    )
    if name == "mrm":
        result = problem.mrm(int(entry))
        return float(result.tau_total), [float(x) for x in result.rate_to_absorbing], False
    result = problem.fpta(int(entry), float(draw))
    return float(result.t_exit), [float(x) for x in result.weights], True


def _kernel_exit(
    prefer: str,
    transient: list[int],
    absorbing: list[int],
    rates: list,
    entry: int,
    uniform: Callable[[], float],
) -> tuple[str, int, float, bool] | None:
    """Pick one absorbing state. MRM is the mean clock. FPTA is one sample."""
    order = (prefer, "fpta" if prefer == "mrm" else "mrm")
    errors: list[str] = []
    for name in order:
        if name not in _EXIT_KERNELS:
            continue
        try:
            draw = float(uniform())
            # FPTA rejects the closed ends of the unit interval.
            draw = min(max(draw, 1e-16), 1.0 - 1e-16)
            got = _invoke_kernel(name, transient, absorbing, rates, entry, draw)
        except Exception as exc:
            errors.append("%s: %s: %s" % (name, type(exc).__name__, exc))
            continue
        if got is None:
            errors.append("%s: missing" % name)
            continue
        clock, weights, time_is_sample = got
        if len(weights) != len(absorbing) or not math.isfinite(clock) or clock < 0.0:
            errors.append("%s: bad result" % name)
            continue
        index = _sample_index(weights, float(uniform()))
        if index is None:
            errors.append("%s: empty weights" % name)
            continue
        return name, int(absorbing[index]), float(clock), bool(time_is_sample)
    if errors:
        logger.warning("amsel exit kernels failed: %s", "; ".join(errors))
    return None


def _direct_escape_partition(
    view: Any, entry: int, e_min_init: float
) -> tuple[list[int], list[int], list] | None:
    """Entry state plus its exit edges, when nothing is below the cutoff.

    The transient set is that one state. Products of barriers at or above
    ``e_min_init`` are absorbing. A non-positive rate is dropped.
    """
    cutoff = float(e_min_init)
    state = view.state_dict.get(int(entry))
    if state is None:
        return None
    absorbing: list[int] = []
    seen: set[int] = set()
    rates: list[tuple[int, int, float]] = []
    for process_key, proc in state.get_process_table().items():
        pid = _product_id(proc, process_key)
        if pid is None or pid == int(entry):
            continue
        if _process_barrier_eV(proc) < cutoff:
            continue
        rate = float(proc.get("rate", 0.0) or 0.0)
        if rate <= 0.0:
            continue
        if pid not in seen:
            seen.add(pid)
            absorbing.append(pid)
        rates.append((int(entry), int(pid), rate))
    if not absorbing or not rates:
        return None
    return [int(entry)], absorbing, rates


def _exit_process(
    state_dict: Mapping[int, Any],
    transient: list[int],
    absorbing_id: int,
    e_min_init: float,
) -> tuple[int, Any, float, float] | None:
    """Fastest exit process from the transient set into ``absorbing_id``."""
    cutoff = float(e_min_init)
    best: tuple[float, int, Any, float] | None = None
    for number in transient:
        state = state_dict.get(int(number))
        if state is None:
            continue
        for process_key, proc in state.get_process_table().items():
            if _product_id(proc, process_key) != int(absorbing_id):
                continue
            barrier = _process_barrier_eV(proc)
            if barrier < cutoff:
                continue
            rate = float(proc.get("rate", 0.0) or 0.0)
            if best is None or rate > best[0]:
                best = (rate, int(number), process_key, barrier)
    if best is None:
        return None
    rate, number, pid, barrier = best
    return number, pid, barrier, rate


def _amsel_floats(config: Any) -> dict[str, Any]:
    return {
        "e_min_init": float(getattr(config, "amsel_e_min_init", 0.5)),
        "e_min_step": float(getattr(config, "amsel_e_min_step", 0.05)),
        "e_min_floor": float(getattr(config, "amsel_e_min_floor", 0.05)),
        "cv_threshold": float(getattr(config, "amsel_cv_threshold", 10.0)),
        "on_error": str(getattr(config, "amsel_on_error", "fallback_single")),
    }


def amsel_discover_exit(
    entry_state: Any,
    resolve_state: Callable[[int], Any] | None,
    config: Any,
    uniform: Callable[[], float] | None = None,
) -> AmselExit | None:
    """Log ``amsel discover_decide status`` and maybe return an MRM or FPTA exit.

    The log line is written for every call, including a table with one
    state and a missing ``amsel`` package. An in-basin edge is a barrier
    strictly below ``e_min_init``. With such an edge, an exit is returned
    when discover accepts the basin, retightens it, or splits it. With no
    such edge, a barrier at or above the cutoff is a direct exit. The
    transient set is the entry state, the products are absorbing, and the
    logged status is ``accepted``. A product column of -1 uses a 32-bit
    absorbing label. The process id is unchanged.

    ``debug_use_mean_time`` selects MRM. Otherwise the clock is one FPTA
    sample. The other kernel is used when the selected one is missing.
    A missing package logs ``unavailable`` and returns no exit.
    """
    if uniform is None:
        import random

        uniform = random.random
    knobs = _amsel_floats(config)
    view = basin_view_from_state(entry_state, resolve_state, knobs["e_min_init"])
    decision = discover_decide_for_superbasin(
        view,
        entry_state,
        e_min_init=knobs["e_min_init"],
        e_min_step=knobs["e_min_step"],
        e_min_floor=knobs["e_min_floor"],
        cv_threshold=knobs["cv_threshold"],
        on_error=knobs["on_error"],
    )
    entry = int(entry_state.number)
    # A lone exit has no fast edge for discover to call a basin. The
    # absorbing states are the products, and the transient set is the entry.
    # A missing package stays ``unavailable`` and does not hop.
    direct = None
    if str(decision.get("status", "unavailable")) != "unavailable" and not _has_in_basin_edge(
        view, knobs["e_min_init"]
    ):
        direct = _direct_escape_partition(view, entry, knobs["e_min_init"])
    if direct is not None:
        transient, absorbing, rates = direct
        decision = {
            "available": True,
            "status": "accepted",
            "primary_transient": list(transient),
            "raw": ("accepted", list(transient), list(absorbing), list(rates)),
        }
    status = _log_status(decision)
    if not status_pair(status, bool(decision.get("available"))):
        logger.warning(
            "amsel discover_decide status=%s available=%s is outside the vocabulary",
            status,
            decision.get("available"),
        )
    if status not in _EXIT_STATUSES or not decision.get("available"):
        return None
    if direct is None:
        parts = _partition(decision)
        if parts is None:
            return None
        transient, absorbing, rates = parts
        if entry not in transient:
            transient = [entry] + [n for n in transient if n != entry]
    prefer = "mrm" if bool(getattr(config, "debug_use_mean_time", False)) else "fpta"
    chosen = _kernel_exit(prefer, transient, absorbing, rates, entry, uniform)
    if chosen is None:
        return None
    kernel, absorbing_id, clock, time_is_sample = chosen
    proc = _exit_process(view.state_dict, transient, absorbing_id, knobs["e_min_init"])
    if proc is None:
        logger.warning(
            "amsel %s absorbing state %s has no exit process at or above e_min_init %.3f eV",
            kernel,
            absorbing_id,
            knobs["e_min_init"],
        )
        return None
    exit_number, proc_id, barrier, rate = proc
    logger.info(
        "amsel discover_decide exit kernel=%s product=%s time=%.6e",
        kernel,
        absorbing_id,
        clock,
    )
    return AmselExit(
        exit_number=exit_number,
        proc_id=proc_id,
        product_number=absorbing_id,
        mean_time=clock,
        time_is_sample=time_is_sample,
        barrier=barrier,
        rate=rate,
        kernel=kernel,
    )
