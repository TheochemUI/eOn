"""Optional amsel discover_decide gate before Superbasin.step MCAMC solve (L6).

When enabled and the ``amsel`` Python package is importable, builds the
candidate basin graph from the superbasin process tables and calls
``amsel.discover_decide_status``. Outcomes:

* ``accepted`` / ``retightened`` — proceed with the existing MCAMC ``step``.
* ``split_required`` — restrict MCAMC to the primary transient set returned
  by amsel (primary basin only); log sibling components.
* ``rejected_no_metastable_basin`` — raise ``AmselSuperbasinReject`` so the
  AKMC driver can fall back to ordinary single-state KMC (no superbasin
  acceleration on a non-metastable candidate).

This is Python (eOn's Superbasin controller is Python). There is no C++
``Superbasin::step`` in upstream eOn; eonclient Hessian/saddle jobs remain
separate. Optional soft-dependency: if amsel is missing or the gate is off,
behaviour is unchanged from legacy MCAMC.
"""
from __future__ import annotations

import logging
from typing import Any, Mapping

logger = logging.getLogger("superbasin.amsel_gate")


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
        for _pid, proc in list(procs.items()):
            rate = float(proc.get("rate", 0.0) or 0.0)
            if rate <= 0.0:
                continue
            product = proc.get("product")
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
) -> dict[str, Any]:
    """Run amsel.discover_decide_status; return a structured result dict.

    Keys: status (str), primary_transient (list[int]|None), raw (tuple|None),
    available (bool). If amsel is not importable, available=False and status
    is ``unavailable`` (caller should proceed with legacy MCAMC).
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
    e_init = float(e_min_init)
    if barriers:
        floor = max(barriers) + 0.05
        if floor > e_init:
            logger.info(
                "amsel e_min_init raised from %s to %s (max barrier + 0.05 eV)",
                e_init,
                floor,
            )
            e_init = floor

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
        return {
            "available": True,
            "status": "unavailable",
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


def apply_gate_to_superbasin(
    superbasin: Any,
    entry_state: Any,
    decision: Mapping[str, Any],
) -> str:
    """Mutate superbasin in-place for split_required; raise on reject.

    Returns the status string for logging.
    """
    status = str(decision.get("status", "unavailable"))
    if not decision.get("available") or status in ("unavailable", "accepted", "retightened"):
        return status
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
        return status
    if status in ("rejected_no_metastable_basin", "rejected"):
        raise AmselSuperbasinReject(
            status,
            "amsel discover_decide rejected superbasin %s (entry %s): %s"
            % (superbasin.state_numbers, entry_state.number, decision.get("reason", "")),
        )
    return status
