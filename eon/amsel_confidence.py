"""State confidence from amsel for ``[AKMC] confidence_scheme = amsel``.

``[amsel] state_confidence`` selects the score, and ``amsel.state_confidence``
computes it. A process whose product state also lists this state is known from
the far side. The score compares those processes with the ones this state's own
searches found.
"""

from __future__ import annotations

from typing import Any

STATE_CONFIDENCE_SCHEMES = ("heuristic", "good-turing", "reverse-capped", "two-sided")


def _catalogue(state: Any, ids: Any = None) -> dict[str, list]:
    state.load_process_table()
    procs = state.procs
    counts = state.get_proc_random_count()
    chosen = list(procs) if ids is None else [i for i in ids if i in procs]
    own = [i for i in chosen if counts.get(i, 0) > 0]
    implied = []
    for i in chosen:
        product = procs[i]["product"]
        if product < 0 or product == state.number:
            continue
        other = state.statelist.get_state(product)
        other.load_process_table()
        if any(p["product"] == state.number for p in other.procs.values()):
            implied.append((float(procs[i]["rate"]), float(procs[i]["barrier"])))
    return {
        "hits": [int(counts[i]) for i in own],
        "rates": [float(procs[i]["rate"]) for i in own],
        "barriers": [float(procs[i]["barrier"]) for i in own],
        "implied": implied,
    }


def amsel_state_confidence(
    state: Any, *, scheme: str, barrier_tolerance: float, ids: Any = None
) -> float:
    """Confidence of one state. ``ids`` limits it to these processes."""
    try:
        from amsel import state_confidence
    except ImportError as exc:
        raise RuntimeError(
            "[AKMC] confidence_scheme = amsel needs the amsel package"
        ) from exc
    network = (0, 0)
    if scheme == "two-sided":
        recovered = known = 0
        for number in range(state.statelist.get_num_states()):
            found, total = state_confidence(
                scheme,
                **_catalogue(state.statelist.get_state(number)),
                barrier_tolerance=barrier_tolerance,
            )[1]
            recovered += found
            known += total
        network = (recovered, known)
    return float(
        state_confidence(
            scheme,
            **_catalogue(state, ids),
            network=network,
            barrier_tolerance=barrier_tolerance,
        )[0]
    )
