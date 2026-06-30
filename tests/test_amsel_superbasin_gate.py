"""Unit tests for eon.amsel_superbasin_gate (optional amsel dependency)."""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from eon.amsel_superbasin_gate import (
    AmselSuperbasinReject,
    apply_gate_to_superbasin,
    build_graph_from_superbasin,
    discover_decide_for_superbasin,
)


class _FakeState:
    def __init__(self, number, procs):
        self.number = number
        self._procs = procs

    def get_process_table(self):
        return self._procs


def _fake_superbasin():
    # 0 <-> 1 fast, 0 -> 2 exit (product outside basin for MCAMC would be 2)
    s0 = _FakeState(
        0,
        {
            0: {"rate": 1e10, "product": 1, "barrier": 0.1},
            1: {"rate": 1e8, "product": 2, "barrier": 0.3},
        },
    )
    s1 = _FakeState(
        1,
        {
            0: {"rate": 1e10, "product": 0, "barrier": 0.1},
        },
    )
    sb = SimpleNamespace(
        state_numbers=[0, 1],
        states=[s0, s1],
        state_dict={0: s0, 1: s1},
        id=1,
    )
    return sb


def test_build_graph_from_superbasin_edges():
    sb = _fake_superbasin()
    cands, rates, barriers = build_graph_from_superbasin(sb, 0)
    assert 0 in cands and 1 in cands
    assert any(r[0] == 0 and r[1] == 1 for r in rates)
    assert len(barriers) == len(rates)


def test_discover_decide_for_superbasin_runs_or_unavailable():
    sb = _fake_superbasin()
    entry = SimpleNamespace(number=0)
    decision = discover_decide_for_superbasin(sb, entry)
    assert "status" in decision
    assert "available" in decision
    # With amsel installed in dev env, expect a real status; otherwise unavailable
    if decision["available"] and decision["status"] != "unavailable":
        assert decision["status"] in (
            "accepted",
            "retightened",
            "split_required",
            "rejected_no_metastable_basin",
        )


def test_apply_gate_split_restricts_members():
    sb = _fake_superbasin()
    entry = SimpleNamespace(number=0)
    decision = {
        "available": True,
        "status": "split_required",
        "primary_transient": [0],
        "raw": None,
        "reason": "",
    }
    status = apply_gate_to_superbasin(sb, entry, decision)
    assert status == "split_required"
    assert sb.state_numbers == [0]


def test_apply_gate_reject_raises():
    sb = _fake_superbasin()
    entry = SimpleNamespace(number=0)
    with pytest.raises(AmselSuperbasinReject):
        apply_gate_to_superbasin(
            sb,
            entry,
            {
                "available": True,
                "status": "rejected_no_metastable_basin",
                "primary_transient": None,
                "reason": "test",
            },
        )


def test_split_keeps_entry_state():
    from eon.amsel_superbasin_gate import apply_gate_to_superbasin

    class _St:
        def __init__(self, n):
            self.number = n

        def get_process_table(self):
            return {}

    class _SB:
        def __init__(self):
            self.state_numbers = [0, 1, 2]
            self.state_dict = {0: _St(0), 1: _St(1), 2: _St(2)}
            self.states = [self.state_dict[n] for n in self.state_numbers]
            self.written = False

        def write_data(self):
            self.written = True

    sb = _SB()
    entry = _St(0)
    # primary omits entry 0 — gate must re-add it
    status = apply_gate_to_superbasin(
        sb,
        entry,
        {
            "available": True,
            "status": "split_required",
            "primary_transient": [1, 2],
            "reason": "",
        },
    )
    assert status == "split_required"
    assert 0 in sb.state_numbers
    assert sb.written
