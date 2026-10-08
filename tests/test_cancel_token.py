"""A cancelled token stops the next in-process job."""

from __future__ import annotations

import sys
import types

import pytest

from eon.cancel import Cancelled
from eon.communicator_inprocess import LocalInProcess


class _WireStructure:
    """The wire lists read length, atom ids, and free axes from the structure."""

    atom_ids: list[int] = []
    _free: list[float] = []

    def __len__(self) -> int:
        return 0

    def to_conframe(self) -> str:
        return "frame"


class _PC:
    class JobType:
        Minimization = 1
        Unknown = 0
        Point = 2
        Process_Search = 3
        Saddle_Search = 4

    job = JobType.Point

    def make_potential(self, params):
        return object()


def test_cancel_state_stops_the_following_job(monkeypatch):
    monkeypatch.setattr(
        "eon.communicator_inprocess._require_pyeonclient", lambda: _PC()
    )
    monkeypatch.setattr(
        "eon.communicator_inprocess._params_from_invariants",
        lambda pc, invariants: types.SimpleNamespace(job=_PC.JobType.Point),
    )
    monkeypatch.setattr(
        "eon.communicator_inprocess._structure_from_job",
        lambda job, key="pos.con": object(),
    )
    bridge = types.ModuleType("pyeonclient.bridge")

    class _Matter:
        potential_energy = 0.0
        force_calls = 1
        positions = None

    bridge.structure_to_matter = lambda *args, **kwargs: _Matter()
    bridge.matter_to_structure = lambda matter: _WireStructure()
    monkeypatch.setitem(sys.modules, "pyeonclient.bridge", bridge)

    comm = LocalInProcess("scratch", config=object())
    ran = []

    def _run(*args, **kwargs):
        ran.append(kwargs.get("token"))
        assert comm.cancel_state("a") == 1
        return {
            "matter": _Matter(),
            "energy": -1.0,
            "force_calls": 1,
            "status": 0,
            "job_type": "point",
            "converged": True,
        }

    monkeypatch.setattr("eon.communicator_inprocess._run_inprocess_job", _run)
    jobs = [{"id": "a", "structure": object()}, {"id": "b", "structure": object()}]
    with pytest.raises(Cancelled):
        comm.submit_jobs(jobs, {})
    assert len(ran) == 1
    assert len(comm._finished) == 1
    assert comm.token.cancelled is False
    assert comm.cancel_state("a") == 0

    ran.clear()
    comm.submit_jobs([jobs[0]], {})
    assert len(ran) == 1


def test_idle_cancel_does_not_block_later_jobs(monkeypatch):
    monkeypatch.setattr(
        "eon.communicator_inprocess._require_pyeonclient", lambda: _PC()
    )
    monkeypatch.setattr(
        "eon.communicator_inprocess._params_from_invariants",
        lambda pc, invariants: types.SimpleNamespace(job=_PC.JobType.Point),
    )
    monkeypatch.setattr(
        "eon.communicator_inprocess._structure_from_job",
        lambda job, key="pos.con": object(),
    )
    bridge = types.ModuleType("pyeonclient.bridge")

    class _Matter:
        potential_energy = 0.0
        force_calls = 1
        positions = None

    bridge.structure_to_matter = lambda *args, **kwargs: _Matter()
    bridge.matter_to_structure = lambda matter: _WireStructure()
    monkeypatch.setitem(sys.modules, "pyeonclient.bridge", bridge)

    comm = LocalInProcess("scratch", config=object())
    assert comm.cancel_state(3) == 0
    ran = []

    def _run(*args, **kwargs):
        ran.append(kwargs.get("token"))
        return {
            "matter": _Matter(),
            "energy": -1.0,
            "force_calls": 1,
            "status": 0,
            "job_type": "point",
            "converged": True,
        }

    monkeypatch.setattr("eon.communicator_inprocess._run_inprocess_job", _run)
    jobs = [{"id": "a", "structure": object()}, {"id": "b", "structure": object()}]
    comm.submit_jobs(jobs, {})
    assert len(ran) == 2
    assert comm.token.cancelled is False
