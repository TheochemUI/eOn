"""The MPI communicator harvests returned job directories one by one."""
from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np

from eon.communicator import MPI


class _FakeComm:
    """Queues job paths the way client ranks send them on tag 0."""

    def __init__(self, paths):
        self.paths = list(paths)

    def Iprobe(self, source, tag, status):
        return bool(self.paths)

    def Recv(self, buf_and_type, source, tag):
        buf = buf_and_type[0]
        raw = self.paths.pop(0).encode()
        buf[: len(raw)] = np.frombuffer(raw, dtype="S1")


def _fake_mpi4py(monkeypatch):
    mpi = types.SimpleNamespace(
        Status=lambda: types.SimpleNamespace(source=1),
        ANY_SOURCE=-1,
        CHARACTER=object(),
    )
    pkg = types.ModuleType("mpi4py")
    pkg.MPI = mpi
    monkeypatch.setitem(sys.modules, "mpi4py", pkg)
    monkeypatch.setitem(sys.modules, "mpi4py.MPI", mpi)


def test_missing_returned_directory_is_skipped(tmp_path, monkeypatch):
    _fake_mpi4py(monkeypatch)
    scratch = tmp_path / "scratch"
    results = tmp_path / "results"
    scratch.mkdir()
    results.mkdir()
    job = scratch / "0_1"
    job.mkdir()
    (job / "results.dat").write_text("0 termination_reason\n")

    comm = MPI.__new__(MPI)
    comm.scratchpath = str(scratch)
    comm.bundle_size = 1
    comm.config = types.SimpleNamespace(debug_keep_all_results=False)
    comm.comm = _FakeComm([str(scratch / "0_0"), str(job)])

    got = list(comm.get_results(str(results), lambda name: True))

    assert [r["name"] for r in got] == ["0_1"]
    assert (Path(results) / "0_1" / "results.dat").is_file()
    assert not comm.comm.paths
