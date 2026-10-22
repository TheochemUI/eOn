"""script_job_ids survives a failed write and a truncated file."""

from __future__ import annotations

import pickle
from types import SimpleNamespace

import pytest

from eon.communicator import Script


def _comm(tmp_path):
    scripts = tmp_path / "scripts"
    scripts.mkdir(exist_ok=True)
    for name in ("queued_jobs.sh", "cancel_job.sh", "submit_job.sh"):
        path = scripts / name
        if not path.exists():
            path.write_text("#!/bin/sh\n")
    return Script(
        str(tmp_path / "scratch"),
        1,
        "eon",
        str(scripts),
        "queued_jobs.sh",
        "cancel_job.sh",
        "submit_job.sh",
        config=SimpleNamespace(
            debug_keep_all_results=False,
            debug_results_path="debug",
            path_root=str(tmp_path),
            path_pot="",
        ),
    )


def test_truncated_job_id_file_starts_empty(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    (scratch / "script_job_ids").write_bytes(b"\x80")
    comm = _comm(tmp_path)
    assert comm.jobids == {}


def test_failed_job_id_write_keeps_previous(tmp_path, monkeypatch):
    comm = _comm(tmp_path)
    comm.jobids = {1: "0"}
    comm.save_jobids()

    def boom(*_args, **_kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(pickle, "dump", boom)
    comm.jobids = {2: "1"}
    with pytest.raises(RuntimeError, match="disk full"):
        comm.save_jobids()
    assert _comm(tmp_path).jobids == {1: "0"}
