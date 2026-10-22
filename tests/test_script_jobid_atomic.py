"""A killed job-id rewrite leaves the previous map in place."""

import pickle

import pytest

from eon.communicator import Script


def test_killed_jobid_rewrite_keeps_the_previous_map(tmp_path, monkeypatch):
    path = tmp_path / "script_job_ids"
    with path.open("wb") as handle:
        pickle.dump({9: "old"}, handle)
    comm = Script.__new__(Script)
    comm.job_id_path = str(path)
    comm.jobids = {1: "4"}

    def killed(_src, _dst):
        raise OSError("killed")

    monkeypatch.setattr("eon.fileio.os.replace", killed)
    with pytest.raises(OSError, match="killed"):
        comm.save_jobids()
    with path.open("rb") as handle:
        assert pickle.load(handle) == {9: "old"}
