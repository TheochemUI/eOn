"""A finished job replaces an existing results directory."""

from __future__ import annotations

from types import SimpleNamespace

from eon.communicator import Local


def _local(tmp_path, keep_debug=False):
    return Local(
        str(tmp_path / "scratch"),
        "/usr/bin/true",
        1,
        1,
        config=SimpleNamespace(
            debug_keep_all_results=keep_debug,
            debug_results_path="debug",
            path_root=str(tmp_path),
            path_pot="",
        ),
    )


def _plant(tmp_path):
    scratch = tmp_path / "scratch"
    job = scratch / "0_0"
    job.mkdir(parents=True)
    (job / "results.dat").write_text("new\n")
    old = tmp_path / "results" / "0_0"
    old.mkdir(parents=True)
    (old / "results.dat").write_text("old\n")
    return job


def test_finished_job_replaces_existing_results(tmp_path):
    _plant(tmp_path)
    comm = _local(tmp_path)
    list(comm.get_results(str(tmp_path / "results"), lambda name: True))
    dest = tmp_path / "results" / "0_0"
    assert (dest / "results.dat").read_text() == "new\n"
    assert not (dest / "0_0").exists()


def test_debug_copy_replaces_an_existing_directory(tmp_path):
    _plant(tmp_path)
    debug = tmp_path / "debug" / "0_0"
    debug.mkdir(parents=True)
    (debug / "stale.txt").write_text("stale\n")
    comm = _local(tmp_path, keep_debug=True)
    list(comm.get_results(str(tmp_path / "results"), lambda name: True))
    assert (debug / "results.dat").read_text() == "new\n"
    assert not (debug / "stale.txt").exists()
    assert (tmp_path / "results" / "0_0" / "results.dat").read_text() == "new\n"
