"""Communicator bundles files and the base class raises the builtin error."""

from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from eon.communicator import (
    Communicator,
    CommunicatorError,
    Local,
    get_communicator,
    reset_communicators,
)


def test_base_methods_raise_builtin_not_implemented(tmp_path):
    comm = Communicator(tmp_path, 1, config=SimpleNamespace())
    with pytest.raises(NotImplementedError) as caught:
        comm.submit_jobs([], {})
    assert type(caught.value) is NotImplementedError
    with pytest.raises(NotImplementedError):
        comm.get_results(tmp_path)
    with pytest.raises(NotImplementedError):
        comm.cancel_state(0)


def test_single_results_file_is_not_a_bundle(tmp_path):
    job = tmp_path / "0_1"
    job.mkdir()
    (job / "results.dat").write_text("1 termination_reason\n")
    (job / "saddle.con").write_text("saddle\n")
    comm = Communicator(tmp_path, 1, config=SimpleNamespace(path_pot=""))
    size, bundled = comm.get_bundle_size(str(job))
    assert size == 1
    assert bundled is False
    kept = list(comm.unbundle(tmp_path, lambda name: name == "0_1"))
    assert len(kept) == 1
    slot = kept[0][0]
    assert slot["name"] == "0_1"
    assert slot["results.dat"].getvalue() == "1 termination_reason\n"
    assert slot["saddle.con"].getvalue() == "saddle\n"


def test_empty_client_directory_is_skipped(tmp_path):
    job = tmp_path / "0_2"
    job.mkdir()
    comm = Communicator(tmp_path, 1, config=SimpleNamespace(path_pot=""))
    assert list(comm.unbundle(tmp_path, lambda name: True)) == []


def test_make_bundles_writes_the_invariant_and_the_displacement(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    comm = Communicator(scratch, 1, config=SimpleNamespace(path_pot=""))
    paths = list(
        comm.make_bundles(
            [{"id": "0_3", "displacement.con": StringIO("disp\n")}],
            {"pos.con": (StringIO("pos\n"), 0o644)},
        )
    )
    job = Path(paths[0])
    assert (job / "pos.con").read_text() == "pos\n"
    assert (job / "displacement.con").read_text() == "disp\n"


def test_missing_local_client_raises(tmp_path):
    with pytest.raises(CommunicatorError, match="Can't find client"):
        Local(tmp_path, str(tmp_path / "missing-client"), 1, 1, config=SimpleNamespace())


def test_unknown_communicator_type_raises():
    reset_communicators()
    with pytest.raises(ValueError):
        get_communicator(SimpleNamespace(comm_type="nope"))


def test_local_queue_starts_empty(tmp_path):
    scratch = tmp_path / "scratch"
    results = tmp_path / "results"
    scratch.mkdir()
    results.mkdir()
    comm = Local(scratch, "/usr/bin/true", 1, 1, config=SimpleNamespace(debug_keep_all_results=False))
    assert comm.get_queue_size() == 0
    assert list(comm.get_results(results, lambda name: True)) == []
    job = tmp_path / "failed-job"
    job.mkdir()
    finished = SimpleNamespace(returncode=1, pid=1)
    assert comm.check_job((finished, str(job))) is None
    comm.cleanup()
