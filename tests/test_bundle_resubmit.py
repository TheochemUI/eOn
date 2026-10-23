"""Submitting the same bundle id again replaces the job directory."""

from io import StringIO
from pathlib import Path
from types import SimpleNamespace

from eon.communicator import Communicator


def test_second_bundle_of_the_same_id_replaces_the_job(tmp_path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    comm = Communicator(scratch, 1, config=SimpleNamespace(path_pot=""))
    first = {"id": "job 1", "displacement.con": StringIO("first\n")}
    list(comm.make_bundles([first], {"pos.con": (StringIO("pos\n"), 0o644)}))
    second = {"id": "job 1", "displacement.con": StringIO("second\n")}
    paths = list(comm.make_bundles([second], {"pos.con": (StringIO("pos\n"), 0o644)}))
    job = Path(paths[0])
    assert job.name == "job 1"
    assert (job / "displacement.con").read_text() == "second\n"
    assert (job / "pos.con").read_text() == "pos\n"
