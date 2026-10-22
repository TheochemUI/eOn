"""A killed recycling metadata rewrite leaves the previous text."""

from types import SimpleNamespace

import pytest

from eon.recycling import Recycling, SB_Recycling


def _boom(_src, _dst):
    raise OSError("killed")


def test_killed_superbasin_metadata_keeps_the_previous_text(tmp_path, monkeypatch):
    path = tmp_path / "recycling_data.txt"
    path.write_text("OLD\n")
    recycler = SB_Recycling.__new__(SB_Recycling)
    recycler.path = str(tmp_path)
    recycler.sb_state_nums = [0, 1]
    recycler.in_progress = False
    monkeypatch.setattr("eon.fileio.os.replace", _boom)
    with pytest.raises(OSError, match="killed"):
        recycler.write_metadata()
    assert path.read_text() == "OLD\n"


def test_killed_recycling_metadata_keeps_the_previous_text(tmp_path, monkeypatch):
    path = tmp_path / "recycling_metadata.txt"
    path.write_text("OLD\n")
    recycler = Recycling.__new__(Recycling)
    recycler.metadata_path = str(path)
    recycler.ref_state = SimpleNamespace(number=3)
    recycler.process_number = 1
    recycler.num_procs = 2
    recycler.moved = [0]
    recycler.unmoved = [1]
    recycler.process_atoms = [0]
    monkeypatch.setattr("eon.fileio.os.replace", _boom)
    with pytest.raises(OSError, match="killed"):
        recycler.write_recycling_metadata()
    assert path.read_text() == "OLD\n"
