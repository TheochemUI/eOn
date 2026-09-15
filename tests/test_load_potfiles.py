"""load_potfiles skips subdirectories of pot_dir, not CWD names."""

import pytest

pytest.importorskip("readcon")

from eon.fileio import load_potfiles


def test_load_potfiles_skips_subdir_inside_pot_dir(tmp_path):
    pot = tmp_path / "potfiles"
    pot.mkdir()
    (pot / "POTCAR").write_text("ok")
    (pot / "nested").mkdir()
    (pot / "nested" / "x").write_text("no")
    got = load_potfiles(str(pot))
    assert set(got) == {"POTCAR"}
    data, _mode = got["POTCAR"]
    assert data.getvalue() == "ok"


def test_load_potfiles_missing_dir(tmp_path):
    assert load_potfiles(str(tmp_path / "nope")) == {}
