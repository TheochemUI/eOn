"""Bundle slot names reconstruct consecutive wuids from the first job id."""

import pytest

from eon.communicator import bundled_job_name


def test_bundled_job_name_slots():
    assert bundled_job_name("12_40", 0) == "12_40"
    assert bundled_job_name("12_40", 2) == "12_42"


def test_bundled_job_name_passthrough():
    assert bundled_job_name("oddname", 1) == "oddname"


def test_remove_tree_and_empty_parents(tmp_path):
    pytest.importorskip("readcon")
    from eon.fileio import remove_tree_and_empty_parents

    leaf = tmp_path / "a" / "b" / "c"
    leaf.mkdir(parents=True)
    (leaf / "x").write_text("1")
    remove_tree_and_empty_parents(str(leaf))
    assert not leaf.exists()
