"""displace_atom_list = None is an empty list, and -1 stays -1."""

from __future__ import annotations

import pytest

from eon.config import ConfigClass


def _load(tmp_path, atom_list):
    path = tmp_path / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "[Saddle Search]",
                "displace_listed_atom_weight = 1",
                "displace_atom_list = %s" % atom_list,
                "",
            ]
        )
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def test_none_atom_list_is_empty(tmp_path):
    cfg = _load(tmp_path, "None")
    assert cfg.disp_listed_atoms == []


def test_negative_one_atom_list_is_kept(tmp_path):
    cfg = _load(tmp_path, "-1")
    assert cfg.disp_listed_atoms == [-1]


def test_comma_atom_list_is_ints(tmp_path):
    cfg = _load(tmp_path, "1, 2")
    assert cfg.disp_listed_atoms == [1, 2]
