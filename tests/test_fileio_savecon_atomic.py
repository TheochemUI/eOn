"""A failed savecon rewrite leaves the previous .con."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.structure import Structure


def _one():
    structure = Structure(1)
    structure.r[0] = [0.0, 0.0, 0.0]
    structure.names[0] = "H"
    structure.mass[0] = 1.0
    structure.box = np.eye(3) * 10.0
    return structure


def test_failed_savecon_keeps_previous_con(tmp_path, monkeypatch):
    dest = tmp_path / "pos.con"
    dest.write_text("ORIGINAL\n")

    def boom(path, frames):
        Path(path).write_text("partial\n")
        raise RuntimeError("disk full")

    monkeypatch.setattr(io.readcon, "write_con", boom)
    with pytest.raises(RuntimeError, match="disk full"):
        io.savecon(str(dest), _one())
    assert dest.read_text() == "ORIGINAL\n"
