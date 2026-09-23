"""Con path I/O stores the frame in a readcon-db corpus."""

from __future__ import annotations

from pathlib import Path

import pytest
import readcon

readcon_db = pytest.importorskip("readcon_db")

from eon import fileio as io
from eon.concorpus import corpus_dir, stored_frame_text


def _one_frame_fixture() -> Path:
    for path in sorted(Path("tests/data").rglob("*.con")):
        try:
            frames = readcon.read_con(str(path))
        except OSError:
            continue
        if len(frames) == 1:
            return path
    pytest.skip("no single-frame con fixture")


def _stored(path: Path) -> str:
    stored = stored_frame_text(path)
    assert stored
    return stored


def test_savecon_mirrors_into_readcon_db(tmp_path: Path):
    fixture = _one_frame_fixture()
    (tmp_path / "config.ini").write_text("[Main]\n")
    src = tmp_path / "in.con"
    src.write_text(fixture.read_text())
    atoms = io.loadcon(str(src))
    dest = tmp_path / "reactant.con"
    io.savecon(str(dest), atoms)

    stored = _stored(dest)
    original = readcon.read_con(str(dest))[0]
    mirrored = readcon.read_con_string(stored)[0]
    assert len(mirrored.atoms) == len(original.atoms)
    assert corpus_dir(dest) == tmp_path / "readcon.db"


def test_loadcon_mirrors_an_existing_file(tmp_path: Path):
    fixture = _one_frame_fixture()
    dest = tmp_path / "copied.con"
    dest.write_text(fixture.read_text())
    io.loadcon(str(dest))
    stored = _stored(dest)
    assert len(readcon.read_con_string(stored)[0].atoms) == len(
        readcon.read_con(str(dest))[0].atoms
    )
