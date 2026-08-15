"""POSCAR is the one non-readcon structure format: VASP 4 (kdb) and VASP 5 (movies)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.structure import Structure


def _cu_h() -> Structure:
    s = Structure(2)
    s.box = np.diag([10.0, 11.0, 12.0])
    s.names = ["Cu", "H"]
    s.r = np.array([[0.25, 0.0, 0.0], [1.5, 0.5, 0.25]])
    s.free = np.array([1.0, 0.0])
    s.mass = np.array([63.546, 1.008])
    return s


def _vasp4_text() -> str:
    # Line 6 is the counts: the kdb / VASP 4 layout.
    return (
        "Cu H\n"
        "1.0\n"
        "10.0 0.0 0.0\n"
        "0.0 11.0 0.0\n"
        "0.0 0.0 12.0\n"
        "1 1\n"
        "Selective Dynamics\n"
        "Cartesian\n"
        "0.25 0.0 0.0 T T T\n"
        "1.5 0.5 0.25 F F F\n"
    )


def test_tokens_are_ints_detects_counts_versus_species():
    assert io._tokens_are_ints(["1", "2", "3"])
    assert not io._tokens_are_ints(["Cu", "H"])
    assert not io._tokens_are_ints(["Cu", "1"])
    assert not io._tokens_are_ints([])
    assert not io._tokens_are_ints(["1.0"])


def test_loadposcar_reads_vasp4_kdb_layout(tmp_path: Path):
    path = tmp_path / "SADDLE_0"
    path.write_text(_vasp4_text())
    s = io.loadposcar(str(path))
    assert list(s.names) == ["Cu", "H"]
    np.testing.assert_allclose(s.r[0], [0.25, 0.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(s.r[1], [1.5, 0.5, 0.25], atol=1e-12)
    np.testing.assert_array_equal(s.free, [1.0, 0.0])
    np.testing.assert_allclose(s.box, np.diag([10.0, 11.0, 12.0]), atol=1e-12)


def test_saveposcar_loadposcar_roundtrip(tmp_path: Path):
    """eOn must be able to read the VASP 5 POSCAR it writes."""
    original = _cu_h()
    out = tmp_path / "movie.poscar"
    io.saveposcar(str(out), original)
    text = out.read_text().splitlines()
    # VASP 5: species names on line 6, counts on line 7.
    assert text[5].split() == ["Cu", "H"]
    assert text[6].split() == ["1", "1"]

    back = io.loadposcar(str(out))
    assert list(back.names) == ["Cu", "H"]
    np.testing.assert_allclose(back.r, original.r, atol=1e-12)
    np.testing.assert_array_equal(back.free, original.free)
    np.testing.assert_allclose(back.box, original.box, atol=1e-12)


def test_loadposcar_reads_vasp5_with_comment_on_line_1(tmp_path: Path):
    path = tmp_path / "POSCAR"
    path.write_text(
        "some comment that is not species\n"
        "1.0\n"
        "10.0 0.0 0.0\n"
        "0.0 10.0 0.0\n"
        "0.0 0.0 10.0\n"
        "Cu H\n"
        "1 1\n"
        "Cartesian\n"
        "0.0 0.0 0.0\n"
        "1.0 0.0 0.0\n"
    )
    s = io.loadposcar(str(path))
    assert list(s.names) == ["Cu", "H"]
    np.testing.assert_allclose(s.r[1], [1.0, 0.0, 0.0], atol=1e-12)


def test_loadposcars_reads_eon_written_movie(tmp_path: Path):
    a = _cu_h()
    b = _cu_h()
    b.r = b.r + 0.5
    movie = tmp_path / "dynamics.poscar"
    io.saveposcar(str(movie), a, w="w")
    io.saveposcar(str(movie), b, w="a")
    frames = io.loadposcars(str(movie))
    assert len(frames) == 2
    np.testing.assert_allclose(frames[0].r, a.r, atol=1e-12)
    np.testing.assert_allclose(frames[1].r, b.r, atol=1e-12)
    assert list(frames[0].names) == ["Cu", "H"]
    assert list(frames[1].names) == ["Cu", "H"]


def test_loadposcar_rejects_a_species_name_as_a_count(tmp_path: Path):
    """A truncated VASP 5 file (species line, no counts) must not int('Cu')."""
    path = tmp_path / "bad.poscar"
    path.write_text(
        "Cu H\n"
        "1.0\n"
        "10.0 0.0 0.0\n"
        "0.0 10.0 0.0\n"
        "0.0 0.0 10.0\n"
        "Cu H\n"
    )
    with pytest.raises((ValueError, IndexError)):
        io.loadposcar(str(path))
