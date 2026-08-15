"""Golden / integration tests: eon.fileio through live readcon on fixtures.

Proves loadcon/savecon/loadcons call the readcon package (not a private parser)
and that geometry / species / free flags / box round-trip correctly.
"""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import numpy as np
import pytest
import readcon

from eon import fileio as io

DATA = Path(__file__).resolve().parent / "data"


def _all_con_fixtures():
    paths = sorted(DATA.rglob("*.con"))
    assert paths, "expected tests/data/**/*.con fixtures"
    return paths


def _readable_con_fixtures():
    """Fixtures that the live readcon package can parse (skip known bad NEB files)."""
    good = []
    for path in _all_con_fixtures():
        try:
            frames = readcon.read_con(str(path))
            if frames:
                good.append(path)
        except OSError:
            continue
    assert good, "no readable .con fixtures under tests/data"
    return good


def test_fileio_imports_and_calls_readcon():
    """Structural: the shipped .con path binds the real readcon APIs (not a private parser).

    Frame construction lives in eon.structure and the file-level readers and
    writers in eon.fileio, so the pair is checked together: asserting every
    call sits in one module would only pin where the split happens to fall.
    """
    from eon import structure as structmod

    fileio_src = Path(io.__file__).read_text(encoding="utf-8")
    structure_src = Path(structmod.__file__).read_text(encoding="utf-8")
    src = fileio_src + structure_src

    assert "import readcon" in fileio_src
    assert "import readcon" in structure_src
    for api in (
        "readcon.read_con",
        "readcon.read_con_string",
        "readcon.write_con",
        "readcon.write_con_string",
        "readcon.Atom",
        "readcon.ConFrame",
    ):
        assert api in src, f"the .con path must use {api}"
    # .con stays on the core reader; chemfiles is the XYZ/PDB/GRO edge.
    assert "readcon.read_con" in fileio_src
    assert "readcon.read_chemfiles_first" in fileio_src
    # Not a hand-rolled tokenizer loop for coordinates
    assert "Coordinates of component" not in src


def _version_tuple(ver: str) -> tuple[int, int, int]:
    """(major, minor, patch) from a version string, ignoring any pre/post suffix."""
    parts = []
    for chunk in ver.split(".")[:3]:
        digits = re.match(r"\d+", chunk)
        parts.append(int(digits.group()) if digits else 0)
    parts += [0] * (3 - len(parts))
    return tuple(parts)


def test_readcon_version_and_chemfiles_status():
    """Record readcon pin; eOn server I/O does not require chemfiles extra."""
    ver = getattr(readcon, "__version__", None)
    assert ver is not None
    # eOn pins readcon >=0.14,<0.15 in pyproject.toml, pyproject-pyeonclient.toml
    # and pixi.toml, matching subprojects/readcon-core.wrap v0.14.0. Client and
    # server must share CON spec 3.
    parts = _version_tuple(ver)
    assert parts >= (0, 14, 0), ver
    assert parts < (0, 15, 0), ver
    has_cf = readcon.has_chemfiles_support()
    # Default install (no [chemfiles] extra): False. Document either way.
    assert isinstance(has_cf, bool)
    # Server path must work without chemfiles
    assert callable(readcon.read_con)
    assert callable(readcon.write_con)


def test_readcon_core_not_chemfiles_for_server_io():
    """The .con path uses readcon-core even when chemfiles is installed."""
    src = inspect.getsource(io.loadcon)
    assert "read_chemfiles" not in src
    fixture = DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con"
    frames = readcon.read_con(str(fixture))
    assert len(frames) >= 1
    assert len(frames[0]) > 0


def test_loadxyz_requires_chemfiles_extra():
    """Without the chemfiles extra, loadxyz names the missing extra."""
    if readcon.has_chemfiles_support():
        pytest.skip("readcon-chemfiles is installed")
    with pytest.raises(RuntimeError, match="readcon-chemfiles"):
        io.loadxyz("water.xyz")


@pytest.mark.parametrize("path", _readable_con_fixtures(), ids=lambda p: str(p.relative_to(DATA)))
def test_loadcon_matches_live_readcon_frame(path: Path):
    """Primary observable: eon.fileio.loadcon == readcon frame geometry."""
    frames = readcon.read_con(str(path))
    assert frames
    ref = frames[0]
    atoms = io.loadcon(str(path))

    assert len(atoms) == len(ref)
    # Species
    assert list(atoms.names) == [a.symbol for a in ref.atoms]
    # Coordinates
    ref_r = np.array([[a.x, a.y, a.z] for a in ref.atoms], dtype=float)
    np.testing.assert_allclose(atoms.r, ref_r, rtol=0, atol=1e-9)
    # Mass
    ref_mass = np.array(
        [a.mass if a.mass is not None else 0.0 for a in ref.atoms], dtype=float
    )
    np.testing.assert_allclose(atoms.mass, ref_mass, rtol=0, atol=1e-9)
    # Free / fixed: per-axis, 0.0 if that CON flag is set
    for i, a in enumerate(ref.atoms):
        expect = [0.0 if flag else 1.0 for flag in a.fixed]
        np.testing.assert_array_equal(atoms.free[i], expect)
    # Box from cell lengths + angles (alpha, beta, gamma)
    from eon.geometry.cell import box_to_length_angle

    box_lens, box_angs = box_to_length_angle(atoms.box)
    np.testing.assert_allclose(box_lens, np.array(list(ref.cell), dtype=float), rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(box_angs, np.array(list(ref.angles), dtype=float), rtol=1e-6, atol=1e-6)


def test_loadcon_server_fixture_first_atom_content():
    """Content check on a known fixture (not merely 'ran')."""
    path = DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con"
    atoms = io.loadcon(str(path))
    # Cross-check first atom against live readcon (source of truth for the file)
    ref = readcon.read_con(str(path))[0].atoms[0]
    assert len(atoms) == 343
    assert atoms.names[0] == "Pt" == ref.symbol
    np.testing.assert_allclose(
        atoms.r[0], [ref.x, ref.y, ref.z], rtol=0, atol=1e-9
    )
    # Explicit content from fixture (readcon-parsed): first Pt near heptamer layer
    assert atoms.r[0, 2] > 14.0
    assert float(atoms.mass[0]) == pytest.approx(195.084, abs=1e-6)


@pytest.mark.parametrize(
    "path",
    [
        DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con",
        DATA / "client" / "oxadiazole" / "pos.con",
        DATA / "client" / "mini_1Pt_oneLayer" / "pos.con",
        DATA / "client" / "one_Pt_on_frozenSurface" / "pos.con",
    ],
    ids=lambda p: p.parent.name + "/" + p.name,
)
def test_loadcon_savecon_roundtrip(path: Path, tmp_path: Path):
    """Round-trip positions, species, free flags, box via shipped savecon/loadcon."""
    original = io.loadcon(str(path))
    out = tmp_path / "roundtrip.con"
    io.savecon(str(out), original)
    again = io.loadcon(str(out))

    assert len(again) == len(original)
    assert list(again.names) == list(original.names)
    np.testing.assert_allclose(again.r, original.r, rtol=0, atol=1e-6)
    np.testing.assert_allclose(again.mass, original.mass, rtol=0, atol=1e-6)
    np.testing.assert_array_equal(again.free.astype(int), original.free.astype(int))
    np.testing.assert_allclose(again.box, original.box, rtol=1e-6, atol=1e-6)
    np.testing.assert_array_equal(again.free, original.free)


def test_loadcons_multi_frame_via_readcon(tmp_path: Path):
    """loadcons returns one Atoms per frame written with multi-frame write_con."""
    path = DATA / "client" / "oxadiazole" / "pos.con"
    a = io.loadcon(str(path))
    b = a.copy()
    b.r = b.r + 0.1
    multi = tmp_path / "multi.con"
    # Build multi-frame file through readcon (same path savecon append uses)
    f1 = io._atoms_to_frame(a)
    f2 = io._atoms_to_frame(b)
    readcon.write_con(str(multi), [f1, f2])
    frames = io.loadcons(str(multi))
    assert len(frames) == 2
    assert len(frames[0]) == len(a)
    np.testing.assert_allclose(frames[0].r, a.r, atol=1e-6)
    np.testing.assert_allclose(frames[1].r, b.r, atol=1e-6)


def test_savecon_append_mode(tmp_path: Path):
    path = DATA / "client" / "oxadiazole" / "pos.con"
    a = io.loadcon(str(path))
    out = tmp_path / "append.con"
    io.savecon(str(out), a, w="w")
    shifted = a.copy()
    shifted.r = shifted.r + 0.25
    io.savecon(str(out), shifted, w="a")
    frames = io.loadcons(str(out))
    assert len(frames) == 2
    np.testing.assert_allclose(frames[1].r, shifted.r, atol=1e-6)


def test_savecon_append_does_not_reread_the_movie(tmp_path: Path, monkeypatch):
    """Append must serialize the new frame only, not parse N existing frames."""
    path = DATA / "client" / "oxadiazole" / "pos.con"
    a = io.loadcon(str(path))
    out = tmp_path / "append.con"
    io.savecon(str(out), a, w="w")
    before = out.read_bytes()

    def boom(*_args, **_kwargs):
        raise AssertionError("append must not re-parse the movie")

    monkeypatch.setattr(readcon, "read_con", boom)
    shifted = a.copy()
    shifted.r = shifted.r + 0.25
    io.savecon(str(out), shifted, w="a")
    monkeypatch.undo()

    after = out.read_bytes()
    assert after.startswith(before)
    assert len(after) > len(before)
    frames = io.loadcons(str(out))
    assert len(frames) == 2
    np.testing.assert_allclose(frames[1].r, shifted.r, atol=1e-6)
