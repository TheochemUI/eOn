"""Client and server CON bytes against the restart golden."""

from __future__ import annotations

import inspect
import math
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon import process_catalog
from eon.concorpus import (
    corpus_buffered_frames,
    flush_corpus_appends,
    note_corpus_append,
    traj_id,
)
from eon.structure import Structure

_PI = 3.14159265358979323846
_BETA = float.fromhex("0x1.2c00000000001p+6")
_GOLDEN = (
    Path(__file__).resolve().parents[1]
    / "client"
    / "unit_tests"
    / "data"
    / "restart"
)
_PATH_ID = 14164577628502097180


def _subnormal() -> float:
    return float(np.array([1], dtype=np.uint64).view(np.float64)[0])


def _restart_box() -> np.ndarray:
    alpha = 80.0 * _PI / 180.0
    beta = _BETA * _PI / 180.0
    gamma = 100.0 * _PI / 180.0
    box = np.zeros((3, 3), dtype=float)
    box[0, 0] = 1.0
    box[1, 0] = math.cos(gamma)
    box[1, 1] = math.sin(gamma)
    box[2, 0] = math.cos(beta)
    box[2, 1] = (math.cos(alpha) - box[1, 0] * box[2, 0]) / box[1, 1]
    box[2, 2] = math.sqrt(1.0 - box[2, 0] ** 2 - box[2, 1] ** 2)
    box[0, :] *= 10.0
    box[1, :] *= 11.0
    box[2, :] *= 12.0
    return box


def _restart_structure(*, velocities: bool) -> Structure:
    sub = _subnormal()
    structure = Structure(4)
    structure.box = _restart_box()
    structure.names = ["Cu", "H", "Cu", "H"]
    structure.mass[:] = [63.546, 1.008, 63.546, 1.008]
    structure.r[:, 0] = np.array([-0.0, 1.0, 2.0, 3.0])
    structure.r[:, 1] = [0.5, 1.5, 2.5, 3.5]
    structure.r[0, 2] = sub
    structure.r[1:, 2] = [1.5, 2.5, 3.5]
    structure.free[0, 0] = 0.0
    structure.con_energy = -1.25
    if velocities:
        structure.velocities = np.array(
            [
                [0.0, -0.0, sub],
                [0.01, -0.02, 0.03],
                [-0.04, 0.05, -0.06],
                [0.07, 0.08, -0.09],
            ],
            dtype=float,
        )
    else:
        structure.forces = np.array(
            [
                [0.0, 0.1, -0.2],
                [0.3, -0.0, 0.4],
                [-0.5, 0.2, 0.1],
                [0.05, -0.3, 0.25],
            ],
            dtype=float,
        )
    return structure


def _u64(values: np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def test_savecon_matches_restart_golden(tmp_path: Path):
    golden = _GOLDEN / "restart_exact.con"
    out = tmp_path / "restart_exact.con"
    io.savecon(str(out), _restart_structure(velocities=False))
    assert out.read_bytes() == golden.read_bytes()


def test_savecon_matches_velocity_golden(tmp_path: Path):
    golden = _GOLDEN / "restart_exact_vel.con"
    out = tmp_path / "restart_exact_vel.con"
    io.savecon(str(out), _restart_structure(velocities=True))
    assert out.read_bytes() == golden.read_bytes()


def test_from_conframe_xyz_keeps_signed_zero_and_subnormal():
    import readcon

    frame = readcon.read_con(str(_GOLDEN / "restart_exact.con"))[0]
    structure = Structure.from_conframe(frame)
    got = _u64(structure.r)
    expect = _u64(_restart_structure(velocities=False).r)
    assert got.tolist() == expect.tolist()
    src = inspect.getsource(Structure.from_conframe)
    assert "frame.xyz" in src
    assert "coords_array" not in src


def test_loaded_golden_rewrites_the_same_bytes(tmp_path: Path):
    golden = (_GOLDEN / "restart_exact.con").read_bytes()
    structure = io.loadcon(str(_GOLDEN / "restart_exact.con"))
    out = tmp_path / "again.con"
    io.savecon(str(out), structure)
    assert out.read_bytes() == golden


def test_two_cu_masses_collapse_to_the_first(caplog: pytest.LogCaptureFixture):
    structure = Structure(2)
    structure.names = ["Cu", "Cu"]
    structure.mass[:] = [63.5, 62.0]
    structure.box = np.eye(3) * 10.0
    with caplog.at_level("INFO", logger="eon.structure"):
        frame = structure.to_conframe()
    assert frame.atoms[0].mass == 63.5
    assert frame.atoms[1].mass == 63.5
    assert any("Collapsed mass" in rec.message for rec in caplog.records)


def test_traj_id_is_the_path_only():
    path = "/var/eon/restart_exact.con"
    assert traj_id(path) == _PATH_ID
    assert traj_id(path, text="ignored") == _PATH_ID


def test_corpus_appends_flush_in_groups_of_64(tmp_path: Path, monkeypatch):
    calls: list[tuple[str, list[str]]] = []

    def _fake(resolved: str, frames: list[str]) -> None:
        calls.append((resolved, list(frames)))

    flush_corpus_appends()
    monkeypatch.setattr("eon.concorpus._extend_batch", _fake)
    path = tmp_path / "movie.con"
    for _ in range(63):
        note_corpus_append(path, "frame\n")
    assert corpus_buffered_frames(path) == 63
    assert calls == []
    note_corpus_append(path, "frame\n")
    assert corpus_buffered_frames(path) == 0
    assert len(calls) == 1
    assert len(calls[0][1]) == 64
    note_corpus_append(path, "frame\n")
    assert corpus_buffered_frames(path) == 1
    flush_corpus_appends()
    assert corpus_buffered_frames(path) == 0
    assert len(calls) == 2
    assert len(calls[1][1]) == 1


def test_reading_the_corpus_flushes(tmp_path: Path, monkeypatch):
    from eon.concorpus import stored_frame_text

    calls: list[int] = []

    def _fake(_resolved: str, frames: list[str]) -> None:
        calls.append(len(frames))

    flush_corpus_appends()
    monkeypatch.setattr("eon.concorpus._extend_batch", _fake)
    path = tmp_path / "one.con"
    note_corpus_append(path, "frame\n")
    assert corpus_buffered_frames(path) == 1
    stored_frame_text(path)
    assert calls == [1]
    assert corpus_buffered_frames(path) == 0


def test_catalog_omits_precision_17():
    src = inspect.getsource(process_catalog._saddle_text_with_mode)
    assert "write_con_string([frame])" in src
    assert ", 17)" not in src


def test_append_once_then_read(tmp_path: Path):
    import readcon

    structure = Structure(1)
    structure.names = ["Cu"]
    structure.mass[:] = [63.546]
    structure.box = np.eye(3) * 10.0
    structure.r[0] = [1.0, 2.0, 3.0]
    out = tmp_path / "movie.con"
    io.savecon(str(out), structure, w="w")
    structure.r[0, 0] = 1.5
    io.savecon(str(out), structure, w="a")
    frames = readcon.read_con(str(out))
    assert len(frames) == 2
    assert frames[1].atoms[0].x == pytest.approx(1.5)
