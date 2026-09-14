"""prng.pkl must follow ConfigClass.path_root, not the process CWD."""
from __future__ import annotations

import os

import numpy as np

from eon import fileio as io


def test_save_prng_state_writes_given_path_not_cwd(tmp_path, monkeypatch):
    cwd = tmp_path / "cwd"
    root = tmp_path / "sim"
    cwd.mkdir()
    root.mkdir()
    monkeypatch.chdir(cwd)

    dest = io.prng_state_path(type("Cfg", (), {"path_root": str(root)})())
    np.random.seed(42)
    expected = np.random.random()
    np.random.seed(42)
    io.save_prng_state(dest)

    assert not (cwd / "prng.pkl").exists()
    assert os.path.isfile(dest)

    np.random.seed(0)
    assert np.random.random() != expected
    io.get_prng_state(dest)
    assert np.random.random() == expected


def test_prng_state_path_uses_path_root(tmp_path):
    cfg = type("Cfg", (), {"path_root": str(tmp_path / "root")})()
    assert io.prng_state_path(cfg) == os.path.join(str(tmp_path / "root"), "prng.pkl")
