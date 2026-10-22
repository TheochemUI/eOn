"""A failed prng.pkl write must leave the previous state."""

from __future__ import annotations

import pickle

import pytest

from eon import fileio as io


def test_failed_prng_dump_keeps_previous_bytes(tmp_path, monkeypatch):
    path = tmp_path / "prng.pkl"
    io.save_prng_state(str(path))
    original = path.read_bytes()
    assert original

    def boom(*_args, **_kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(pickle, "dump", boom)
    with pytest.raises(RuntimeError, match="disk full"):
        io.save_prng_state(str(path))
    assert path.read_bytes() == original
    with path.open("rb") as fh:
        pickle.load(fh)
