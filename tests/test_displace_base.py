"""Displacement helpers raise the builtin exception and stay finite."""

from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from eon.displace import Displace


class _Reactant:
    def __init__(self) -> None:
        self.r = np.zeros((2, 3))

    def atom_is_free(self):
        return np.array([True, True])

    def copy(self):
        other = _Reactant()
        other.r = self.r.copy()
        return other


def _displace() -> Displace:
    reactant = _Reactant()
    config = SimpleNamespace(
        disp_min_norm=0.0,
        displace_1d=False,
        random_mode=False,
        comp_brute_neighbors=False,
    )
    disp = Displace(reactant, 0.1, 3.0, None, config=config)
    disp.neighbors_list = [[1], [0]]
    return disp


def test_base_displace_raises_builtin_not_implemented():
    disp = Displace(
        SimpleNamespace(r=np.zeros((1, 3))),
        0.1,
        3.0,
        None,
        config=SimpleNamespace(),
    )
    with pytest.raises(NotImplementedError) as caught:
        disp.make_displacement()
    assert type(caught.value) is NotImplementedError


def test_list_epicenter_debug_log_stays_finite(monkeypatch):
    disp = _displace()
    disp.void_bias_fraction = 0.0
    log = logging.getLogger("displace")
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)
    monkeypatch.setattr(np.random, "normal", lambda scale, size: np.ones(size))
    try:
        atoms, mode = disp.get_displacement([0])
    finally:
        log.removeHandler(handler)
        log.setLevel(logging.WARNING)
    assert np.isfinite(atoms.r).all()
    assert np.isfinite(mode).all()
    assert np.isclose(np.linalg.norm(mode), 1.0)


def test_cancelled_void_bias_keeps_a_finite_mode(monkeypatch):
    disp = _displace()
    monkeypatch.setattr(
        "eon.displace.atoms.neighbor_list_vectors",
        lambda reactant, radius, brute: [np.zeros((1, 3)), np.zeros((1, 3))],
    )
    monkeypatch.setattr(np.random, "normal", lambda scale, size: np.ones(size))
    _atoms, mode = disp.get_displacement(0)
    assert np.isfinite(mode).all()
    assert np.isclose(np.linalg.norm(mode), 1.0)
