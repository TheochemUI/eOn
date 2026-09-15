"""rot_match Kabsch fallback: a rigid translation is the same structure."""

import numpy as np
import pytest

pytest.importorskip("vesin")

from eon.atoms import rot_match
from eon.structure import Structure


def _triangle():
    p = Structure(3)
    p.r = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    p.box = np.eye(3) * 20.0
    p.names = ["H", "H", "H"]
    return p


def test_translated_copy_matches():
    a = _triangle()
    b = a.copy()
    b.r = a.r + np.array([2.5, -1.0, 0.4])
    assert rot_match(a, b, 1e-8)


def test_stretched_copy_does_not_match():
    a = _triangle()
    b = a.copy()
    b.r = a.r * 1.5
    assert not rot_match(a, b, 0.05)
