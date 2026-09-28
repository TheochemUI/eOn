"""atoms.match removes a rigid translation the way Matter::compare does."""
from __future__ import annotations

import numpy as np

from eon import atoms
from eon.structure import Structure


def _cell(fixed_first: bool = False) -> Structure:
    s = Structure(3)
    s.r = np.array([[1.0, 1.0, 1.0], [2.4, 1.1, 0.9], [1.2, 2.5, 1.3]])
    s.box = np.eye(3) * 6.0
    s.mass = np.ones(3)
    free = np.ones(3)
    if fixed_first:
        free[0] = 0.0
    s.free = free
    s.names = ["Si"] * 3
    return s


def _shifted(s: Structure) -> Structure:
    # A rigid drift that carries atoms across the periodic boundary.
    t = s.copy()
    t.r = np.mod(s.r + np.array([5.3, -0.7, 2.9]), 6.0)
    return t


def test_rigid_translation_matches_only_with_removal():
    a = _cell()
    b = _shifted(a)
    assert not atoms.match(a, b, 0.1, 3.3, False)
    assert atoms.match(a, b, 0.1, 3.3, False, remove_translation=True)
    assert atoms.match(a, b, 0.1, 3.3, True, use_identical=True, remove_translation=True)


def test_translation_is_kept_with_a_fixed_atom():
    a = _cell(fixed_first=True)
    b = _shifted(a)
    assert not atoms.match(a, b, 0.1, 3.3, False, remove_translation=True)


def test_a_different_configuration_still_differs():
    a = _cell()
    b = _shifted(a)
    b.r[2] += np.array([0.5, 0.0, 0.0])
    assert not atoms.match(a, b, 0.1, 3.3, False, remove_translation=True)
