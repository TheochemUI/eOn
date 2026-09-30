"""eon.atoms loads and compares structures when minimage is not installed.

minimage is optional: the pbc helper falls back to numpy, and only the
readcon_ops callers need the Cell.wrap_many binding. Importing the module
has to work without it, since the aKMC server imports it at start-up.
"""

import importlib
import sys

import numpy as np


def test_atoms_imports_and_matches_without_minimage(monkeypatch):
    monkeypatch.setitem(sys.modules, "minimage", None)  # import raises
    monkeypatch.delitem(sys.modules, "eon.atoms", raising=False)
    atoms = importlib.import_module("eon.atoms")
    from eon.structure import Structure

    a = Structure(2)
    a.r = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    a.box = np.eye(3) * 10.0
    a.names = ["Al", "Al"]
    a.mass = np.array([26.98, 26.98])
    b = a.copy()
    b.r = a.r[::-1].copy()  # same structure, atoms swapped
    assert atoms.identical(a, b, 0.1)
    b.r[0, 0] += 0.5
    assert not atoms.identical(a, b, 0.1)
