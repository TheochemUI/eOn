"""Python CNA labels match EpiCenters.cpp (0 fcc, 1 hcp, 2 other)."""

import numpy as np
import pytest

pytest.importorskip("vesin")

from eon import atoms
from eon.structure import Structure


def _isolated_pair():
    p = Structure(2)
    p.r = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    p.box = np.eye(3) * 30.0
    p.names = ["Al", "Al"]
    return p


def _fcc_al(nrep=3, a=4.05):
    basis = np.array(
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
    )
    coords = [
        (np.array([i, j, k], dtype=float) + b) * a
        for i in range(nrep)
        for j in range(nrep)
        for k in range(nrep)
        for b in basis
    ]
    p = Structure(len(coords))
    p.r = np.asarray(coords)
    p.box = np.eye(3) * (nrep * a)
    p.names = ["Al"] * len(coords)
    return p


def test_isolated_atoms_are_other():
    labels = atoms.cna(_isolated_pair(), cutoff=3.0)
    assert set(labels.tolist()) == {atoms.CNA_OTHER}
    assert atoms.not_HCP_or_FCC(_isolated_pair(), cutoff=3.0) == [0, 1]


def test_periodic_fcc_is_label_zero():
    labels = atoms.cna(_fcc_al(), cutoff=3.2)
    assert int(np.min(labels)) == atoms.CNA_FCC
    assert int(np.max(labels)) == atoms.CNA_FCC
    assert atoms.not_HCP_or_FCC(_fcc_al(), cutoff=3.2) == []
