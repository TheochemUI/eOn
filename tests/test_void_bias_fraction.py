"""A configured void bias of zero leaves the Gaussian mode unchanged."""

from types import SimpleNamespace

import numpy as np

from eon.displace import Displace


class _Reactant:
    def __init__(self):
        self.r = np.zeros((2, 3))

    def atom_is_free(self):
        return np.array([True, True])

    def copy(self):
        other = _Reactant()
        other.r = self.r.copy()
        return other


def test_zero_void_bias_keeps_the_gaussian_mode(monkeypatch):
    reactant = _Reactant()
    config = SimpleNamespace(
        disp_min_norm=0.0,
        displace_1d=False,
        random_mode=False,
        comp_brute_neighbors=False,
        void_bias_fraction=0.0,
    )
    disp = Displace(reactant, 0.1, 3.0, None, config=config)
    disp.neighbors_list = [[1], [0]]
    monkeypatch.setattr(
        "eon.displace.atoms.neighbor_list_vectors",
        lambda reactant, radius, brute: [
            [np.array([1.0, 0.0, 0.0])],
            [np.array([-1.0, 0.0, 0.0])],
        ],
    )
    monkeypatch.setattr(np.random, "normal", lambda scale, size: np.ones(size))
    _atoms, mode = disp.get_displacement(0)
    unbiased = np.zeros_like(reactant.r)
    unbiased[0] = np.ones(3)
    unbiased[1] = np.ones(3)
    unbiased /= np.linalg.norm(unbiased)
    assert np.allclose(mode, unbiased)
