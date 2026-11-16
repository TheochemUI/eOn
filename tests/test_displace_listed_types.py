"""displace_all_listed covers a type list, not only an atom list."""

from types import SimpleNamespace

import numpy as np

from eon.displace import DisplacementManager
from eon.structure import Structure


def _reactant():
    atoms = Structure(3)
    atoms.names = ["Al", "Al", "Pt"]
    atoms.mass[:] = 26.982
    atoms.box = np.diag([30.0, 30.0, 30.0])
    atoms.r[:] = [
        [0.0, 0.0, 0.0],
        [12.0, 0.0, 0.0],
        [0.0, 12.0, 0.0],
    ]
    return atoms


def _config():
    return SimpleNamespace(
        disp_magnitude=0.1,
        disp_radius=0.0,
        disp_min_norm=0.0,
        displace_1d=False,
        random_mode=False,
        comp_brute_neighbors=False,
        comp_neighbor_cutoff=3.5,
        comp_use_covalent=False,
        comp_covalent_scale=1.3,
        disp_max_coord=12,
        displace_all_listed=True,
        disp_listed_atoms=[],
        disp_listed_types=["Al"],
        displace_random_weight=0.0,
        displace_listed_atom_weight=0.0,
        displace_listed_type_weight=1.0,
        displace_under_coordinated_weight=0.0,
        displace_least_coordinated_weight=0.0,
        displace_not_FCC_HCP_weight=0.0,
        displace_not_TCP_BCC_weight=0.0,
        displace_not_TCP_weight=0.0,
        displace_water_weight=0.0,
    )


def test_all_listed_types_move(monkeypatch):
    monkeypatch.setattr(np.random, "random_sample", lambda: 0.5)
    monkeypatch.setattr(np.random, "randint", lambda *args, **kwargs: 0)
    monkeypatch.setattr(np.random, "normal", lambda scale, size: np.ones(size))
    manager = DisplacementManager(_reactant(), None, config=_config())
    _atoms, mode = manager.make_displacement()
    assert np.linalg.norm(mode[0]) > 0.0
    assert np.linalg.norm(mode[1]) > 0.0
    assert np.linalg.norm(mode[2]) == 0.0
