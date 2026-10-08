"""Each displacement kind returns a finite mode on a small cluster."""

from types import SimpleNamespace

import numpy as np

from eon.displace import (
    DisplacementManager,
    Leastcoordinated,
    ListedAtoms,
    ListedTypes,
    NotFCCorHCP,
    NotTCP,
    NotTCPorBCC,
    Random,
    Undercoordinated,
    Water,
)
from eon.structure import Structure


def _cluster():
    atoms = Structure(4)
    atoms.names = ["Pt", "Pt", "Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[:] = [
        [0.0, 0.0, 0.0],
        [2.6, 0.0, 0.0],
        [1.3, 2.2, 0.0],
        [1.3, 0.8, 2.1],
    ]
    return atoms


def _cfg():
    return SimpleNamespace(
        disp_magnitude=0.1,
        disp_radius=4.0,
        disp_min_norm=0.0,
        displace_1d=False,
        random_mode=False,
        comp_brute_neighbors=False,
        comp_neighbor_cutoff=3.5,
        comp_use_covalent=False,
        comp_covalent_scale=1.3,
        disp_max_coord=12,
        displace_all_listed=False,
        disp_listed_atoms=[0, 1],
        disp_listed_types=["Pt"],
        displace_random_weight=1.0,
        displace_listed_atom_weight=0.0,
        displace_listed_type_weight=0.0,
        displace_under_coordinated_weight=0.0,
        displace_least_coordinated_weight=0.0,
        displace_not_FCC_HCP_weight=0.0,
        displace_not_TCP_BCC_weight=0.0,
        displace_not_TCP_weight=0.0,
        displace_water_weight=0.0,
    )


def _finite(pair):
    atoms, mode = pair
    assert np.isfinite(atoms.r).all()
    assert np.isfinite(mode).all()
    assert np.isclose(np.linalg.norm(mode), 1.0)


def test_each_epicenter_kind_is_finite():
    np.random.seed(1)
    reactant = _cluster()
    cfg = _cfg()
    kinds = [
        Random(reactant, 0.1, 4.0, config=cfg),
        ListedAtoms(reactant, 0.1, 4.0, config=cfg),
        ListedTypes(reactant, 0.1, 4.0, config=cfg),
        Undercoordinated(reactant, 12, 0.1, 4.0, cutoff=3.5, config=cfg),
        Leastcoordinated(reactant, 0.1, 4.0, cutoff=3.5, config=cfg),
        NotFCCorHCP(reactant, 0.1, 4.0, cutoff=3.5, config=cfg),
        NotTCPorBCC(reactant, 0.1, 4.0, cutoff=3.5, config=cfg),
        NotTCP(reactant, 0.1, 4.0, cutoff=3.5, config=cfg),
    ]
    for kind in kinds:
        _finite(kind.make_displacement())


def test_manager_random_weight_is_finite():
    np.random.seed(2)
    reactant = _cluster()
    _atoms, mode = DisplacementManager(reactant, None, _cfg()).make_displacement()
    assert np.isfinite(mode).all()
    assert np.isclose(np.linalg.norm(mode), 1.0)


def test_water_rigid_shift_stays_finite():
    np.random.seed(3)
    water = Structure(3)
    water.names = ["H", "H", "O"]
    water.mass[:] = [1.0, 1.0, 16.0]
    water.box = np.diag([15.0, 15.0, 15.0])
    water.r[:] = [[0.0, 0.8, 0.0], [0.0, -0.8, 0.0], [0.0, 0.0, 0.1]]
    moved, displacement = Water(water, 0.05, 0.1).make_displacement()
    assert np.isfinite(moved.r).all()
    assert np.isfinite(displacement).all()
    # The three atoms share one translation, so the displacement is not zero.
    assert np.linalg.norm(displacement) > 0.0
