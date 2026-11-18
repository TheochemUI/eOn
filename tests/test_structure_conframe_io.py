"""``Structure`` <-> ``readcon.ConFrame`` conversion on cells of real size.

Every ``.con`` the server reads or writes goes through these two methods.
The tests hold the output to the per-atom definitions in
``driver_oracles`` and bound the interpreter work: it must not grow with
the number of atoms, because a few thousand atoms are parsed per search.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

readcon = pytest.importorskip("readcon")

from driver_oracles import (  # noqa: E402
    count_events,
    ref_con_text,
    ref_from_conframe,
    total,
)

from eon.structure import Structure  # noqa: E402


def make_structure(n, *, forces=False, velocities=False, shuffle=True, seed=3):
    """Two interleaved species, ids permuted, a third of the atoms frozen."""
    rng = np.random.RandomState(seed)
    s = Structure(n)
    s.box = np.diag([20.0, 21.0, 22.0])
    s.r = rng.uniform(0.0, 20.0, size=(n, 3))
    s.names = ["Cu" if i % 2 == 0 else "Ag" for i in range(n)]
    s.mass = np.array([63.546 if i % 2 == 0 else 107.8682 for i in range(n)])
    free = np.ones((n, 3))
    free[::3] = 0.0
    free[1::7, 1] = 0.0
    s.free = free
    if shuffle:
        s.atom_ids = (rng.permutation(n) + 1).astype(np.uint64)
    if forces:
        s.forces = rng.normal(size=(n, 3))
    if velocities:
        s.velocities = rng.normal(size=(n, 3))
    s.con_energy = -12.5
    return s


def reorder(s):
    """The Structure a write and read of *s* gives: ascending ids."""
    order = np.argsort(s.atom_ids, kind="stable")
    t = s.copy()
    t.r = s.r[order]
    t.free = s.free[order]
    t.mass = s.mass[order]
    t.atom_ids = s.atom_ids[order]
    t.names = [s.names[j] for j in order]
    return t


def assert_same_structure(a, b):
    assert len(a) == len(b)
    assert a.names == b.names
    for field in ("r", "box", "mass", "atom_ids", "file_to_struct"):
        x, y = getattr(a, field), getattr(b, field)
        assert x.dtype == y.dtype, field
        assert np.array_equal(x, y), field
    assert a._free.dtype == b._free.dtype
    assert np.array_equal(a._free, b._free)
    for field in ("forces", "velocities"):
        x, y = getattr(a, field), getattr(b, field)
        assert (x is None) == (y is None), field
        if x is not None:
            assert np.array_equal(x, y), field
    assert a.con_energy == b.con_energy


CASES = {
    "plain": dict(),
    "forces": dict(forces=True),
    "velocities": dict(velocities=True),
    "both": dict(forces=True, velocities=True),
    "sorted_ids": dict(shuffle=False),
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_from_conframe_matches_the_per_atom_definition(case):
    text = ref_con_text(make_structure(300, **CASES[case]))
    frame = readcon.read_con_string(text)[0]
    assert_same_structure(Structure.from_conframe(frame), ref_from_conframe(frame))


def test_from_conframe_reads_back_what_to_conframe_wrote():
    s = make_structure(300, forces=True, velocities=True)
    frame = readcon.read_con_string(readcon.write_con_string([s.to_conframe()]))[0]
    back = Structure.from_conframe(frame)
    want = reorder(s)
    assert back.names == want.names
    assert np.array_equal(back.atom_ids, want.atom_ids)
    assert np.allclose(back.r, want.r, rtol=0.0, atol=1e-12)
    assert np.array_equal(back._free, want._free)
    assert np.array_equal(back.mass, want.mass)
    assert np.allclose(back.forces, s.forces[np.argsort(s.atom_ids)], atol=1e-12)
    assert np.allclose(
        back.velocities, s.velocities[np.argsort(s.atom_ids)], atol=1e-12
    )


@pytest.mark.parametrize("case", sorted(CASES))
def test_to_conframe_text_matches_the_per_atom_definition(case):
    s = make_structure(300, **CASES[case])
    assert readcon.write_con_string([s.to_conframe()]) == ref_con_text(s)


def test_to_conframe_collapses_mass_of_a_species_like_the_definition():
    s = make_structure(40)
    s.mass[4] = 70.0
    assert readcon.write_con_string([s.to_conframe()]) == ref_con_text(s)


def test_to_conframe_of_one_species_and_one_mass_matches_the_definition():
    s = make_structure(40)
    s.names = ["Cu"] * 40
    s.mass = np.full(40, 63.546)
    assert readcon.write_con_string([s.to_conframe()]) == ref_con_text(s)
    s.atom_ids = np.arange(40, 0, -1, dtype=np.uint64)
    assert readcon.write_con_string([s.to_conframe()]) == ref_con_text(s)
    s.atom_ids = np.ones(40, dtype=np.uint64)
    assert readcon.write_con_string([s.to_conframe()]) == ref_con_text(s)


def test_to_conframe_without_atom_ids_numbers_them_from_one():
    s = make_structure(40, shuffle=False)
    s.atom_ids = np.arange(5, dtype=np.uint64)
    assert readcon.write_con_string([s.to_conframe()]) == ref_con_text(s)


def _from_conframe_events(n):
    frame = readcon.read_con_string(
        ref_con_text(make_structure(n, forces=True, velocities=True))
    )[0]
    with count_events() as counts:
        Structure.from_conframe(frame)
    return counts


def test_from_conframe_interpreter_work_does_not_grow_with_atoms():
    small, large = _from_conframe_events(64), _from_conframe_events(1024)
    assert total(large) <= total(small), {
        k: (small[k], large[k]) for k in large if large[k] != small[k]
    }


def _to_conframe_events(n):
    s = make_structure(n, forces=True, velocities=True)
    with count_events() as counts:
        s.to_conframe()
    return counts


def test_to_conframe_makes_one_call_per_atom_and_nothing_else_that_grows():
    small, large = _to_conframe_events(64), _to_conframe_events(1024)
    grew = {k: (small[k], large[k]) for k in large if large[k] != small[k]}
    assert set(grew) <= {"c:readcon.Atom", "c:Atom"}, grew
