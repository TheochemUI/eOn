"""Batched geometry helpers against their one-pair-at-a-time definitions.

``pbc``, ``neighbor_list``, ``neighbor_list_vectors`` and
``coordination_numbers`` run on every displacement and every structure
comparison, on cells of a few thousand atoms. The outputs are held to the
per-pair definitions in ``driver_oracles`` bit for bit, and the
interpreter work per call is held independent of the number of atoms.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

minimage = pytest.importorskip("minimage")
linkcell = pytest.importorskip("linkcell")

from driver_oracles import (  # noqa: E402
    count_events,
    ref_neighbor_list_vectors,
    ref_neighbor_lists,
    ref_pair_lists,
    ref_pbc_rows,
)

from eon.geometry import neighbors  # noqa: E402
from eon.geometry.cell import length_angle_to_box  # noqa: E402
from eon.geometry.pbc import pbc  # noqa: E402
from eon.structure import Structure  # noqa: E402

A0 = 3.615


def fcc(cells, jitter=0.05, seed=5, vacancy=True):
    pts = []
    for i in range(cells):
        for j in range(cells):
            for k in range(cells):
                for b in ((0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)):
                    pts.append([(i + b[0]) * A0, (j + b[1]) * A0, (k + b[2]) * A0])
    r = np.array(pts)
    if vacancy:
        r = np.delete(r, len(r) // 2, axis=0)
    rng = np.random.RandomState(seed)
    r = r + rng.normal(scale=jitter, size=r.shape)
    s = Structure(len(r))
    s.r = r % (cells * A0)
    s.box = np.eye(3) * cells * A0
    s.names = ["Cu"] * len(r)
    s.mass = np.full(len(r), 63.546)
    return s


BOXES = {
    "cubic": np.eye(3) * 32.535,
    "orthorhombic": np.diag([32.535, 10.845, 14.46]),
    "edges_differ_in_the_last_bit": np.diag([10.845, 11.5, 9.7]),
    "from_con_angles": length_angle_to_box(
        np.array([A0 * 9, A0 * 9, A0 * 9]), np.array([90.0, 90.0, 90.0])
    ),
    "from_con_angles_unequal": length_angle_to_box(
        np.array([7.3, 8.1, 6.6]), np.array([90.0, 90.0, 90.0])
    ),
    "triclinic": np.array([[10.0, 0.0, 0.0], [2.5, 9.0, 0.0], [1.0, 1.5, 8.0]]),
}


def same_bits(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return (
        a.shape == b.shape
        and a.dtype == b.dtype
        and np.array_equal(a, b)
        and np.array_equal(np.signbit(a), np.signbit(b))
    )


def rows_for(box, n=4000, reach=1.0, seed=11):
    """Random rows, rows on and beside every half-box tie, and zero rows."""
    rng = np.random.RandomState(seed)
    extent = np.abs(box).sum(axis=0)
    rows = rng.uniform(-reach * extent, reach * extent, size=(n, 3))
    half = np.diag(box) / 2.0
    ties = []
    for axis in range(3):
        for sign in (1.0, -1.0):
            for nudge in (0.0, 1e-13, -1e-13, 1e-9, -1e-9):
                row = rng.uniform(-0.4, 0.4, size=3) * np.diag(box)
                row[axis] = sign * half[axis] + nudge
                ties.append(row)
            edge = half[axis]
            for steps in range(-3, 4):
                row = rng.uniform(-0.4, 0.4, size=3) * np.diag(box)
                for _ in range(abs(steps)):
                    edge = np.nextafter(edge, np.inf if steps > 0 else 0.0)
                row[axis] = sign * edge
                ties.append(row)
                edge = half[axis]
    zeros = [[0.0, 0.0, 0.0], [-0.0, 0.0, -0.0]]
    return np.vstack([rows, ties, zeros])


@pytest.mark.parametrize("name", sorted(BOXES))
@pytest.mark.parametrize("reach", [1.0, 2.5])
def test_pbc_rows_match_minimage_displacement(name, reach):
    box = BOXES[name]
    rows = rows_for(box, reach=reach)
    assert same_bits(pbc(rows, box), ref_pbc_rows(rows, box))


def test_pbc_single_row_and_empty_input():
    box = BOXES["cubic"]
    row = np.array([20.0, -3.0, 17.0])
    assert same_bits(pbc(row, box), ref_pbc_rows(row, box)[0])
    assert pbc(np.zeros((0, 3)), box).shape == (0, 3)


def test_pbc_with_a_free_axis_leaves_it_alone():
    box = BOXES["cubic"]
    rows = np.array([[20.0, 20.0, 20.0], [-17.0, 3.0, 18.0]])
    out = pbc(rows, box, periodic=[True, True, False])
    assert np.array_equal(out[:, 2], rows[:, 2])
    assert not np.array_equal(out[:, :2], rows[:, :2])


def _pbc_events(box, n):
    rows = rows_for(box, n=n)
    with count_events() as counts:
        pbc(rows, box)
    return counts


@pytest.mark.parametrize(
    "name",
    [
        "cubic",
        "orthorhombic",
        "edges_differ_in_the_last_bit",
        "from_con_angles",
        "from_con_angles_unequal",
    ],
)
def test_pbc_on_an_orthogonal_cell_makes_no_per_row_calls(name):
    box = BOXES[name]
    small, large = _pbc_events(box, 64), _pbc_events(box, 4096)
    grew = {k: (small[k], large[k]) for k in large if large[k] > small[k]}
    assert not grew, grew


def test_pair_lists_match_the_set_definition_including_images():
    for cells, cutoff in ((4, 3.0), (2, 4.2), (1, 5.0)):
        s = fcc(cells)
        i, j, _shift, _d2 = linkcell.pairs_within(
            np.ascontiguousarray(s.r), np.ascontiguousarray(s.box), cutoff
        )
        i, j = np.asarray(np.from_dlpack(i)), np.asarray(np.from_dlpack(j))
        got = neighbors._pair_lists(len(s), i, j)
        assert got == ref_pair_lists(len(s), i, j)
        assert all(type(v) is int for row in got for v in row)


def test_pair_lists_of_nothing():
    empty = np.zeros(0, dtype=np.int64)
    assert neighbors._pair_lists(3, empty, empty) == [[], [], []]


@pytest.mark.parametrize("cells,cutoff", [(4, 3.0), (5, 5.0)])
def test_neighbor_list_matches_linkcell_definition(cells, cutoff):
    s = fcc(cells)
    assert neighbors.neighbor_list(s, cutoff) == ref_neighbor_lists(s, cutoff)


@pytest.mark.parametrize("cells,cutoff", [(4, 3.0), (5, 5.0)])
def test_neighbor_list_vectors_match_the_per_pair_definition(cells, cutoff):
    s = fcc(cells)
    got = neighbors.neighbor_list_vectors(s, cutoff)
    want = ref_neighbor_list_vectors(s, cutoff)
    assert len(got) == len(want) == len(s)
    for center, (g, w) in enumerate(zip(got, want)):
        assert len(g) == len(w), center
        for a, b in zip(g, w):
            assert same_bits(a, b), center


def test_neighbor_list_vectors_of_an_isolated_atom():
    s = Structure(1)
    s.box = np.eye(3) * 100.0
    s.r = np.array([[1.0, 2.0, 3.0]])
    s.names = ["Cu"]
    s.mass = np.array([63.546])
    assert neighbors.neighbor_list_vectors(s, 3.0) == [[]]


def test_coordination_numbers_are_the_neighbor_counts_as_plain_ints():
    s = fcc(4)
    got = neighbors.coordination_numbers(s, 3.0)
    assert got == [len(row) for row in ref_neighbor_lists(s, 3.0)]
    assert all(type(c) is int for c in got)
    assert min(got) < 12 <= max(got)


def test_open_axis_neighbors_skip_the_wrap_on_that_axis():
    s = fcc(3, jitter=0.1)
    s.periodic = np.array([True, True, False])
    cutoff = 3.2
    got = neighbors.neighbor_list(s, cutoff)
    box = np.diag(s.box)
    want = []
    for a in range(len(s)):
        d = s.r - s.r[a]
        d[:, :2] -= box[:2] * np.round(d[:, :2] / box[:2])
        dist = np.sqrt((d**2).sum(axis=1))
        want.append([int(b) for b in np.flatnonzero(dist < cutoff) if b != a])
    assert got == want


def _list_events(fn, cells):
    s = fcc(cells, vacancy=False)
    with count_events() as counts:
        fn(s, 3.0)
    return counts


@pytest.mark.parametrize(
    "fn",
    [
        neighbors.neighbor_list,
        neighbors.neighbor_list_vectors,
        neighbors.coordination_numbers,
    ],
    ids=lambda f: f.__name__,
)
def test_interpreter_work_does_not_grow_with_atoms(fn):
    small, large = _list_events(fn, 3), _list_events(fn, 6)
    grew = {k: (small[k], large[k]) for k in large if large[k] > small[k]}
    assert not grew, grew
