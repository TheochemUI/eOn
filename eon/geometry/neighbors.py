"""Neighbor lists.

A periodic cutoff list is :func:`linkcell.pairs_within`: one row per
atom-image, the caller's shift ``S``, and a squared distance strictly
below the cutoff squared. :func:`neighbor_list` unique-indexes those
rows. An open axis stays on :class:`vesin.NeighborList`, which accepts
per-axis periodicity. ``knearest`` is the k-nearest list, not a cutoff
list. The *brute* flag is retained for API compatibility and does not
select another algorithm.
"""

from __future__ import annotations

from typing import List, Sequence, Union

import numpy as np
from vesin import NeighborList as VesinNeighborList

from eon.geometry.pbc import pbc

PeriodicSpec = Union[bool, Sequence[bool], np.ndarray]

# Structure-like: needs .r, .box, __len__, and optionally .free
StructureLike = object


def _positions_box(p: StructureLike):
    r = np.asarray(p.r, dtype=float)
    box = np.asarray(p.box, dtype=float)
    if r.ndim != 2 or r.shape[1] != 3:
        raise ValueError(f"positions must be (N, 3), got {r.shape}")
    if box.shape != (3, 3):
        raise ValueError(f"box must be (3, 3), got {box.shape}")
    return r, box


def _pair_lists(n: int, i: np.ndarray, j: np.ndarray) -> List[List[int]]:
    """Convert vesin pair indices to eOn adjacency lists (plain Python ints).

    Vesin may return the same neighbor index multiple times under different
    periodic shifts when the cutoff reaches multiple images. eOn historically
    stores unique atom indices only (no multi-image multiplicity).
    """
    nl_sets = [set() for _ in range(n)]
    for a, b in zip(i.tolist(), j.tolist()):
        ai, bi = int(a), int(b)
        if ai == bi:
            continue
        nl_sets[ai].add(bi)
    return [sorted(s) for s in nl_sets]


def _periodic_flags(p: StructureLike, periodic: PeriodicSpec | None) -> PeriodicSpec:
    """Resolve ``periodic`` from the explicit argument or ``p.periodic``."""
    if periodic is not None:
        return periodic
    flags = getattr(p, "periodic", None)
    if flags is None:
        flags = getattr(p, "pbc", None)
    if flags is None:
        return True
    return flags


def _all_periodic(flags: PeriodicSpec) -> bool:
    """True when every axis wraps. ``pairs_within`` has no open axis."""
    if isinstance(flags, (bool, np.bool_)):
        return bool(flags)
    arr = np.asarray(flags, dtype=bool).reshape(-1)
    if arr.size == 1:
        return bool(arr[0])
    return arr.size == 3 and bool(arr.all())


def _linkcell_pairs(r: np.ndarray, box: np.ndarray, cutoff: float, *, half: bool = False):
    """``(i, j, S, dist2)`` from :func:`linkcell.pairs_within`, or ``None``."""
    try:
        import linkcell
    except ImportError:
        return None
    pairs_within = getattr(linkcell, "pairs_within", None)
    if pairs_within is None:
        return None
    xyz = np.ascontiguousarray(r, dtype=np.float64)
    cell = np.ascontiguousarray(box, dtype=np.float64)
    i, j, shift, dist2 = pairs_within(xyz, cell, float(cutoff), half=half)
    return (
        np.asarray(np.from_dlpack(i), dtype=np.int64),
        np.asarray(np.from_dlpack(j), dtype=np.int64),
        np.asarray(np.from_dlpack(shift), dtype=np.int32),
        np.asarray(np.from_dlpack(dist2), dtype=np.float64),
    )


def neighbor_list(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,  # noqa: ARG001 — API compat; vesin always used
    periodic: PeriodicSpec | None = None,
) -> List[List[int]]:
    """Return neighbors within *cutoff* for each atom (PBC, full undirected list).

    Parameters
    ----------
    p
        Structure with ``.r`` (N,3) and ``.box`` (3,3). Optional
        ``.periodic`` (bool or length-3 bools) is used when *periodic*
        is omitted.
    cutoff
        Pair cutoff distance.
    brute
        Ignored; kept so callers using ``config.comp_brute_neighbors`` need no change.
    periodic
        A single bool applies to all axes; a length-3 sequence is
        per-axis. Default is ``p.periodic`` or all-periodic. Fully
        periodic cells use :func:`linkcell.pairs_within`. An open axis
        uses vesin.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    if n == 0:
        return []
    if cutoff <= 0:
        return [[] for _ in range(n)]
    flags = _periodic_flags(p, periodic)
    if _all_periodic(flags):
        packed = _linkcell_pairs(r, box, cutoff, half=False)
        if packed is not None:
            return _pair_lists(n, packed[0], packed[1])
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j = calc.compute(r, box, periodic=flags, quantities="ij")
    return _pair_lists(n, np.asarray(i), np.asarray(j))


def brute_neighbor_list(p: StructureLike, cutoff: float) -> List[List[int]]:
    """Alias of :func:`neighbor_list` (historical name)."""
    return neighbor_list(p, cutoff, brute=True)


def neighbor_list_vectors(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,
) -> List[List[np.ndarray]]:
    """Neighbor list with minimum-image vectors from center → neighbor.

    Unique-index adjacency is the historical eOn contract. Vectors are
    one MIC wrap of ``r[j] - r[center]`` per neighbour index, computed
    in one :func:`eon.geometry.pbc` call (minimage ``wrap_many`` when
    installed). Multi-image pair lists belong on
    :func:`neighbor_list_pairs`.
    """
    nl = neighbor_list(p, cutoff, brute=brute)
    r, box = _positions_box(p)
    ibox = np.linalg.inv(box)
    pairs = [(center, j) for center, neighs in enumerate(nl) for j in neighs]
    if not pairs:
        return [[] for _ in nl]
    diffs = np.empty((len(pairs), 3), dtype=float)
    for k, (center, j) in enumerate(pairs):
        diffs[k] = r[j] - r[center]
    wrapped = np.atleast_2d(
        pbc(diffs, box, ibox, periodic=_periodic_flags(p, None))
    )
    out: List[List[np.ndarray]] = [[] for _ in nl]
    for k, (center, _j) in enumerate(pairs):
        out[center].append(np.asarray(wrapped[k], dtype=float))
    return out


def neighbor_list_pairs(
    p: StructureLike,
    cutoff: float,
    periodic: PeriodicSpec | None = None,
):
    """Cutoff pair list with cell shifts (ASE/tonari ``ijS``).

    Returns ``(i, j, S)`` with one row per atom-image pair, ordered by
    ``(i, j, S)``. Unlike :func:`neighbor_list`, this does not
    unique-index or apply a minimum-image reduction. Displacement is
    ``r[j] - r[i] + S @ box``. A fully periodic cell uses
    :func:`linkcell.pairs_within`. An open axis uses vesin. A squared
    distance on a linkcell row is strictly below the cutoff squared.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    empty = (
        np.zeros(0, dtype=np.int64),
        np.zeros(0, dtype=np.int64),
        np.zeros((0, 3), dtype=np.int32),
    )
    if n == 0 or cutoff <= 0:
        return empty
    flags = _periodic_flags(p, periodic)
    if _all_periodic(flags):
        packed = _linkcell_pairs(r, box, cutoff, half=False)
        if packed is not None:
            i, j, shift, _dist2 = packed
            if i.size == 0:
                return empty
            order = np.lexsort((shift[:, 2], shift[:, 1], shift[:, 0], j, i))
            return i[order], j[order], shift[order]
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j, shift = calc.compute(r, box, periodic=flags, quantities="ijS")
    i = np.asarray(i)
    j = np.asarray(j)
    shift = np.asarray(shift)
    if i.size == 0:
        return empty
    order = np.lexsort((shift[:, 2], shift[:, 1], shift[:, 0], j, i))
    return i[order], j[order], shift[order]


def neighbor_list_linkcell(
    p: StructureLike,
    cutoff: float,
) -> List[List[int]]:
    """Unique-index adjacency from :func:`linkcell.pairs_within`.

    The cell is periodic. ``knearest`` is a different query and is not
    called. An empty list is returned when linkcell has no
    ``pairs_within``.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    if n == 0:
        return []
    if cutoff <= 0 or n == 1:
        return [[] for _ in range(n)]
    packed = _linkcell_pairs(r, box, cutoff, half=False)
    if packed is None:
        return [[] for _ in range(n)]
    return _pair_lists(n, packed[0], packed[1])


def coordination_numbers(
    p: StructureLike, cutoff: float, brute: bool = False
) -> List[int]:
    return [len(l) for l in neighbor_list(p, cutoff, brute)]


def least_coordinated(
    p: StructureLike, cutoff: float, brute: bool = False
) -> List[int]:
    """Indices of free atoms with the lowest coordination number."""
    from eon.structure import as_atom_free

    cn = coordination_numbers(p, cutoff, brute)
    if not cn:
        return []
    maxcoord = max(cn)
    mincoord = min(cn)
    free = as_atom_free(getattr(p, "free", np.ones(len(cn))))
    while mincoord <= maxcoord:
        least = [i for i in range(len(cn)) if cn[i] <= mincoord and free[i]]
        if least:
            return least
        mincoord += 1
    return []
