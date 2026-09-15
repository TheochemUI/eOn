"""Neighbor lists via vesin (PBC-aware).

Historically eOn used a Python sweep-and-prune or an O(N²) brute loop.
Both paths now go through :class:`vesin.NeighborList`, which is the supported
geometry kernel for pair finding. The *brute* flag is retained for API
compatibility but no longer selects a different algorithm.
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
    """Resolve vesin ``periodic`` from the explicit argument or ``p.periodic``."""
    if periodic is not None:
        return periodic
    flags = getattr(p, "periodic", None)
    if flags is None:
        flags = getattr(p, "pbc", None)
    if flags is None:
        return True
    return flags


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
        Passed to vesin. A single bool applies to all axes; a length-3
        sequence is per-axis. Default is ``p.periodic`` or all-periodic.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    if n == 0:
        return []
    if cutoff <= 0:
        return [[] for _ in range(n)]
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j = calc.compute(
        r, box, periodic=_periodic_flags(p, periodic), quantities="ij"
    )
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
    """Vesin pair list with cell shifts (ASE/tonari ``ijS``).

    Returns ``(i, j, S)`` with one row per atom-image pair. Unlike
    :func:`neighbor_list`, this does not unique-index or apply a
    minimum-image reduction. Displacement is
    ``r[j] - r[i] + S @ box``.
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
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j, S = calc.compute(
        r, box, periodic=_periodic_flags(p, periodic), quantities="ijS"
    )
    return np.asarray(i), np.asarray(j), np.asarray(S)


def neighbor_list_linkcell(
    p: StructureLike,
    cutoff: float,
    k: int | None = None,
) -> List[List[int]]:
    """Neighbor list via :func:`linkcell.knearest`, filtered to *cutoff*.

    Production :func:`neighbor_list` stays on vesin. This path exists so
    the two kernels can be compared on the same Structure.
    """
    import linkcell

    r, box = _positions_box(p)
    n = r.shape[0]
    if n == 0:
        return []
    if cutoff <= 0 or n == 1:
        return [[] for _ in range(n)]
    kk = n - 1 if k is None else int(k)
    if kk < 1:
        return [[] for _ in range(n)]
    xyz = np.ascontiguousarray(r, dtype=np.float64)
    cell = np.ascontiguousarray(box, dtype=np.float64)
    raw = linkcell.knearest(xyz, cell, kk)
    if isinstance(raw, tuple) and len(raw) == 2:
        nn = np.from_dlpack(raw[0])
        d2 = np.from_dlpack(raw[1])
    else:
        nn = np.from_dlpack(raw)
        import minimage

        mi = minimage.Cell.from_vesin(box.tolist())
        d2 = np.empty(nn.shape, dtype=np.float64)
        for i in range(n):
            for j in range(kk):
                idx = int(nn[i, j])
                if idx < 0:
                    d2[i, j] = np.nan
                else:
                    d2[i, j] = mi.dist2(r[i].tolist(), r[idx].tolist())
    cut2 = float(cutoff) * float(cutoff)
    out: List[List[int]] = []
    for i in range(n):
        neigh = []
        for j in range(kk):
            idx = int(nn[i, j])
            if idx < 0 or idx == i:
                continue
            if float(d2[i, j]) <= cut2:
                neigh.append(idx)
        out.append(sorted(set(neigh)))
    return out


def coordination_numbers(
    p: StructureLike, cutoff: float, brute: bool = False
) -> List[int]:
    return [len(l) for l in neighbor_list(p, cutoff, brute)]


def least_coordinated(
    p: StructureLike, cutoff: float, brute: bool = False
) -> List[int]:
    """Indices of free atoms with the lowest coordination number."""
    cn = coordination_numbers(p, cutoff, brute)
    if not cn:
        return []
    maxcoord = max(cn)
    mincoord = min(cn)
    free = np.asarray(getattr(p, "free", np.ones(len(cn))), dtype=bool)
    while mincoord <= maxcoord:
        least = [i for i in range(len(cn)) if cn[i] <= mincoord and free[i]]
        if least:
            return least
        mincoord += 1
    return []
