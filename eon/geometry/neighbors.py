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


def neighbor_list(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,  # noqa: ARG001 — API compat; vesin always used
) -> List[List[int]]:
    """Return neighbors within *cutoff* for each atom (PBC, full undirected list).

    Parameters
    ----------
    p
        Structure with ``.r`` (N,3) and ``.box`` (3,3).
    cutoff
        Pair cutoff distance.
    brute
        Ignored; kept so callers using ``config.comp_brute_neighbors`` need no change.
    """
    r, box = _positions_box(p)
    n = r.shape[0]
    if n == 0:
        return []
    if cutoff <= 0:
        return [[] for _ in range(n)]
    calc = VesinNeighborList(cutoff=float(cutoff), full_list=True)
    i, j = calc.compute(r, box, periodic=True, quantities="ij")
    return _pair_lists(n, np.asarray(i), np.asarray(j))


def brute_neighbor_list(p: StructureLike, cutoff: float) -> List[List[int]]:
    """Alias of :func:`neighbor_list` (historical name)."""
    return neighbor_list(p, cutoff, brute=True)


def neighbor_list_vectors(
    p: StructureLike,
    cutoff: float,
    brute: bool = False,
) -> List[List[np.ndarray]]:
    """Neighbor list with minimum-image vectors from center → neighbor."""
    nl = neighbor_list(p, cutoff, brute=brute)
    r, box = _positions_box(p)
    ibox = np.linalg.inv(box)
    out: List[List[np.ndarray]] = []
    for center, neighs in enumerate(nl):
        vecs = []
        for j in neighs:
            vecs.append(pbc(r[j] - r[center], box, ibox))
        out.append(vecs)
    return out


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
    nn, d2 = linkcell.knearest(xyz, cell, kk)
    nn = np.from_dlpack(nn)
    d2 = np.from_dlpack(d2)
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
