"""Pre-#368 scalar atoms helpers, extracted from commit ee3012498 (parent of #368).

These are the production implementations that shipped before vectorization.
Golden tests compare live eon.atoms outputs to these functions on the same inputs.
Do not edit algorithmically without re-extracting from git history.
"""
from __future__ import annotations

import numpy

from eon.atoms import elements, pbc, per_atom_norm


def free_r_scalar(atoms) -> numpy.ndarray:
    """Historical Atoms.free_r (ee3012498).

    ``int(...)`` on the free count is required under current NumPy: ``sum`` on a
    float free array yields ``numpy.float64``, which is not a valid shape dim.
    Semantics match the pre-#368 loop (truthy free flags).
    """
    free = numpy.asarray(atoms.free)
    if free.ndim == 2:
        atom_free = free.any(axis=1)
    else:
        atom_free = numpy.asarray(free, dtype=bool)
    nfree = int(sum(atom_free))
    temp = numpy.zeros((nfree, 3))
    index = 0
    for i in range(len(atoms.r)):
        if atom_free[i]:
            temp[index] = atoms.r[i]
            index += 1
    return temp


def get_process_atoms_scalar(r, p, epsilon_r=0.2, nshells=1):
    """Historical get_process_atoms (ee3012498)."""
    mobileAtoms = []
    r2p = per_atom_norm(p.r - r.r, r.box)
    for i in range(len(r)):
        if r2p[i] > epsilon_r:
            mobileAtoms.append(i)
    if len(mobileAtoms) == 0:
        mobileAtoms.append(list(r2p).index(max(r2p)))
    neighborAtoms = []
    for atom in mobileAtoms:
        r1 = elements[r.names[atom]]["radius"]
        for i in range(len(r)):
            if i in mobileAtoms or i in neighborAtoms:
                continue
            r2 = elements[r.names[i]]["radius"]
            maxDist = (r1 + r2) * (1.0 + 0.2) * nshells
            if numpy.linalg.norm(pbc(p.r[atom] - p.r[i], p.box)) < maxDist:
                neighborAtoms.append(i)
    return mobileAtoms + neighborAtoms


def brute_neighbor_list_scalar(p, cutoff):
    """Historical brute_neighbor_list (ee3012498)."""
    nl = []
    ibox = numpy.linalg.inv(p.box)
    for a in range(len(p)):
        nl.append([])
        for b in range(len(p)):
            if b != a:
                dist = numpy.linalg.norm(pbc(p.r[a] - p.r[b], p.box, ibox))
                if dist < cutoff:
                    nl[a].append(b)
    return nl
