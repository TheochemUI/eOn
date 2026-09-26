"""Lists the wire Geometry can carry beyond positions and the cell.

readcon ``ConFrame`` is the on-disk codec. ``Geometry`` in
``schema/eon_job_result.capnp`` is the runtime contract. These lists are
the append-only fields that keep forces, atom ids, and per-axis constraints
on that contract. An empty list means the field is absent.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from eon.structure import Structure


def geometry_wire_lists(
    structure: Structure, forces: Any = None
) -> dict[str, list[Any]]:
    """Return ``forces``, ``atomId``, and ``fixedAxes`` for one Structure."""
    n = len(structure)
    if forces is None:
        force_list: list[float] = []
    else:
        arr = np.asarray(forces, dtype=np.float64).reshape(-1)
        force_list = [] if arr.size == 0 else arr.tolist()
    ids = np.asarray(structure.atom_ids, dtype=np.uint64).reshape(-1)
    atom_id = ids.tolist() if ids.shape[0] == n else []
    free = np.asarray(structure._free, dtype=float)
    if free.size == 0:
        bits: list[int] = []
    else:
        if free.ndim == 1:
            free = np.repeat(free.reshape(-1, 1), 3, axis=1)
        bits = []
        for i in range(n):
            word = 0
            for axis in range(3):
                if float(free[i, axis]) <= 0.5:
                    word |= 1 << axis
            bits.append(word)
    return {"forces": force_list, "atomId": atom_id, "fixedAxes": bits}
