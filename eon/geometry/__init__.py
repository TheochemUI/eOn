"""Geometry kernels for the eOn Python server.

- :mod:`eon.geometry.cell` — cell length/angle transforms
- :mod:`eon.geometry.pbc` — minimum-image PBC and per-atom norms
- :mod:`eon.geometry.neighbors` — neighbor lists via **vesin**
  (optional :func:`neighbor_list_linkcell` for linkcell/minimage checks)
- :mod:`eon.geometry.process` — process-atom selection for recycling
"""

from eon.geometry.cell import box_to_length_angle, length_angle_to_box
from eon.geometry.neighbors import (
    brute_neighbor_list,
    coordination_numbers,
    least_coordinated,
    neighbor_list,
    neighbor_list_linkcell,
    neighbor_list_vectors,
)
from eon.geometry.pbc import pbc, per_atom_norm, per_atom_norm_gen
from eon.geometry.process import get_process_atoms

__all__ = [
    "box_to_length_angle",
    "length_angle_to_box",
    "pbc",
    "per_atom_norm",
    "per_atom_norm_gen",
    "neighbor_list",
    "neighbor_list_linkcell",
    "neighbor_list_vectors",
    "brute_neighbor_list",
    "coordination_numbers",
    "least_coordinated",
    "get_process_atoms",
]
