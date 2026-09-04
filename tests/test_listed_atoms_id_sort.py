"""displace_atom_list is CON file-order; Structure rows are atom_id order."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from eon.displace import ListedAtoms
from eon.structure import Structure, file_rows_to_structure_rows


def test_file_rows_to_structure_rows_maps_past_id_sort():
    assert file_rows_to_structure_rows([100, 101, 0, 1], [0, 1]) == [2, 3]
    assert file_rows_to_structure_rows([0, 1, 2, 3], [0, 1]) == [0, 1]


def test_listed_atoms_remaps_file_order_when_raw_rows_are_frozen():
    import readcon

    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 100, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [False, False, False], 101, 63.5),
            readcon.Atom("Cu", 2.0, 0.0, 0.0, [True, True, True], 0, 63.5),
            readcon.Atom("Cu", 3.0, 0.0, 0.0, [True, True, True], 1, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    # File rows 0,1 are free (ids 100,101). After sort, Structure rows
    # 0,1 are ids 0,1 (frozen).
    cfg = SimpleNamespace(disp_listed_atoms=[0, 1], random_mode=False)
    la = ListedAtoms(p, config=cfg)
    assert la.listed_atoms == [2, 3]
