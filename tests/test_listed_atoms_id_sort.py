"""displace_atom_list is CON file-order; Structure rows are atom_id order."""

from __future__ import annotations

from types import SimpleNamespace


import pytest

from eon.displace import DisplaceError, ListedAtoms
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


def test_listed_atoms_minus_one_means_all_free():
    import readcon

    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 0, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [True, True, True], 1, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    cfg = SimpleNamespace(disp_listed_atoms=[-1], random_mode=False)
    la = ListedAtoms(p, config=cfg)
    free = [i for i, ok in enumerate(p.atom_is_free()) if ok]
    assert la.listed_atoms == free
    assert len(la.listed_atoms) >= 1


def test_listed_atoms_scalar_minus_one_means_all_free():
    import readcon

    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 0, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [True, True, True], 1, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    cfg = SimpleNamespace(disp_listed_atoms=-1, random_mode=False)
    la = ListedAtoms(p, config=cfg)
    assert la.listed_atoms == [i for i, ok in enumerate(p.atom_is_free()) if ok]


def test_listed_atoms_mixed_file_rows_all_remap():
    import readcon

    # File: free id=2, frozen id=0, free id=3, frozen id=1.
    # After sort: frozen 0, frozen 1, free 2, free 3.
    # listed [0, 2] is file-order free atoms. A Structure-row first
    # pass would keep only row 2 and drop file row 0.
    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 2, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [True, True, True], 0, 63.5),
            readcon.Atom("Cu", 2.0, 0.0, 0.0, [False, False, False], 3, 63.5),
            readcon.Atom("Cu", 3.0, 0.0, 0.0, [True, True, True], 1, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    cfg = SimpleNamespace(disp_listed_atoms=[0, 2], random_mode=False)
    la = ListedAtoms(p, config=cfg)
    assert la.listed_atoms == [2, 3]


def test_listed_atoms_script_indices_are_structure_rows():
    import readcon

    # File: free id=100, frozen id=0. After sort, Structure row 0 is
    # frozen and row 1 is free. A kmc script that analyzed
    # savecon(Structure) prints [1] for the free atom. Remapping that
    # as original file-order would send it to the frozen atom.
    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 100, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [True, True, True], 0, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    assert list(p.atom_ids) == [0, 100]
    cfg = SimpleNamespace(
        disp_listed_atoms=[1],
        random_mode=False,
        disp_listed_from_script=True,
    )
    la = ListedAtoms(p, config=cfg)
    assert la.listed_atoms == [1]


def test_listed_atoms_script_frozen_structure_row_raises():
    import readcon

    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 100, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [True, True, True], 0, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    cfg = SimpleNamespace(
        disp_listed_atoms=[0],
        random_mode=False,
        disp_listed_from_script=True,
    )
    with pytest.raises(DisplaceError, match="all frozen"):
        ListedAtoms(p, config=cfg)


def test_listed_atoms_frozen_file_row_does_not_take_other_free():
    import readcon

    # File: free id=1, frozen id=0. After sort, Structure row 1 is free.
    # listed=[1] is the frozen file row. Remap must not keep row 1.
    frame = readcon.ConFrame(
        cell=(10.0, 10.0, 10.0),
        angles=(90.0, 90.0, 90.0),
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 1, 63.5),
            readcon.Atom("Cu", 1.0, 0.0, 0.0, [True, True, True], 0, 63.5),
        ],
    )
    p = Structure.from_conframe(frame)
    cfg = SimpleNamespace(disp_listed_atoms=[1], random_mode=False)
    with pytest.raises(DisplaceError, match="all frozen"):
        ListedAtoms(p, config=cfg)
