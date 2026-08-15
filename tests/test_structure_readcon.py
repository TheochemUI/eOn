"""Structure ↔ readcon.ConFrame is the canonical I/O path (no ASE)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import readcon

from eon.structure import Structure, Atoms, structure_order
from eon import fileio as io

DATA = Path(__file__).resolve().parent / "data"


def interleaved_cu_h() -> Structure:
    """Cu H Cu H, the layout a .con file cannot hold as written.

    Every per-atom field is a different function of the index, so a
    permutation anywhere shows up as a mismatch rather than cancelling out:
    x is the index, the free flags do not follow the species, and the ids
    number the atoms in this order.
    """
    s = Structure(4)
    s.box = np.diag([12.0, 13.0, 14.0])
    s.names = ["Cu", "H", "Cu", "H"]
    s.r = np.array(
        [[0.0, 0.5, 0.5], [1.0, 1.5, 1.5], [2.0, 2.5, 2.5], [3.0, 3.5, 3.5]]
    )
    s.free = np.array([1.0, 1.0, 0.0, 1.0])
    s.mass = np.array([63.546, 1.008, 63.546, 1.008])
    return s


def test_con_angles_are_alpha_beta_gamma():
    """CON header line 4 is alpha beta gamma, not gamma beta alpha."""
    frame = readcon.ConFrame(
        cell=[10.0, 11.0, 12.0],
        angles=[80.0, 85.0, 95.0],
        atoms=[readcon.Atom("Cu", 0.0, 0.0, 0.0, [False, False, False], 1, 63.5)],
        prebox_header=["t", ""],
    )
    s = Structure.from_conframe(frame)
    a, b, c = s.box

    def ang(u, v):
        return float(
            np.degrees(
                np.arccos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            )
        )

    assert ang(b, c) == pytest.approx(80.0, abs=1e-6)
    assert ang(a, c) == pytest.approx(85.0, abs=1e-6)
    assert ang(a, b) == pytest.approx(95.0, abs=1e-6)
    back = s.to_conframe()
    np.testing.assert_allclose(back.angles, [80.0, 85.0, 95.0], atol=1e-6)


def test_atoms_alias_is_structure():
    assert Atoms is Structure


def test_per_axis_con_mask_survives_roundtrip():
    """CON column 4 is per-axis; x-only fixed must not become fully frozen."""
    frame = readcon.ConFrame(
        cell=[10.0, 10.0, 10.0],
        angles=[90.0, 90.0, 90.0],
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [True, False, False], 1, 63.5),
            readcon.Atom("Cu", 1.5, 0.0, 0.0, [False, True, False], 2, 63.5),
        ],
        prebox_header=["t", ""],
    )
    s = Structure.from_conframe(frame)
    np.testing.assert_array_equal(s.free[0], [0.0, 1.0, 1.0])
    np.testing.assert_array_equal(s.free[1], [1.0, 0.0, 1.0])
    back = s.to_conframe()
    assert list(back.atoms[0].fixed) == [True, False, False]
    assert list(back.atoms[1].fixed) == [False, True, False]


def test_from_to_conframe_roundtrip_free_fixed():
    frame = readcon.ConFrame(
        cell=[10.0, 10.0, 10.0],
        angles=[90.0, 90.0, 90.0],
        atoms=[
            readcon.Atom("Cu", 0.0, 0.0, 0.0, [True, True, True], 1, 63.5),
            readcon.Atom("Cu", 1.5, 0.0, 0.0, [False, False, False], 2, 63.5),
        ],
        prebox_header=["t", ""],
    )
    s = Structure.from_conframe(frame)
    assert len(s) == 2
    np.testing.assert_array_equal(s.free[0], [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(s.free[1], [1.0, 1.0, 1.0])
    assert s.names == ["Cu", "Cu"]
    back = s.to_conframe()
    assert list(back.atoms[0].fixed) == [True, True, True]
    assert list(back.atoms[1].fixed) == [False, False, False]
    np.testing.assert_allclose(
        [[a.x, a.y, a.z] for a in back.atoms], s.r, atol=1e-12
    )


def test_loadcon_returns_structure_via_readcon():
    path = DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con"
    s = io.loadcon(str(path))
    assert isinstance(s, Structure)
    ref = readcon.read_con(str(path))[0]
    assert len(s) == len(ref)
    assert s.names[0] == ref.atoms[0].symbol
    np.testing.assert_allclose(s.r[0], [ref.atoms[0].x, ref.atoms[0].y, ref.atoms[0].z])


def test_fileio_does_not_import_ase():
    import eon.fileio as fio
    import inspect

    src = inspect.getsource(fio)
    assert "import ase" not in src
    assert "from ase" not in src
    assert "import readcon" in src


def test_new_structure_numbers_its_atoms_from_one():
    assert list(Structure(3).atom_ids) == [1, 2, 3]
    assert Structure(3).atom_ids.dtype == np.uint64
    assert list(Structure(0).atom_ids) == []


def test_copy_carries_the_atom_ids():
    s = interleaved_cu_h()
    s.atom_ids = np.array([11, 22, 33, 44], dtype=np.uint64)
    assert list(s.copy().atom_ids) == [11, 22, 33, 44]


def test_append_takes_one_past_the_largest_id():
    s = Structure(0)
    s.append([0.0, 0.0, 0.0], 1.0, "Cu", 63.546)
    s.append([1.0, 0.0, 0.0], 1.0, "H", 1.008)
    assert list(s.atom_ids) == [1, 2]
    s.atom_ids = np.array([5, 9], dtype=np.uint64)
    s.append([2.0, 0.0, 0.0], 1.0, "H", 1.008)
    assert list(s.atom_ids) == [5, 9, 10]


def test_append_honours_an_explicit_id():
    s = Structure(0)
    s.append([0.0, 0.0, 0.0], 1.0, "Cu", 63.546, atom_id=41)
    assert list(s.atom_ids) == [41]


def test_structure_order_leaves_ascending_ids_alone():
    ids = np.array([0, 1, 2, 3], dtype=np.uint64)
    np.testing.assert_array_equal(structure_order(ids), [0, 1, 2, 3])


def test_structure_order_inverts_a_grouping():
    # Cu H Cu H written out as Cu Cu H H carries these ids.
    ids = np.array([1, 3, 2, 4], dtype=np.uint64)
    np.testing.assert_array_equal(structure_order(ids), [0, 2, 1, 3])


def test_structure_order_keeps_file_order_when_ids_repeat():
    # A writer that numbers each species from 1 leaves no permutation.
    ids = np.array([1, 2, 1, 2], dtype=np.uint64)
    np.testing.assert_array_equal(structure_order(ids), [0, 1, 2, 3])
    zeros = np.zeros(4, dtype=np.uint64)
    np.testing.assert_array_equal(structure_order(zeros), [0, 1, 2, 3])


def test_structure_order_handles_the_empty_and_single_cases():
    np.testing.assert_array_equal(structure_order(np.zeros(0, np.uint64)), [])
    np.testing.assert_array_equal(structure_order(np.array([7], np.uint64)), [0])


def test_written_con_groups_the_species(tmp_path: Path):
    """The hazard itself: the file cannot hold the in-memory order."""
    s = interleaved_cu_h()
    out = tmp_path / "interleaved.con"
    io.savecon(str(out), s)

    frame = readcon.read_con(str(out))[0]
    assert [a.symbol for a in frame.atoms] == ["Cu", "Cu", "H", "H"]
    # Column 5 records where each atom sat before the grouping.
    assert [a.atom_id for a in frame.atoms] == [1, 3, 2, 4]


def test_interleaved_species_survive_a_con_write_and_read(tmp_path: Path):
    s = interleaved_cu_h()
    out = tmp_path / "interleaved.con"
    io.savecon(str(out), s)
    back = io.loadcon(str(out))

    assert len(back) == len(s)
    assert list(back.names) == ["Cu", "H", "Cu", "H"]
    assert list(back.atom_ids) == [1, 2, 3, 4]
    np.testing.assert_allclose(back.r, s.r, rtol=0, atol=1e-9)
    np.testing.assert_array_equal(back.free, s.free)
    np.testing.assert_allclose(back.mass, s.mass, rtol=0, atol=1e-9)


def test_an_index_into_r_addresses_the_same_atom_after_a_reload(tmp_path: Path):
    """mode.dat is positional; the con it accompanies must agree with it."""
    s = interleaved_cu_h()
    mode = np.array([[float(i), 0.0, 0.0] for i in range(len(s))])

    out = tmp_path / "pos.con"
    modefile = tmp_path / "mode.dat"
    io.savecon(str(out), s)
    io.save_mode(str(modefile), mode)

    back = io.loadcon(str(out))
    back_mode = io.load_mode(str(modefile))

    for i in range(len(back)):
        # Row i of mode.dat belongs to the atom that x == i identifies.
        assert back.r[i][0] == pytest.approx(back_mode[i][0], abs=1e-9)
        assert back.names[i] == s.names[i]


def test_save_mode_zeros_fixed_axes(tmp_path: Path):
    mode = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    free = np.array([[0.0, 1.0, 1.0], [1.0, 1.0, 1.0]])
    out = tmp_path / "mode.dat"
    io.save_mode(str(out), mode, free)
    back = io.load_mode(str(out))
    np.testing.assert_allclose(back, [[0.0, 2.0, 3.0], [4.0, 5.0, 6.0]])


def test_savecon_preserves_the_ids_a_file_arrived_with(tmp_path: Path):
    path = DATA / "server" / "Pt_Heptamer_oneLayer" / "pos.con"
    original_ids = [a.atom_id for a in readcon.read_con(str(path))[0].atoms]

    s = io.loadcon(str(path))
    out = tmp_path / "again.con"
    io.savecon(str(out), s)

    assert [a.atom_id for a in readcon.read_con(str(out))[0].atoms] == original_ids


def test_structure_like_without_ids_writes_sequential_ids():
    """fileio accepts a duck-typed configuration; it numbers one from 1."""

    class Bare:
        def __init__(self):
            self.box = np.diag([10.0, 10.0, 10.0])
            self.names = ["Cu", "H"]
            self.r = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            self.free = np.array([1.0, 1.0])
            self.mass = np.array([63.546, 1.008])

        def __len__(self):
            return 2

    frame = io._atoms_to_frame(Bare())
    assert [a.atom_id for a in frame.atoms] == [1, 2]
    assert [a.symbol for a in frame.atoms] == ["Cu", "H"]


def test_saveposcar_body_matches_its_own_declared_counts(tmp_path: Path):
    """The header declares per-type counts; the coordinate block must agree.

    Parsed without assuming where the counts line sits, so this stays true
    whichever POSCAR convention saveposcar settles on.
    """
    s = interleaved_cu_h()
    out = tmp_path / "POSCAR"
    io.saveposcar(str(out), s)

    lines = out.read_text().splitlines()
    marker = next(
        i for i, line in enumerate(lines) if line.strip().lower().startswith("selective")
    )
    types = lines[0].split()
    counts = [int(c) for c in lines[marker - 1].split()]
    assert types == ["Cu", "H"]
    assert counts == [2, 2]

    body = lines[marker + 2 : marker + 2 + sum(counts)]
    assert len(body) == len(s)
    read_names = [t for t, c in zip(types, counts) for _ in range(c)]
    read_x = [float(line.split()[0]) for line in body]

    # x doubles as each atom's identity in the fixture.
    for name, x in zip(read_names, read_x):
        assert s.names[int(round(x))] == name


def test_poscar_selective_dynamics_is_per_axis(tmp_path: Path):
    s = Structure(1)
    s.box = np.diag([10.0, 10.0, 10.0])
    s.names = ["Cu"]
    s.r = np.array([[1.0, 2.0, 3.0]])
    s.free = np.array([[0.0, 1.0, 1.0]])
    s.mass = np.array([63.546])
    out = tmp_path / "POSCAR"
    io.saveposcar(str(out), s)
    flags = out.read_text().strip().splitlines()[-1].split()[3:6]
    assert flags == ["F", "T", "T"]
