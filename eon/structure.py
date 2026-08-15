"""Structure storage for the eOn Python server.

Canonical on-disk form is ``readcon.ConFrame`` (``.con`` via the readcon package).
In-memory working form is :class:`Structure`: contiguous numpy arrays for geometry
algorithms (PBC, neighbor lists, displacements). This replaces the historical
aselite-era mini-Atoms type without requiring ASE.

``Atoms`` is an alias of :class:`Structure` for call-site compatibility.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Union

import numpy as np
import readcon

from eon.geometry.cell import box_to_length_angle, length_angle_to_box


def coerce_free(value, n: int) -> np.ndarray:
    """Normalize a free-flag array to shape ``(n, 3)``.

    A length-``n`` vector is the whole-atom form and is broadcast onto every
    Cartesian axis. ``1.0`` is free, ``0.0`` is fixed.
    """
    arr = np.asarray(value, dtype=float)
    if n == 0:
        return np.zeros((0, 3), dtype=float)
    if arr.shape == (n, 3):
        return np.array(arr, dtype=float, copy=True)
    if arr.shape == (n,):
        return np.repeat(arr.reshape(n, 1), 3, axis=1)
    if arr.shape == (3,) and n == 1:
        return np.array(arr, dtype=float, copy=True).reshape(1, 3)
    raise ValueError(f"free must have shape ({n},) or ({n}, 3), got {arr.shape}")


def as_atom_free(free) -> np.ndarray:
    """Per-atom free flag: True if any Cartesian axis is free."""
    f = np.asarray(free, dtype=float)
    if f.ndim == 1:
        return f > 0.5
    if f.ndim == 2 and f.shape[-1] == 3:
        return (f > 0.5).any(axis=1)
    raise ValueError(f"free must be shape (N,) or (N, 3), got {f.shape}")


def structure_order(atom_ids: np.ndarray) -> np.ndarray:
    """Permutation taking file order to :class:`Structure` order.

    ``Structure`` orders its atoms by ascending ``atom_id``; a ``.con`` file
    orders them by species block, because CON header lines 7 and 8 are a type
    count and per-type counts. Ascending ids therefore survive the grouping the
    writer applies and are what puts the atoms back where an index into
    :attr:`Structure.r` finds the same atom it found before the write.

    Ids that repeat carry no permutation to invert (a writer that numbers each
    species from 1 produces those), so file order stands.
    """
    ids = np.asarray(atom_ids)
    n = ids.shape[0]
    if n < 2:
        return np.arange(n)
    if bool(np.all(ids[1:] > ids[:-1])):
        return np.arange(n)
    if np.unique(ids).shape[0] != n:
        return np.arange(n)
    return np.argsort(ids, kind="stable")


class Structure:
    """Mutable atomic configuration (numpy-backed).

    Atoms sit in ascending :attr:`atom_ids` order, which is the order a
    from-scratch ``Structure`` builds them in and the order every index-coupled
    sidecar (``mode.dat``, ``masses.dat``, ``direction.dat``, the Hessian,
    ``Prefactor.movedAtoms``) is written in. ``.con`` files store atoms grouped
    by species, so :meth:`from_conframe` undoes the grouping through the ids
    rather than taking file order as given.

    Attributes
    ----------
    r : (N, 3) float
        Cartesian positions.
    free : (N, 3) float
        Per-axis free flags (1.0 free, 0.0 fixed). A length-N assignment
        broadcasts onto all three axes. CON column 4 is this mask.
    box : (3, 3) float
        Cell matrix with lattice vectors as **rows** (ASE/vesin convention).
    names : list[str]
        Chemical symbols, length N.
    mass : (N,) float
        Atomic masses.
    atom_ids : (N,) uint64
        Per-atom identity, CON column 5. ``1..N`` for a Structure built from
        scratch; whatever the file carried for one read off disk.
    """

    __slots__ = ("r", "_free", "box", "names", "mass", "atom_ids")

    def __init__(self, n_atoms: int = 0):
        self.r = np.zeros((n_atoms, 3), dtype=float)
        self._free = np.ones((n_atoms, 3), dtype=float)
        self.box = np.zeros((3, 3), dtype=float)
        self.names: List[str] = [""] * n_atoms
        self.mass = np.zeros(n_atoms, dtype=float)
        self.atom_ids = np.arange(1, n_atoms + 1, dtype=np.uint64)

    @property
    def free(self) -> np.ndarray:
        return self._free

    @free.setter
    def free(self, value) -> None:
        self._free = coerce_free(value, len(self))

    def __len__(self) -> int:
        return int(self.r.shape[0])

    def copy(self) -> "Structure":
        p = Structure(len(self))
        p.r = self.r.copy()
        p.free = self.free.copy()
        p.box = self.box.copy()
        p.names = list(self.names)
        p.mass = self.mass.copy()
        p.atom_ids = self.atom_ids.copy()
        return p

    def ids_or_sequential(self) -> np.ndarray:
        """Atom ids, or ``1..N`` when the id array is out of step with the geometry.

        Callers that grow ``r`` / ``names`` / ``mass`` field by field leave the
        ids behind; a short or long id array carries no identity to preserve,
        so the sequential numbering stands in.
        """
        n = len(self)
        ids = np.asarray(self.atom_ids, dtype=np.uint64).reshape(-1)
        if ids.shape[0] == n:
            return ids
        return np.arange(1, n + 1, dtype=np.uint64)

    def atom_is_free(self) -> np.ndarray:
        """True for atoms that have at least one free Cartesian axis."""
        return as_atom_free(self._free)

    def free_r(self) -> np.ndarray:
        """Positions of atoms that have at least one free axis."""
        return self.r[self.atom_is_free()]

    def free_mask(self) -> np.ndarray:
        return np.asarray(self._free, dtype=float) > 0.5

    def fixed_mask(self) -> np.ndarray:
        return ~self.free_mask()

    def append(self, r, free, name, mass, atom_id: Optional[int] = None) -> None:
        """Add one atom at the end.

        An atom_id of None takes one past the largest in use, keeping the ids
        distinct and ascending so the appended atom stays last across a save
        and load cycle.
        """
        if atom_id is None:
            atom_id = int(self.atom_ids.max()) + 1 if len(self.atom_ids) else 1
        self.r = np.append(self.r, [r], 0)
        free_row = np.asarray(free, dtype=float).reshape(-1)
        if free_row.size == 1:
            free_row = np.repeat(free_row, 3)
        elif free_row.size != 3:
            raise ValueError("free must be a scalar or length-3")
        if self._free.shape[0] == 0:
            self._free = free_row.reshape(1, 3)
        else:
            self._free = np.vstack((self._free, free_row.reshape(1, 3)))
        self.names.append(name)
        self.mass = np.append(self.mass, mass)
        self.atom_ids = np.concatenate(
            (self.atom_ids, np.array([atom_id], dtype=np.uint64))
        )

    # --- readcon ConFrame bridge (canonical I/O type) ---

    @classmethod
    def from_conframe(cls, frame: "readcon.ConFrame") -> "Structure":
        """Build a Structure from a readcon ConFrame (live API).

        The frame holds atoms grouped by species; :func:`structure_order` puts
        them back in ``atom_id`` order.
        """
        n = len(frame)
        p = cls(n)
        boxlengths = np.asarray(list(frame.cell), dtype=float)
        boxangles = np.asarray(list(frame.angles), dtype=float)
        p.box = length_angle_to_box(boxlengths, boxangles)
        # Prefer bulk coords when available
        try:
            coords = np.asarray(frame.coords_array(), dtype=float)
            if coords.shape == (n, 3):
                p.r = coords.copy()
            else:
                raise ValueError("coords shape")
        except Exception:
            for i, atom in enumerate(frame.atoms):
                p.r[i] = [atom.x, atom.y, atom.z]
        try:
            ids = np.asarray(frame.atom_ids_array(), dtype=np.uint64)
            if ids.shape != (n,):
                raise ValueError("atom_ids shape")
        except Exception:
            ids = np.array(
                [atom.atom_id for atom in frame.atoms], dtype=np.uint64
            ).reshape(n)
        p.atom_ids = ids
        for i, atom in enumerate(frame.atoms):
            p.names[i] = atom.symbol
            p.mass[i] = atom.mass if atom.mass is not None else 0.0
            fixed = atom.fixed
            if fixed is None:
                p._free[i] = (1.0, 1.0, 1.0)
            else:
                p._free[i] = tuple(0.0 if flag else 1.0 for flag in fixed)
        order = structure_order(p.atom_ids)
        if not np.array_equal(order, np.arange(n)):
            p.r = p.r[order]
            p.free = p.free[order]
            p.mass = p.mass[order]
            p.atom_ids = p.atom_ids[order]
            p.names = [p.names[j] for j in order]
        return p

    def to_conframe(
        self,
        prebox_header: Optional[Sequence[str]] = None,
        postbox_header: Optional[Sequence[str]] = None,
    ) -> "readcon.ConFrame":
        """Convert to a readcon ConFrame for writing.

        Atoms go out in Structure order carrying their own ids. The writer
        groups them by species; the ids are what :meth:`from_conframe` reads
        the grouping back out of.
        """
        lengths, angles = box_to_length_angle(self.box)
        atom_ids = self.ids_or_sequential()
        atom_list = []
        for i in range(len(self)):
            fixed = [bool(self._free[i, a] < 0.5) for a in range(3)]
            # Column 4 encodes x-only as 1; readcon-core decodes 1 as
            # fully fixed. Refuse the write rather than emit a file
            # whose constraints differ from the ones in memory.
            if fixed[0] and not fixed[1] and not fixed[2]:
                raise ValueError(
                    f"atom {i} is fixed in x only, which CON column 4 "
                    "writes as 1 and readers decode as fully fixed; "
                    "constrain another axis or leave the atom free"
                )
            atom_list.append(
                readcon.Atom(
                    symbol=self.names[i],
                    x=float(self.r[i][0]),
                    y=float(self.r[i][1]),
                    z=float(self.r[i][2]),
                    fixed=fixed,
                    atom_id=int(atom_ids[i]),
                    mass=float(self.mass[i]),
                )
            )
        return readcon.ConFrame(
            cell=list(lengths),
            angles=list(angles),
            atoms=atom_list,
            prebox_header=list(prebox_header)
            if prebox_header is not None
            else ["Generated by eOn", ""],
            postbox_header=list(postbox_header) if postbox_header is not None else None,
        )


# Back-compat name used throughout the server
Atoms = Structure


def structure_to_ase(structure: "Structure", *, pbc: bool = True):
    """Optional ASE export (requires ase + pyeonclient.ase_bridge)."""
    from pyeonclient.ase_bridge import structure_to_ase as _f

    return _f(structure, pbc=pbc)


def ase_to_structure(atoms):
    """Optional ASE import (requires ase + pyeonclient.ase_bridge)."""
    from pyeonclient.ase_bridge import ase_to_structure as _f

    return _f(atoms)


# Methods on Structure for ergonomics
def _structure_to_ase(self, *, pbc: bool = True):
    return structure_to_ase(self, pbc=pbc)


Structure.to_ase = _structure_to_ase  # type: ignore[attr-defined]


def _structure_from_ase(cls, atoms):
    return ase_to_structure(atoms)


Structure.from_ase = classmethod(_structure_from_ase)  # type: ignore[attr-defined, assignment]
