"""Reference implementations and call counters for the server-driver tests.

Each ``ref_*`` function is the plain per-element definition of what the
production helper of the same name computes: one loop iteration per atom or
per pair, no batching. The production helpers are batched; the tests hold
them to these definitions bit for bit.

The counters watch the interpreter, not the helpers, so a batched helper
cannot satisfy them by being renamed.
"""

from __future__ import annotations

import contextlib
import sys
from collections import Counter

import numpy as np


@contextlib.contextmanager
def count_events():
    """Count Python calls and C calls made inside the block.

    Yields a ``Counter`` keyed by ``"py:<qualname>"`` for Python functions
    and ``"c:<module>.<name>"`` for builtin and extension callables.
    """
    counts: Counter = Counter()

    def profiler(frame, event, arg):
        if event == "call":
            code = frame.f_code
            counts["py:" + getattr(code, "co_qualname", code.co_name)] += 1
        elif event == "c_call":
            owner = getattr(arg, "__module__", None) or type(
                getattr(arg, "__self__", None)
            ).__name__
            counts["c:%s.%s" % (owner, arg.__name__)] += 1

    previous = sys.getprofile()
    sys.setprofile(profiler)
    try:
        yield counts
    finally:
        sys.setprofile(previous)


def total(counts: Counter, *, prefix: str = "") -> int:
    return sum(v for k, v in counts.items() if k.startswith(prefix))


def ref_from_conframe(frame):
    """One atom at a time: the definition ``Structure.from_conframe`` batches."""
    from eon.structure import (
        Structure,
        length_angle_to_box,
        structure_order,
    )

    n = len(frame)
    p = Structure(n)
    p.box = length_angle_to_box(
        np.asarray(list(frame.cell), dtype=float),
        np.asarray(list(frame.angles), dtype=float),
    )
    p.r = np.array(np.asarray(frame.xyz, dtype=float), dtype=float, copy=True)
    p.atom_ids = np.asarray(frame.atom_ids_array(), dtype=np.uint64)
    for i, atom in enumerate(frame.atoms):
        p.names[i] = atom.symbol
        p.mass[i] = atom.mass if atom.mass is not None else 0.0
        fixed = atom.fixed
        if fixed is None:
            p._free[i] = (1.0, 1.0, 1.0)
        else:
            p._free[i] = tuple(0.0 if flag else 1.0 for flag in fixed)
    order = structure_order(p.atom_ids)
    file_to_struct = np.empty(n, dtype=np.int64)
    file_to_struct[order] = np.arange(n, dtype=np.int64)
    p.file_to_struct = file_to_struct
    identity = np.array_equal(order, np.arange(n))
    if not identity:
        p.r = p.r[order]
        p.free = p.free[order]
        p.mass = p.mass[order]
        p.atom_ids = p.atom_ids[order]
        p.names = [p.names[j] for j in order]

    def section(attrs):
        rows = []
        present = False
        for atom in frame.atoms:
            vals = [getattr(atom, name) for name in attrs]
            if any(val is not None for val in vals):
                present = True
            rows.append([0.0 if val is None else float(val) for val in vals])
        if not present:
            return None
        arr = np.array(rows, dtype=float)
        return arr if identity else arr[order]

    p.forces = section(("fx", "fy", "fz"))
    p.velocities = section(("vx", "vy", "vz"))
    meta = getattr(frame, "metadata", None)
    if meta is not None:
        try:
            energy = meta["energy"]
        except Exception:
            energy = None
        if energy is not None:
            p.con_energy = float(energy)
    return p


def ref_con_text(structure) -> str:
    """The con text of *structure*, atom by atom through ``readcon.Atom``."""
    import readcon

    from eon.structure import box_to_length_angle

    lengths, angles = box_to_length_angle(structure.box)
    ids = structure.ids_or_sequential()
    first_mass: dict = {}
    atom_list = []
    for i in range(len(structure)):
        fixed = [bool(structure._free[i, a] < 0.5) for a in range(3)]
        symbol = structure.names[i]
        mass = float(structure.mass[i])
        seen = first_mass.get(symbol)
        if seen is None:
            first_mass[symbol] = mass
        elif seen != mass:
            mass = seen
        kwargs = {}
        if structure.forces is not None:
            row = np.asarray(structure.forces[i], dtype=float).reshape(3)
            kwargs.update(fx=float(row[0]), fy=float(row[1]), fz=float(row[2]))
        if structure.velocities is not None:
            row = np.asarray(structure.velocities[i], dtype=float).reshape(3)
            kwargs.update(vx=float(row[0]), vy=float(row[1]), vz=float(row[2]))
        atom_list.append(
            readcon.Atom(
                symbol=symbol,
                x=float(structure.r[i][0]),
                y=float(structure.r[i][1]),
                z=float(structure.r[i][2]),
                fixed=fixed,
                atom_id=int(ids[i]),
                mass=mass,
                **kwargs,
            )
        )
    frame = readcon.ConFrame(
        cell=list(lengths),
        angles=list(angles),
        atoms=atom_list,
        prebox_header=["Generated by eOn", ""],
        postbox_header=None,
    )
    if structure.con_energy is not None:
        frame.set_energy(float(structure.con_energy))
    return readcon.write_con_string([frame])


def ref_pbc_rows(rows, box) -> np.ndarray:
    """``minimage.Cell.displacement`` of each row from the origin."""
    import minimage

    cell = minimage.Cell.from_vesin(np.asarray(box, dtype=float).tolist())
    rows = np.atleast_2d(np.asarray(rows, dtype=float))
    out = np.empty_like(rows)
    for k, row in enumerate(rows):
        out[k] = cell.displacement([0.0, 0.0, 0.0], row.tolist())
    return out


def ref_pair_lists(n, i, j):
    """Sorted unique non-self neighbours per atom, as plain ints."""
    sets = [set() for _ in range(n)]
    for a, b in zip(i.tolist(), j.tolist()):
        if int(a) != int(b):
            sets[int(a)].add(int(b))
    return [sorted(s) for s in sets]


def ref_neighbor_lists(p, cutoff):
    """Cutoff adjacency from the raw ``linkcell.pairs_within`` rows."""
    import linkcell

    r = np.ascontiguousarray(p.r, dtype=np.float64)
    box = np.ascontiguousarray(p.box, dtype=np.float64)
    i, j, _shift, _dist2 = linkcell.pairs_within(r, box, float(cutoff))
    return ref_pair_lists(
        len(p), np.asarray(np.from_dlpack(i)), np.asarray(np.from_dlpack(j))
    )


def ref_neighbor_list_vectors(p, cutoff):
    """Minimum-image vector to each neighbour, one pair at a time."""
    out = []
    for center, neighbours in enumerate(ref_neighbor_lists(p, cutoff)):
        out.append(
            [ref_pbc_rows(p.r[j] - p.r[center], p.box)[0] for j in neighbours]
        )
    return out


def ref_get_displacement(disp, atom_index):
    """``Displace.get_displacement`` as written one atom at a time.

    Draws from ``np.random`` in the order the method documents: one normal
    triple per displaced free atom, then the void-bias force from the
    neighbour vectors.
    """
    reactant = disp.reactant
    config = disp.config
    assert not config.random_mode
    nl = ref_neighbor_lists(reactant, disp.radius)
    if hasattr(atom_index, "__getitem__"):
        neighbours = [nl[i] for i in range(len(nl)) if i in atom_index]
        neighbours = sum(neighbours, [])
        neighbours = np.array(list(set(neighbours)), dtype=int)
        displaced = np.append(atom_index, neighbours)
    else:
        displaced = [atom_index] + nl[atom_index]
    norm = 0.0
    displacement = np.zeros(reactant.r.shape)
    while norm <= config.disp_min_norm:
        displacement = np.zeros(reactant.r.shape)
        for k in range(len(displaced)):
            if not reactant.atom_is_free()[displaced[k]]:
                continue
            displacement[displaced[k]] = np.random.normal(scale=disp.std_dev, size=3)
            if config.displace_1d:
                displacement[displaced[k]] = displacement[displaced[k]] * [1, 0, 0]
        norm = np.linalg.norm(displacement)
    vectors = ref_neighbor_list_vectors(reactant, disp.radius)
    force = np.zeros(reactant.r.shape)
    for atom in displaced:
        for vec in vectors[atom]:
            mag = np.linalg.norm(vec) + 1e-6
            force[atom] += -vec / (mag**3)
    void = force / np.linalg.norm(force)
    displacement = (
        disp.void_bias_fraction * norm * void
        + (1 - disp.void_bias_fraction) * displacement
    )
    moved = reactant.copy()
    moved.r += displacement
    displacement /= np.linalg.norm(displacement)
    return moved, displacement / np.linalg.norm(displacement)
