"""linkcell k-nearest pairs agree with vesin inside the cutoff.

The wire lists carry atom ids and per-axis fixed bits that ConFrame already
stores and the older Geometry struct dropped.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from eon.geometry import neighbor_list, neighbor_list_linkcell
from eon.geometry.wire import geometry_wire_lists
from eon.structure import Structure


def _line() -> Structure:
    p = Structure(4)
    p.r = np.array(
        [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [2.4, 0.0, 0.0], [3.6, 0.0, 0.0]],
        dtype=float,
    )
    p.box = np.eye(3) * 40.0
    p.free = np.array(
        [[1.0, 1.0, 1.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
    )
    p.names = ["H"] * 4
    p.mass = np.ones(4)
    p.periodic = np.array([True, True, True])
    p.atom_ids = np.array([4, 1, 3, 2], dtype=np.uint64)
    return p


def test_linkcell_pairs_match_vesin_inside_cutoff():
    p = _line()
    cutoff = 1.5
    vesin = {i: set(neigh) for i, neigh in enumerate(neighbor_list(p, cutoff))}
    linked = {
        i: set(neigh) for i, neigh in enumerate(neighbor_list_linkcell(p, cutoff))
    }
    assert linked == vesin


def test_geometry_wire_lists_keep_ids_forces_and_fixed_axes():
    forces = np.arange(12, dtype=float).reshape(4, 3)
    lists = geometry_wire_lists(_line(), forces=forces)
    assert lists["atomId"] == [4, 1, 3, 2]
    assert lists["forces"] == list(range(12))
    assert lists["fixedAxes"] == [0, 0b001, 0b010, 0b100]
    schema = Path(__file__).resolve().parents[1] / "schema" / "eon_job_result.capnp"
    text = schema.read_text(encoding="utf-8")
    assert "forces @8" in text
    assert "atomId @9" in text
    assert "fixedAxes @10" in text
