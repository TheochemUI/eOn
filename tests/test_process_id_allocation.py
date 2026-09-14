"""Process ids must be content-addressed (xxh64), not len(procs) / max+1.

Duplicate processtable rows collapse in a dict keyed by id. A counter based
on len(procs) freezes and reuses ids (eOn-cj3o); max(id)+1 still couples
identity to registration order. xxh64 of the process payload does not.
"""

from __future__ import annotations

import os
from unittest import mock

import pytest

from eon.akmcstate import AKMCState
from eon.process_id import allocate_unique_process_id, process_id_from_parts


def test_process_id_from_parts_is_stable_and_positive():
    a = process_id_from_parts(b"saddle", b"barrier=0.1")
    b = process_id_from_parts(b"saddle", b"barrier=0.1")
    c = process_id_from_parts(b"saddle", b"barrier=0.2")
    assert a == b
    assert a != c
    assert a >= 0
    assert a < 2**63


def test_allocate_unique_salts_when_occupied():
    base = process_id_from_parts(b"payload")
    occupied = {base: {}}
    free = allocate_unique_process_id(occupied, b"payload")
    assert free != base
    assert free == process_id_from_parts(b"payload", salt=1)


def _make_akmc_state(tmp_path, table_body: str = ""):
    from eon import fileio as io

    state_dir = tmp_path / "state_0"
    state_dir.mkdir()
    (state_dir / "procdata").mkdir()
    proctable = state_dir / "processtable"
    proctable.write_text(AKMCState.processtable_header + table_body)

    state = object.__new__(AKMCState)
    state.config = mock.Mock()
    state.statelist = mock.Mock()
    state.statelist.kT = 0.025
    state.path = str(state_dir)
    state.number = 0
    state.procs = None
    state.proc_repeat_count = None
    state.procdata_path = os.path.join(state.path, "procdata")
    state.reactant_path = os.path.join(state.path, "reactant.con")
    state.proctable_path = str(proctable)
    state.search_result_path = os.path.join(state.path, "search_results.txt")
    state.tar_path = os.path.join(state.path, "procdata.tar")
    state.info = io.ini(os.path.join(state.path, "info"))
    state.info.set("MetaData", "kT", 0.025)
    return state


# Shape of the Cu V/SIA observation: 8 rows, 3 distinct ids.
_DUP_TABLE = """\
      0          1.00000 1.00000e+12         1          0.50000       1.00000e+12  0.74164  3.47609e-01       0
      1          1.10000 1.00000e+12         2          0.60000       1.00000e+12  0.70427  1.47483e+00       0
      2          1.20000 1.00000e+12         3          0.70000       1.00000e+12  0.70086  1.68270e+00       0
      2          0.90000 1.00000e+12         4          0.10000       1.00000e+12  0.10270  1.88217e+10       0
      2          1.15000 1.00000e+12         5          0.65000       1.00000e+12  0.68879  2.68439e+00       0
      2          0.91000 1.00000e+12         6          0.11000       1.00000e+12  0.10171  1.95599e+10       0
      2          1.14000 1.00000e+12         7          0.64000       1.00000e+12  0.65673  9.27689e+00       0
      2          0.92000 1.00000e+12         8          0.12000       1.00000e+12  0.10271  1.88176e+10       0
"""


def test_get_num_procs_is_distinct_count_not_row_count(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    assert state.get_num_procs() == 3
    assert sorted(state.procs.keys()) == [0, 1, 2]
    assert state.procs[2]["barrier"] == pytest.approx(0.10271)


def test_allocate_does_not_reuse_existing_ids(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    # Old bug: next id was 3 = len(procs), reusing nothing but freeze was at 2.
    # Content hash must never land on 0,1,2 for this unrelated payload, or if
    # it does, salt until free.
    pid = state.allocate_process_id(b"akmc-forward", b"saddle-xyz", b"barrier=0.55")
    assert pid not in (0, 1, 2)
    assert pid not in state.procs


def test_allocate_stable_for_same_payload(tmp_path):
    state = _make_akmc_state(tmp_path, "")
    a = state.allocate_process_id(b"akmc-forward", b"saddle-A", b"barrier=0.10")
    # Force reload empty table still yields same hash id.
    state.procs = None
    b = state.allocate_process_id(b"akmc-forward", b"saddle-A", b"barrier=0.10")
    assert a == b


def test_append_refuses_duplicate_id(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    state.load_process_table()
    with pytest.raises(RuntimeError, match="refusing to clobber process id 2"):
        state.append_process_table(
            id=2,
            saddle_energy=1.0,
            prefactor=1e12,
            product=9,
            product_energy=0.0,
            product_prefactor=1e12,
            barrier=0.5,
            rate=1.0,
            repeats=0,
        )


def test_append_with_allocated_id_does_not_clobber(tmp_path):
    state = _make_akmc_state(tmp_path, _DUP_TABLE)
    nid = state.allocate_process_id(b"akmc-forward", b"new-saddle", b"barrier=0.20")
    assert nid not in state.procs
    state.append_process_table(
        id=nid,
        saddle_energy=1.0,
        prefactor=1e12,
        product=99,
        product_energy=0.0,
        product_prefactor=1e12,
        barrier=0.2,
        rate=1e5,
        repeats=0,
    )
    assert state.get_num_procs() == 4
    assert 2 in state.procs
    assert state.procs[nid]["product"] == 99


def test_forward_and_reverse_tags_differ():
    saddle = b"same-saddle-geometry"
    fwd = process_id_from_parts(b"akmc-forward", saddle, b"barrier=0.10")
    rev = process_id_from_parts(
        b"akmc-reverse", saddle, b"from-state-0", b"forward-proc-1"
    )
    assert fwd != rev


def test_recycling_indexes_sparse_proc_ids(tmp_path):
    """Recycling must map ordinal process_number -> real process table keys."""
    from eon.recycling import Recycling

    # Minimal fake state with sparse content-like keys
    class FakeState:
        def __init__(self, number, procs):
            self.number = number
            self.path = str(tmp_path / f"s{number}")
            os.makedirs(self.path, exist_ok=True)
            self.procs = procs
            self._loaded = True

        def load_process_table(self):
            return

        def get_reactant(self):
            import numpy as np

            class _R:
                def __init__(self):
                    self.r = np.zeros((3, 3))
                    self.box = np.eye(3) * 10.0

                def __len__(self):
                    return 3

                def copy(self):
                    o = _R()
                    o.r = self.r.copy()
                    return o

            return _R()

        def get_process_saddle(self, pid):
            import numpy as np
            from types import SimpleNamespace
            assert pid in self.procs, f"recycling used ordinal as id: {pid}"
            r = np.zeros((3, 3))
            return SimpleNamespace(r=r)

        def get_process_mode(self, pid):
            import numpy as np
            assert pid in self.procs
            return np.zeros((3, 3))

    sparse = {10**12 + 7: {}, 99: {}, 5: {}}
    ref = FakeState(0, sparse)
    cur = FakeState(1, {})
    states = mock.Mock()
    config = mock.Mock()
    config.comp_eps_r = 0.2
    config.recycling_active_region = 5.0
    config.saddle_method = "min_mode"

    with mock.patch("eon.recycling.atoms.get_process_atoms", return_value=[0, 1, 2]):
        with mock.patch("eon.recycling.atoms.per_atom_norm", return_value=[0.0, 0.0, 0.0]):
            rec = Recycling(states, ref, cur, move_distance=0.1, config=config, save=False)
    assert rec.num_procs == 3
    assert rec.proc_ids == sorted(sparse.keys())
    # make_suggestion must call get_process_saddle with real keys
    rec.make_suggestion()
    assert rec.process_number == 1
