"""Small branches that still sit in the coverage table."""

import pickle
import sys
from io import StringIO
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from eon import atoms
from eon import fileio as io
from eon.amsel_superbasin_gate import amsel_discover_exit
from eon.explorer import ServerMinModeExplorer
from eon.process_catalog import _mode_cosine, _open_store, query
from eon.structure import Structure
from tests.test_library_bodies import _Comm, _atoms, _config, _states


def test_appended_frames_and_a_rejected_box(tmp_path):
    left = _atoms(2.5)
    right = _atoms(3.0)
    path = tmp_path / "pair.con"
    io.savecon(str(path), left)
    io.savecon(str(path), right, w="a")
    loaded = io.loadcons(str(path))
    assert len(loaded) == 2
    gz = tmp_path / "pair.con.gz"
    io.savecon(str(gz), left)
    io.savecon(str(gz), right, w="a")
    assert len(io.loadcons(str(gz))) == 2
    empty = tmp_path / "empty.con"
    empty.write_bytes(b"")
    io.savecon(str(empty), left, w="a")
    assert len(io.loadcons(str(empty))) == 1

    shifted = left.copy()
    shifted.box[0, 0] += 1.0
    assert atoms.identical(left, shifted, 0.1) is False
    assert atoms.match(left, Structure(1), 0.1, 3.3, True) is False
    assert atoms.atomic_number("X") == 0
    assert atoms.symbol_for_z(0) == "Xx"
    assert atoms.elements[200]["symbol"] == "Xx"
    io.savecon(str(tmp_path / "a.con"), left)
    io.savecon(str(tmp_path / "b.con"), left)
    assert (
        atoms.point_energy_match(
            str(tmp_path / "a.con"),
            -1.0,
            str(tmp_path / "b.con"),
            -3.0,
            0.1,
            0.2,
            3.3,
        )
        is False
    )
    moved = left.copy()
    moved.r += 0.4
    assert atoms.match(
        left,
        moved,
        0.1,
        3.3,
        False,
        remove_translation=True,
    )


def test_a_broken_store_and_a_reloaded_explorer(tmp_path, monkeypatch):
    cfg, states, state, _product, _proc, reactant = _states(tmp_path)
    cfg.kdb_path = str(tmp_path / "kdb")
    cfg.kdb_scratch_path = str(tmp_path / "scratch")
    module = ModuleType("amsel")

    class _BrokenStore:
        def __init__(self, path):
            raise OSError("store down")

    module.KdbStore = _BrokenStore
    monkeypatch.setitem(sys.modules, "amsel", module)
    assert _open_store(cfg) is None
    assert query(state, cfg) is False
    assert _mode_cosine(np.zeros((0, 3)), Structure(0), Structure(0)) == 0.0

    cfg.akmc_server_side_process_search = True
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_type = "local_inprocess"
    cfg.comm_job_buffer_size = 0
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(tmp_path)
    (tmp_path / "explorer.pickle").write_bytes(pickle.dumps({"search_id": 7}))
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: _Comm())
    explorer = ServerMinModeExplorer(states, state, state, config=cfg)
    assert explorer.search_id == 7

    io.savecon(str(tmp_path / "reactant.con"), reactant)
    import readcon

    frames = readcon.read_con(str(tmp_path / "reactant.con"))
    rebuilt = Structure.from_conframe(frames[0])
    assert len(rebuilt) == len(reactant)
    assert rebuilt.names[0] == "Pt"


def test_a_bad_exit_kernel_returns_nothing(tmp_path, monkeypatch):
    cfg, states, state, _product, _proc, _reactant = _states(tmp_path)
    cfg.amsel_e_min_init = 0.01
    cfg.debug_use_mean_time = False
    module = ModuleType("amsel")

    def discover_decide_status(*_args, **_kwargs):
        return ("accepted", [0], [1], [])

    def fpta(_transient, absorbing, _rates, _entry, _draw):
        return float("nan"), [1.0]

    def mrm(_transient, absorbing, _rates, _entry):
        return -1.0, [1.0] * len(absorbing), None

    module.discover_decide_status = discover_decide_status
    module.fpta = fpta
    module.mrm = mrm
    monkeypatch.setitem(sys.modules, "amsel", module)
    assert amsel_discover_exit(state, states.get_state, cfg, uniform=lambda: 0.2) is None
