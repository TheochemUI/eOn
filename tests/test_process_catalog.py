"""Process catalog: frames in readcon-db, rows in amsel.KdbStore."""

from __future__ import annotations

import inspect
import logging
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon import process_catalog as catalog
from eon.explorer import MinModeExplorer
from eon.structure import Structure

ROOT = Path(__file__).resolve().parents[1]


def _cfg(tmp_path: Path, **overrides):
    class Cfg:
        kdb_on = True
        kdb_only = False
        kdb_nf = 0.15
        kdb_dc = 0.3
        kdb_mac = 0.7
        kdb_path = tmp_path / "kdb"
        kdb_scratch_path = tmp_path / "kdbscratch"
        main_temperature = 300.0
        recycling_on = False
        saddle_method = "min_mode"

    cfg = Cfg()
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def _pair(shift):
    atoms = Structure(2)
    atoms.names = ["Pt", "Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([20.0, 20.0, 20.0])
    atoms.r[0] = [1.0, 2.0, 3.0]
    atoms.r[1] = [2.4, 2.0, 3.0]
    atoms.r[0] = atoms.r[0] + np.asarray(shift, dtype=float)
    return atoms


def _write_process(root: Path, reactant, saddle, product, mode, barrier=0.3007):
    root.mkdir(parents=True, exist_ok=True)
    proc = root / "procdata"
    proc.mkdir(exist_ok=True)
    io.savecon(str(root / "reactant.con"), reactant)
    io.savecon(str(proc / "reactant_1.con"), reactant)
    io.savecon(str(proc / "saddle_1.con"), saddle)
    io.savecon(str(proc / "product_1.con"), product)
    io.save_mode(str(proc / "mode_1.dat"), mode)

    class State:
        number = int(root.name.split("_")[-1]) if root.name.startswith("state_") else 0
        path = str(root)
        reactant_path = str(root / "reactant.con")
        procs = {
            1: {
                "barrier": barrier,
                "prefactor": 1.0e12,
            }
        }

        def get_reactant(self):
            return io.loadcon(self.reactant_path)

        def proc_reactant_path(self, pid):
            return str(proc / f"reactant_{pid}.con")

        def proc_saddle_path(self, pid):
            return str(proc / f"saddle_{pid}.con")

        def proc_product_path(self, pid):
            return str(proc / f"product_{pid}.con")

        def get_process_mode(self, pid):
            return io.load_mode(str(proc / f"mode_{pid}.dat"))

    return State()


def test_old_kdb_module_is_gone():
    assert not (ROOT / "eon" / "eon_kdb.py").exists()
    tree = "\n".join(path.read_text() for path in (ROOT / "eon").glob("*.py"))
    assert "from kdb import" not in tree
    assert "get_params" not in (ROOT / "eon" / "process_catalog.py").read_text()
    assert "Python module kdb not found" not in tree
    assert "rmtree" not in (ROOT / "eon" / "process_catalog.py").read_text()


def test_insert_is_not_gated_on_confidence():
    explore = (ROOT / "eon" / "explorer.py").read_text()
    assert "Adding relevant processes to kinetic database" not in explore
    assert "eon_kdb" not in explore
    state = (ROOT / "eon" / "akmcstate.py").read_text()
    add = state.split("def add_process", 1)[1].split("\n    def ", 1)[0]
    assert "process_catalog.insert" in add
    assert "akmc_confidence" not in add


def test_corpus_methods_take_no_barrier_or_mode():
    readcon_db = pytest.importorskip("readcon_db")
    from eon import concorpus

    for fn in (
        readcon_db.ConCorpus.append_trajectory_str,
        concorpus.store_frame_text,
        concorpus.load_frame_text,
    ):
        names = inspect.signature(fn).parameters
        assert "barrier" not in names
        assert "mode" not in names


def test_canary_runs_server_and_client():
    script = (ROOT / "eon" / "tests" / "canary" / "kdb" / "test.py").read_text()
    ini = (ROOT / "eon" / "tests" / "canary" / "kdb" / "config.ini").read_text()
    assert "akmc.py" not in script
    assert "eon-server" in script
    assert "no search has type kdb" in script
    assert "client_path = eonclient" in ini
    assert "client/client" not in ini


def test_failed_query_logs_nf_and_does_not_mark_the_state(tmp_path, caplog, monkeypatch):
    monkeypatch.setattr(catalog, "_open_store", lambda config: None)
    cfg = _cfg(tmp_path)
    caplog.set_level(logging.INFO, logger="kdb")
    ok = catalog.query(type("S", (), {"number": 4})(), cfg)
    assert ok is False
    assert "nf=0.15" in caplog.text
    assert "path=" in caplog.text
    assert str(cfg.kdb_path) in caplog.text
    assert "Python module kdb not found" not in caplog.text
    assert not (tmp_path / "kdbscratch" / "queried").exists()


def test_kdb_only_does_not_make_a_random_search(tmp_path):
    cfg = _cfg(tmp_path, kdb_only=True)
    reactant = _pair([0, 0, 0])
    state = _write_process(
        tmp_path / "state_0",
        reactant,
        reactant,
        reactant,
        np.zeros((2, 3)),
    )

    class Boom:
        def make_displacement(self):
            raise AssertionError("random displacement")

    expl = MinModeExplorer.__new__(MinModeExplorer)
    expl.config = cfg
    expl.state = state
    expl.reactant = reactant
    expl.displace = Boom()
    displacement, mode, kind = MinModeExplorer.generate_displacement(expl)
    assert displacement is None
    assert mode is None
    assert kind == "kdb-empty"


def test_first_find_is_suggested_before_confidence(tmp_path, caplog):
    pytest.importorskip("amsel")
    pytest.importorskip("readcon_db")
    from amsel import KdbStore

    (tmp_path / "config.ini").write_text("[Main]\njob = akmc\n")
    reactant = _pair([0.0, 0.0, 0.0])
    saddle = _pair([0.5, 0.0, 0.0])
    product = _pair([1.0, 0.0, 0.0])
    mode = saddle.r - reactant.r
    state0 = _write_process(tmp_path / "state_0", reactant, saddle, product, mode)
    other = _pair([5.0, 0.0, 0.0])
    state1 = _write_process(
        tmp_path / "state_1",
        other,
        other,
        other,
        np.zeros((2, 3)),
        barrier=1.0,
    )
    cfg = _cfg(tmp_path)
    assert catalog.insert(state0, 1, cfg)
    store = KdbStore(str(cfg.kdb_path))
    hits = store.lookup(catalog.env_hash(reactant))
    assert len(hits) == 1
    assert hits[0].barrier_ev == pytest.approx(0.3007)
    assert bytes(hits[0].saddle_con) == b""
    assert len(bytes(hits[0].saddle_frame_key)) == 12
    assert list(hits[0].mode)[0] == pytest.approx(0.5)
    traj_id, frame_idx = catalog.unpack_frame_key(bytes(hits[0].saddle_frame_key))
    from eon.concorpus import corpus_dir, load_frame_text

    text = load_frame_text(corpus_dir(tmp_path / "state_0" / "reactant.con"), traj_id, frame_idx)
    assert text
    stored = io.loadcon(io.StringIO(text))
    assert np.allclose(stored.r, saddle.r)

    caplog.set_level(logging.INFO, logger="kdb")
    assert catalog.query(state0, cfg)
    assert catalog.query(state1, cfg)
    saddle_file = tmp_path / "kdbscratch" / "kdbmatches" / "state_0" / "SADDLE_0"
    assert saddle_file.is_file()
    assert "nf=0.15" in caplog.text
    assert (tmp_path / "kdb" / "data.mdb").is_file()

    class Displace:
        def __init__(self):
            self.called = False

        def make_displacement(self):
            self.called = True
            return reactant, np.zeros((2, 3))

    expl = MinModeExplorer.__new__(MinModeExplorer)
    expl.config = cfg
    expl.state = state0
    expl.reactant = reactant
    expl.displace = Displace()
    caplog.clear()
    caplog.set_level(logging.INFO)
    displacement, mode_out, kind = MinModeExplorer.generate_displacement(expl)
    assert kind == "kdb"
    assert expl.displace.called is False
    assert "Made a KDB suggestion" in caplog.text
    assert "0.3007" in caplog.text
    assert "readcon.db" in caplog.text
    assert mode_out.shape == (2, 3)
    assert np.allclose(displacement.r, saddle.r, atol=1e-5)
    assert not saddle_file.is_file()
    again, _, kind_again = MinModeExplorer.generate_displacement(expl)
    assert kind_again == "random"
    assert expl.displace.called is True
    assert again is reactant
