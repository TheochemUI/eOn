"""More shipped helpers: recycling suggestions, energy levels, and con text."""

from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from eon import atoms
from eon import fileio as io
from eon.displace import Displace, Random
from eon.escaperate import load_end_state_table, save_end_state_table, step
from eon.process_catalog import _saddle_text_with_mode, _suggestion_mode
from eon.recycling import Recycling
from eon.structure import Structure
from eon.superbasinscheme import EnergyLevel


def _pair(x1):
    atoms_ = Structure(2)
    atoms_.names = ["Pt", "Pt"]
    atoms_.mass[:] = 195.084
    atoms_.box = np.diag([20.0, 20.0, 20.0])
    atoms_.r[0] = [0.0, 0.0, 0.0]
    atoms_.r[1] = [x1, 0.0, 0.0]
    return atoms_


class _ProcState:
    def __init__(self, number, path, reactant, saddle, mode):
        self.number = number
        self.path = path
        self.reactant = reactant
        self.procs = {3: {"product": 1}}
        self.saddle = saddle
        self.mode = mode

    def load_process_table(self):
        return None

    def get_reactant(self):
        return self.reactant

    def get_process_saddle(self, pid):
        return self.saddle

    def get_process_mode(self, pid):
        return self.mode


def test_recycling_suggestion_moves_the_shifted_atom(tmp_path):
    ref = _pair(2.5)
    curr = _pair(2.5)
    curr.r[0, 0] = 0.4
    saddle = _pair(2.8)
    mode = np.zeros_like(ref.r)
    mode[1, 0] = 0.2
    ref_state = _ProcState(0, tmp_path, ref, saddle, mode)
    current = _ProcState(1, tmp_path / "cur", curr, saddle, mode)
    Path(current.path).mkdir()
    cfg = SimpleNamespace(
        comp_eps_r=0.2,
        recycling_active_region=1,
        comp_brute_neighbors=False,
    )
    recycler = Recycling(None, ref_state, current, 0.2, save=True, config=cfg)
    suggested, suggested_mode = recycler.make_suggestion()
    assert suggested is not None
    assert np.isfinite(suggested.r).all()
    assert np.isfinite(suggested_mode).all()
    assert recycler.process_number == 1
    assert recycler.make_suggestion() == (None, None)


def test_suggestion_mode_follows_the_saddle_frame(tmp_path):
    reactant = _pair(2.5)
    saddle = _pair(2.8)
    mode = saddle.r - reactant.r
    buf = StringIO()
    io.savecon(buf, saddle)
    text = _saddle_text_with_mode(buf.getvalue(), mode)
    assert "Pt" in text
    process = SimpleNamespace(mode=mode.reshape(-1))
    direction = _suggestion_mode(process, reactant, saddle, 0.0, text)
    assert direction.shape == saddle.r.shape
    assert np.isfinite(direction).all()
    assert np.linalg.norm(direction) > 0.0


def test_energy_level_rises_after_a_crossing(tmp_path):
    class _S:
        def __init__(self, number, energy, product, barrier, path):
            self.number = number
            self.energy = energy
            self.path = path
            self.procs = {1: {"product": product, "barrier": barrier}}

        def get_energy(self):
            return self.energy

        def get_process_table(self):
            return self.procs

    class _States:
        def get_num_states(self):
            return 0

    cfg = SimpleNamespace(
        comp_use_identical=False,
        sb_state_file="superbasin",
        sb_max_size=0,
    )
    scheme = EnergyLevel(tmp_path / "basins", _States(), 0.025, 0.05, config=cfg)
    start = _S(0, -1.0, 1, 0.2, tmp_path / "s0")
    end = _S(1, -0.9, 0, 0.2, tmp_path / "s1")
    start.path.mkdir()
    end.path.mkdir()
    scheme.register_transition(start, end)
    assert scheme.global_energy_min == -1.0
    assert scheme.levels[end] > end.get_energy()
    assert scheme.get_energy_increment(-0.9, -1.0, -0.8) > 0.0


def test_end_state_table_round_trip(tmp_path):
    path = tmp_path / "0" / "end_state_table"
    rows = [{"state": 2, "views": 3, "rate": 1.5e6, "time": 4.0, "process_id": 9}]
    save_end_state_table(path, rows)
    loaded = load_end_state_table(path)
    assert loaded[0]["state"] == 2
    assert loaded[0]["process_id"] == 9
    assert loaded[0]["rate"] == 1.5e6


def test_escape_step_records_the_product(tmp_path):
    class _State:
        def __init__(self, number):
            self.number = number
            self.energy = -1.0

        def get_process(self, pid):
            return {"barrier": 0.2, "rate": 1.0}

        def get_energy(self):
            return self.energy

        def zero_time(self):
            return None

    class _States:
        def __init__(self):
            self.product = _State(1)

        def get_product_state(self, reactant, proc_id):
            assert reactant == 0 and proc_id == 4
            return self.product

    cfg = SimpleNamespace(path_results=str(tmp_path))
    current, previous = step(
        0.0,
        _State(0),
        _States(),
        {"process_id": 4, "time": 2.5},
        config=cfg,
    )
    assert current.number == 1
    assert previous.number == 0
    text = Path(tmp_path, "dynamics.txt").read_text()
    assert "4" in text


def test_one_dimensional_random_mode_stays_on_x():
    reactant = _pair(2.5)
    cfg = SimpleNamespace(
        disp_magnitude=0.1,
        disp_radius=4.0,
        disp_min_norm=0.0,
        displace_1d=True,
        random_mode=True,
        comp_brute_neighbors=False,
        comp_neighbor_cutoff=3.5,
        comp_use_covalent=False,
        comp_covalent_scale=1.3,
    )
    np.random.seed(1)
    kind = Random(reactant, 0.1, 4.0, config=cfg)
    _moved, mode = kind.make_displacement()
    assert np.allclose(mode[:, 1:], 0.0)
    assert np.isfinite(mode[:, 0]).all()


def test_empty_hole_reverts_to_the_full_epicenter_list():
    disp = Displace.__new__(Displace)
    disp.hole_epicenters = [9]
    assert disp.filter_epicenters([0, 1]) == [0, 1]
    assert disp.filter_epicenters([9, 1]) == [9]


def test_saved_con_reloads_through_readcon(tmp_path):
    atoms_ = _pair(2.5)
    path = tmp_path / "pair.con"
    io.savecon(str(path), atoms_)
    loaded = io.loadcons(str(path))
    assert len(loaded) == 1
    assert loaded[0].names == ["Pt", "Pt"]
    assert np.allclose(loaded[0].r, atoms_.r)
    neighbors = atoms.neighbor_list(atoms_, 3.5)
    assert 1 in neighbors[0]
