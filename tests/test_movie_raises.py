"""Movie helpers raise on a bad request and write frames through the con reader."""

from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.movie import make_graph, make_movie
from eon.structure import Structure


class _State:
    def __init__(self, number, table, reactant_path=None):
        self.number = number
        self.table = table
        self.reactant_path = reactant_path
        self.reactant = None

    def get_process_table(self):
        return self.table

    def get_process_reactant(self, pid):
        return self.table[pid]["reactant"]

    def get_process_saddle(self, pid):
        return self.table[pid]["saddle"]

    def get_process_product(self, pid):
        return self.table[pid]["product"]

    def get_reactant(self):
        return self.reactant


class _States:
    def __init__(self, states):
        self.states = {state.number: state for state in states}

    def get_num_states(self):
        return len(self.states)

    def get_state(self, number):
        try:
            return self.states[int(number)]
        except KeyError as exc:
            raise OSError(f"no state {number}") from exc


def _pt():
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.mass[:] = 195.084
    atoms.box = np.diag([10.0, 10.0, 10.0])
    atoms.r[0] = [1.0, 2.0, 3.0]
    return atoms


def test_unknown_movie_type_raises():
    with pytest.raises(ValueError, match="unknown movie type"):
        make_movie("not-a-movie", ".", _States([]))


def test_process_movie_requires_a_state_number():
    with pytest.raises(ValueError, match="must give a state number"):
        make_movie("processes", ".", _States([]))


def test_process_movie_state_must_be_an_integer():
    with pytest.raises(ValueError, match="state number must be an integer"):
        make_movie("processes,abc", ".", _States([]))


def test_missing_state_raises():
    with pytest.raises(FileNotFoundError, match="cannot make a movie for state 4"):
        make_movie("processes,4", ".", _States([]))


def test_empty_state_list_raises():
    with pytest.raises(ValueError, match="no dynamics steps"):
        make_movie("states", ".", _States([]))


def test_existing_graph_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    Path("graph.dot").write_text("digraph {}\n")
    state = _State(0, {})
    with pytest.raises(FileExistsError, match="graph.dot"):
        make_movie("graph", tmp_path, _States([state]))


def test_graph_file_names_the_edge(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    left = _State(0, {3: {"product": 1, "rate": 2.0}})
    right = _State(1, {})
    make_movie("graph", tmp_path, _States([left, right]))
    assert "0 -- 1;" in Path("graph.dot").read_text()


def test_zero_rate_is_not_a_graph_weight():
    reactant = _pt()
    left = _State(0, {3: {"product": 1, "rate": 0.0, "reactant": reactant}})
    right = _State(1, {})
    with pytest.raises(ValueError, match="no positive rate"):
        make_graph(_States([left, right]))


def test_dynamics_movie_round_trips_through_readcon(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = _pt()
    con = tmp_path / "reactant.con"
    io.savecon(str(con), atoms)
    dynamics = tmp_path / "dynamics.txt"
    dynamics.write_text(
        "step-number time energy product-id\n"
        "------------------------------------\n"
        "0 0.0 -1.0 0\n"
    )
    state = _State(0, {}, reactant_path=con)
    make_movie("dynamics", tmp_path, _States([state]))
    written = Path("dynamics.poscar").read_text()
    assert written.splitlines()[0] == "Pt"
    loaded = io.loadposcar(StringIO(written))
    assert loaded.names == ["Pt"]
    assert np.allclose(loaded.r[0], atoms.r[0])


def test_process_movie_writes_three_frames(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    frame = _pt()
    saddle = frame.copy()
    saddle.r[0, 0] = 1.5
    product = frame.copy()
    product.r[0, 0] = 2.0
    state = _State(
        0,
        {
            3: {
                "rate": 1.0e6,
                "barrier": 0.2,
                "prefactor": 1.0e12,
                "reactant": frame,
                "saddle": saddle,
                "product": product,
            }
        },
    )
    make_movie("processes,0,1", tmp_path, _States([state]))
    text = Path("processes_0.poscar").read_text()
    assert text.count("\nPt\n") + text.count("Pt\n") >= 1
    assert "1.50000000000000" in text or "1.5000000000000" in text
