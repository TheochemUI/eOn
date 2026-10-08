"""Field ids, atom-list scripts, and movie files the coverage run had not called."""

import logging
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon._utils import ScriptConfig, gen_ids_from_con
from eon.akmcstate import AKMCState
from eon.communicator_inprocess import LocalInProcess, _LazyCon, _job_result
from eon.config import ConfigClass
from eon.displace import DisplacementManager
from eon.explorer import (
    ClientMinModeExplorer,
    ServerMinModeExplorer,
    get_minmodexplorer,
)
from eon.movie import Graph, make_movie, priorityDictionary
from eon.params_ssot import all_field_ids, has_field, normalize_yaml_default
from eon.state import State
from eon.statelist import StateList
from eon.superbasin import Superbasin
from eon.atoms import point_energy_match
from tests.test_library_bodies import _atoms, _config, _states
from tests.test_movie_raises import _State, _States, _pt


def test_field_ids_skip_nested_scalars():
    ids = all_field_ids()
    assert ids
    assert all("." in item for item in ids)
    section, key = ids[0].split(".", 1)
    assert has_field(section, key) is True
    assert has_field(section, "nope") is False
    assert normalize_yaml_default(3) == 3
    assert normalize_yaml_default({"other": 1}) == {"other": 1}


def test_empty_atom_script_returns_no_indices(tmp_path):
    script = tmp_path / "empty.py"
    script.write_text("")
    sconf = ScriptConfig(script, tmp_path / "scratch", tmp_path)
    result = gen_ids_from_con(sconf, _atoms(2.5), logging.getLogger("atom-list"))
    assert result == ""


def test_failing_atom_script_exits(tmp_path):
    script = tmp_path / "bad.py"
    script.write_text("import sys\nsys.exit(1)\n")
    sconf = ScriptConfig(script, tmp_path / "scratch", tmp_path)
    with pytest.raises(SystemExit) as caught:
        gen_ids_from_con(sconf, _atoms(2.5), logging.getLogger("atom-list"))
    assert caught.value.code == 1


def test_missing_script_workdir_exits(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    script = root / "ok.py"
    script.write_text("print('1 2')\n")
    sconf = ScriptConfig(script, root / "scratch", root)
    sconf.root_path = tmp_path / "missing-root"
    with pytest.raises(SystemExit) as caught:
        gen_ids_from_con(sconf, _atoms(2.5), logging.getLogger("atom-list"))
    assert caught.value.code == 1


def test_a_second_config_read_returns_none(tmp_path):
    cfg = _config(tmp_path)
    assert cfg.init(str(tmp_path / "config.ini")) is None


def test_no_config_file_exits(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init("")
    assert caught.value.code == 2


def test_a_foreign_config_continues_when_asked(tmp_path, monkeypatch):
    other = tmp_path / "elsewhere"
    other.mkdir()
    root = tmp_path / "run"
    root.mkdir()
    monkeypatch.chdir(other)
    (other / "config.ini").write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "[Paths]",
                "main_directory = %s" % root,
                "results = %s" % root,
                "states = %s" % (root / "states"),
                "",
            ]
        )
    )
    monkeypatch.setattr("builtins.input", lambda _prompt: "y")
    cfg = ConfigClass()
    cfg.init("")
    assert Path(cfg.path_root).resolve() == root.resolve()


def test_config_restores_a_saved_generator(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "config.ini").write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "random_seed = 4",
                "[Paths]",
                "main_directory = %s" % tmp_path,
                "results = %s" % tmp_path,
                "states = %s" % (tmp_path / "states"),
                "",
            ]
        )
    )
    np.random.seed(7)
    io.save_prng_state(str(tmp_path / "prng.pkl"))
    expected = np.random.random()
    np.random.seed(0)
    cfg = ConfigClass()
    cfg.init(str(tmp_path / "config.ini"))
    assert cfg.main_random_seed == 4
    assert np.random.random() == expected


def test_states_movie_writes_one_frame(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = _pt()
    con = tmp_path / "reactant.con"
    io.savecon(str(con), atoms)
    make_movie("states", tmp_path, _States([_State(0, {}, reactant_path=con)]))
    assert Path("states.poscar").read_text().splitlines()[0] == "Pt"


def test_fastest_movies_name_their_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _cfg, states, _state, _product, _proc, _reactant = _states(tmp_path)
    make_movie("fastestpath", tmp_path, states)
    make_movie("fastestfullpath", tmp_path, states)
    assert "Pt" in Path("fastestpath.poscar").read_text()
    assert "Pt" in Path("fastestfullpath.poscar").read_text()


def test_movie_replaces_an_existing_frame_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    atoms = _pt()
    con = tmp_path / "reactant.con"
    io.savecon(str(con), atoms)
    (tmp_path / "dynamics.txt").write_text(
        "step-number time energy product-id\n"
        "------------------------------------\n"
        "0 0.0 -1.0 0\n"
    )
    Path("dynamics.poscar").write_text("old\n")
    movies = Path("movies")
    movies.mkdir()
    (movies / ("dynamics.poscar.%010d" % 0)).write_text("old\n")
    states = _States([_State(0, {}, reactant_path=con)])
    make_movie("dynamics", tmp_path, states)
    assert Path("dynamics.poscar").read_text().splitlines()[0] == "Pt"
    make_movie("dynamics", tmp_path, states, separate_files=True)
    written = (movies / ("dynamics.poscar.%010d" % 0)).read_text()
    assert written.splitlines()[0] == "Pt"


def test_graph_inserts_a_missing_node_and_rebuilds_priorities():
    graph = Graph("pair")
    assert repr(graph) == "{}"
    left = _State(0, {})
    right = _State(1, {})
    graph.add_edge(right, left, weight=0.5)
    assert graph.neighbors(right) == [left]
    assert graph.neighbors(left) == []
    assert graph.neighbors(_State(2, {})) is None
    assert "0 -- 1;" in graph.dot()
    queue = priorityDictionary()
    queue["a"] = 3.0
    queue["a"] = 1.0
    queue["a"] = 2.0
    assert queue["a"] == 2.0
    assert queue.smallest() == "a"


def test_constructors_require_a_config(tmp_path):
    with pytest.raises(TypeError, match="ConfigClass"):
        State(str(tmp_path), 0, None)
    with pytest.raises(TypeError, match="ConfigClass"):
        AKMCState(str(tmp_path), 0, None)
    with pytest.raises(TypeError, match="ConfigClass"):
        StateList(State, config=None)
    with pytest.raises(TypeError, match="ConfigClass"):
        DisplacementManager(_atoms(2.5), [], None)
    with pytest.raises(ValueError, match="list of states"):
        Superbasin(tmp_path, 1)
    with pytest.raises(TypeError, match="ConfigClass"):
        Superbasin(tmp_path, 1, state_list=[])
    with pytest.raises(TypeError, match="ConfigClass"):
        from eon.explorer import Explorer

        Explorer(config=None)
    with pytest.raises(TypeError, match="ConfigClass"):
        LocalInProcess(str(tmp_path))


def test_server_explorer_and_a_cancelled_result(tmp_path):
    cfg = _config(tmp_path)
    cfg.akmc_server_side_process_search = True
    assert get_minmodexplorer(cfg) is ServerMinModeExplorer
    cfg.akmc_server_side_process_search = False
    assert get_minmodexplorer(cfg) is ClientMinModeExplorer
    cancelled = _job_result(1, 0.1, 2, "point", cancelled=True)
    assert cancelled["termination_reason_text"] == "cancelled"
    lazy = _LazyCon(_pt())
    assert lazy.seek(0) == 0
    assert "Pt" in lazy.read()


def test_inprocess_queue_is_idle(tmp_path):
    cfg = _config(tmp_path)
    comm = LocalInProcess(cfg.path_scratch, config=cfg)
    assert comm.get_queue_size() == 0
    assert comm.get_number_in_progress() == 0
    assert comm.cancel_state(0) == 0


def test_close_energies_of_different_geometries_do_not_match(tmp_path):
    left = _pt()
    right = _pt()
    right.r[0, 0] = 4.0
    left_path = tmp_path / "left.con"
    right_path = tmp_path / "right.con"
    io.savecon(str(left_path), left)
    io.savecon(str(right_path), right)
    assert (
        point_energy_match(
            str(left_path), -1.0, str(right_path), -1.0, 0.01, 0.1, 3.0
        )
        is False
    )
