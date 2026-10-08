"""Table types, a repeated saddle result, and the displacement benchmark."""

import runpy
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pytest

from eon import fileio as io
from eon.explorer import ServerMinModeExplorer
from tests.test_library_bodies import _Comm, _atoms, _buf, _config, _mode, _result, _states


def test_table_types_and_a_masked_mode(tmp_path):
    path = tmp_path / "rates.tbl"
    table = io.Table(str(path), ["state", "rate"])
    table.add_row({"state": 1, "rate": 1.0})
    table.add_row({"state": 2, "rate": 3.0})
    assert table.find_value("rate", min) == 1.0
    assert table.find_row("rate", min)["state"] == 1
    assert table.get_row("state", 9) is None
    with pytest.raises(io.TableException, match="Type mismatch"):
        table.add_row({"state": 3, "rate": 1})
    mode = np.ones((2, 3))
    buf = StringIO()
    io.save_mode(buf, mode, free=np.array([1.0, 0.0]))
    text = buf.getvalue().splitlines()
    assert text[1].startswith("0 ")
    with pytest.raises(OSError, match="Malformed"):
        io.load_mode(StringIO("1.0 0.0 0.0\n1 2\n"))


def test_repeat_result_is_registered_from_the_saddle(tmp_path, monkeypatch):
    cfg, states, state, _product, proc, reactant = _states(tmp_path)
    cfg.akmc_server_side_process_search = True
    cfg.recycling_on = False
    cfg.kdb_on = False
    cfg.comm_type = "local_inprocess"
    cfg.comm_job_buffer_size = 1
    cfg.process_search_minimization_offset = 0.1
    Path(cfg.path_incomplete).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_scratch).mkdir(parents=True, exist_ok=True)
    Path(cfg.path_jobs_in).mkdir(parents=True, exist_ok=True)
    ini = Path(cfg.config_path)
    ini.write_text(ini.read_text() + "\n[Saddle Search]\nmethod = min_mode\n")
    comm = _Comm()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("eon.communicator.get_communicator", lambda config: comm)
    explorer = ServerMinModeExplorer(states, state, state, config=cfg)
    explorer.make_jobs()
    wuid = next(iter(explorer.wuid_to_search_id))
    saddle = io.loadcon(state.proc_saddle_path(proc))
    comm.results = [
        {
            "name": "0_%s" % wuid,
            "results": {
                "job_type": "saddle_search",
                "termination_reason": 0,
                "potential_energy_saddle": saddle_energy(state, proc),
                "potential_energy_reactant": -1.0,
                "barrier_reactant_to_product": state.get_process(proc)["barrier"],
                "total_force_calls": 1,
            },
            "saddle.con": _buf(saddle),
            "mode.dat": _mode(reactant),
        }
    ]
    before = state.get_process(proc)["repeats"]
    registered = explorer.register_results()
    assert registered == 1
    assert state.get_process(proc)["repeats"] > before


def saddle_energy(state, proc):
    return state.get_process(proc)["saddle_energy"]


def test_catalog_insert_survives_a_raised_store(tmp_path, monkeypatch):
    cfg, _states_obj, state, _product, proc, reactant = _states(tmp_path)
    cfg.kdb_on = True

    def _boom(*_args, **_kwargs):
        raise RuntimeError("catalog down")

    monkeypatch.setattr("eon.process_catalog.insert", _boom)
    fresh = _result(reactant, _atoms(2.2), _atoms(3.6))
    fresh["wuid"] = 8
    assert state.add_process(fresh) is not None


def test_displacement_benchmark_reads_a_configuration(tmp_path, monkeypatch, capsys):
    from eon import displace

    root = tmp_path / "bench"
    root.mkdir()
    cfg = _config(root)
    con = root / "reactant.con"
    io.savecon(str(con), _atoms(2.5))
    out = root / "out.con"
    monkeypatch.setattr(
        sys,
        "argv",
        ["displace", str(con), str(out), cfg.config_path],
    )
    runpy.run_path(displace.__file__, run_name="__main__")
    assert out.is_file()
    assert "displacements per second" in capsys.readouterr().out
