"""Explorer reads and writes under the configured run directories."""

import pickle
from types import SimpleNamespace

from eon.explorer import ClientMinModeExplorer, MinModeExplorer, ServerMinModeExplorer
from eon import fileio as io


def _client(root, results):
    explorer = ClientMinModeExplorer.__new__(ClientMinModeExplorer)
    explorer.config = SimpleNamespace(
        comm_job_buffer_size=1,
        comm_type="local",
        process_search_known_saddles=False,
        saddle_method="min_mode",
        recycling_on=False,
        path_root=str(root),
        path_results=str(results),
        config_path=str(root / "config.ini"),
    )
    explorer.comm = SimpleNamespace(
        queued_search_count=lambda: 0,
        submit_jobs=lambda searches, invariants: None,
    )
    explorer.reactant = object()
    explorer.state = SimpleNamespace(
        number=0,
        info=SimpleNamespace(get=lambda *args, **kwargs: ""),
    )
    explorer.job_table = SimpleNamespace(add_row=lambda row: None, write=lambda: None)
    explorer.wuid = 0
    explorer.save_wuid = lambda: None
    explorer.generate_displacement = lambda: (None, None, "dynamics")
    explorer.nrecycled = 0
    return explorer


def test_make_jobs_stats_pos_con_under_path_root(tmp_path, monkeypatch):
    root = tmp_path / "run dir"
    other = tmp_path / "other dir"
    root.mkdir()
    other.mkdir()
    (root / "pos.con").write_text("frame\n")
    (root / "config.ini").write_text("[Main]\njob = process_search\n")
    monkeypatch.chdir(other)
    monkeypatch.setattr(io, "savecon", lambda *args, **kwargs: None)
    monkeypatch.setattr(io, "modify_config", lambda *args, **kwargs: io.StringIO("x"))
    _client(root, tmp_path / "results").make_jobs()


def test_server_explorer_loads_pickle_from_path_root(tmp_path, monkeypatch):
    root = tmp_path / "run dir"
    other = tmp_path / "other dir"
    root.mkdir()
    other.mkdir()
    with (root / "explorer.pickle").open("wb") as handle:
        pickle.dump({"search_id": 7}, handle)
    monkeypatch.chdir(other)
    monkeypatch.setattr(MinModeExplorer, "__init__", lambda *args, **kwargs: None)
    cfg = SimpleNamespace(path_root=str(root), path_results=str(tmp_path / "results"))
    explorer = ServerMinModeExplorer(None, None, None, config=cfg)
    assert explorer.search_id == 7


def test_server_explorer_saves_under_the_run_directories(tmp_path, monkeypatch):
    root = tmp_path / "run dir"
    results = tmp_path / "results dir"
    other = tmp_path / "other dir"
    root.mkdir()
    results.mkdir()
    other.mkdir()
    monkeypatch.chdir(other)
    explorer = ServerMinModeExplorer.__new__(ServerMinModeExplorer)
    explorer.config = SimpleNamespace(path_root=str(root), path_results=str(results))
    explorer.states = object()
    explorer.previous_state = object()
    explorer.state = object()
    explorer.comm = object()
    explorer.job_info = {
        0: {"0_1": {"type": "saddle_search", "status": "running"}},
    }
    explorer.save()
    assert (root / "explorer.pickle").is_file()
    assert (results / "searches.log").is_file()
    assert not (other / "explorer.pickle").exists()
    assert not (other / "searches.log").exists()
