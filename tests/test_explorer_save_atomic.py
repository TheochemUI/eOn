"""A killed explorer save leaves the previous pickle and log in place."""

import os
import pickle
from types import SimpleNamespace

import pytest

from eon.explorer import ServerMinModeExplorer


def _explorer(root, results):
    explorer = ServerMinModeExplorer.__new__(ServerMinModeExplorer)
    explorer.config = SimpleNamespace(path_root=str(root), path_results=str(results))
    explorer.states = object()
    explorer.previous_state = object()
    explorer.state = object()
    explorer.comm = object()
    explorer.search_id = 3
    explorer.wuid_to_search_id = {}
    explorer.process_searches = {}
    explorer.job_info = {}
    return explorer


def _run_dirs(tmp_path, monkeypatch):
    root = tmp_path / "run dir"
    results = tmp_path / "results dir"
    other = tmp_path / "other dir"
    root.mkdir()
    results.mkdir()
    other.mkdir()
    monkeypatch.chdir(other)
    return root, results


def test_killed_pickle_rewrite_keeps_the_previous_explorer(tmp_path, monkeypatch):
    root, results = _run_dirs(tmp_path, monkeypatch)
    previous = {"search_id": 9}
    (root / "explorer.pickle").write_bytes(pickle.dumps(previous))
    (results / "searches.log").write_text("OLD\n")

    def killed(_src, _dst):
        raise OSError("killed")

    monkeypatch.setattr("eon.fileio.os.replace", killed)
    with pytest.raises(OSError, match="killed"):
        _explorer(root, results).save()
    with (root / "explorer.pickle").open("rb") as handle:
        assert pickle.load(handle) == previous
    assert (results / "searches.log").read_text() == "OLD\n"


def test_killed_log_rewrite_keeps_the_previous_log(tmp_path, monkeypatch):
    root, results = _run_dirs(tmp_path, monkeypatch)
    (root / "explorer.pickle").write_bytes(pickle.dumps({"search_id": 9}))
    (results / "searches.log").write_text("OLD\n")
    real_replace = os.replace
    calls = {"n": 0}

    def fail_the_second(src, dst):
        calls["n"] += 1
        if calls["n"] >= 2:
            raise OSError("killed")
        return real_replace(src, dst)

    monkeypatch.setattr("eon.fileio.os.replace", fail_the_second)
    with pytest.raises(OSError, match="killed"):
        _explorer(root, results).save()
    assert (results / "searches.log").read_text() == "OLD\n"
