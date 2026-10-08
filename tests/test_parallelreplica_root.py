"""A missing replica root raises, and the logged speedup is the mean."""

from types import SimpleNamespace

import pytest

from eon.parallelreplica import parallelreplica


def test_missing_root_raises(tmp_path):
    config = SimpleNamespace(path_root=str(tmp_path / "absent"))
    with pytest.raises(FileNotFoundError, match="root directory does not exist"):
        parallelreplica(config)


def test_registered_speedup_is_logged_as_the_mean(tmp_path, monkeypatch, caplog):
    root = tmp_path / "root"
    root.mkdir()
    results = tmp_path / "results"
    results.mkdir()
    state = SimpleNamespace(number=0, get_time=lambda: 0.0)
    states = SimpleNamespace(get_state=lambda number: state)
    monkeypatch.setattr(
        "eon.parallelreplica.get_pr_metadata", lambda config: (0, 0.0, 1)
    )
    monkeypatch.setattr(
        "eon.parallelreplica.get_statelist", lambda config: states
    )
    monkeypatch.setattr(
        "eon.parallelreplica.communicator.get_communicator", lambda config: object()
    )
    monkeypatch.setattr(
        "eon.parallelreplica.register_results",
        lambda *args, **kwargs: (2, None, 8.0),
    )
    monkeypatch.setattr(
        "eon.parallelreplica.make_searches", lambda *args, **kwargs: 3
    )
    monkeypatch.setattr(
        "eon.parallelreplica.io.write_info_txt", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        "eon.parallelreplica.io.save_prng_state", lambda *args, **kwargs: None
    )
    config = SimpleNamespace(path_root=str(root), path_results=str(results))
    caplog.set_level("INFO", logger="pr")
    parallelreplica(config)
    assert any(message == "Average speedup: 4.000000" for message in caplog.messages)
