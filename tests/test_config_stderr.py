"""A missing configuration file is reported on standard error."""

from __future__ import annotations

import pytest

from eon.config import ConfigClass


def test_missing_file_is_reported_on_stderr(tmp_path, capsys):
    missing = tmp_path / "absent.ini"
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(missing))
    assert caught.value.code == 2
    captured = capsys.readouterr()
    assert "does not exist" in captured.err
    assert "does not exist" not in captured.out


def test_missing_config_ini_is_reported_on_stderr(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init("")
    assert caught.value.code == 2
    captured = capsys.readouterr()
    assert "must provide a configuration file" in captured.err
    assert "must provide a configuration file" not in captured.out
