"""Server-side process search errors must exit non-zero."""

from __future__ import annotations

import pytest

from eon.config import ConfigClass


def _ini(tmp_path, name, prefactor, bundle):
    path = tmp_path / name
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "[AKMC]",
                "server_side_process_search = true",
                "[Prefactor]",
                "default_value = %s" % prefactor,
                "[Communicator]",
                "jobs_per_bundle = %s" % bundle,
                "[Paths]",
                "main_directory = %s" % tmp_path,
                "results = %s" % tmp_path,
                "states = %s" % (tmp_path / "states"),
                "",
            ]
        )
    )
    return path


def test_missing_prefactor_exits_nonzero(tmp_path, capsys):
    path = _ini(tmp_path, "server.ini", "0", "1")
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(path))
    assert caught.value.code == 1
    captured = capsys.readouterr()
    assert "default prefactor" in captured.err
    assert "default prefactor" not in captured.out


def test_bundle_size_exits_nonzero(tmp_path, capsys):
    path = _ini(tmp_path, "bundle.ini", "1e12", "2")
    with pytest.raises(SystemExit) as caught:
        ConfigClass().init(str(path))
    assert caught.value.code == 1
    captured = capsys.readouterr()
    assert "bundle size" in captured.err
    assert "bundle size" not in captured.out
