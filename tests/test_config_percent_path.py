"""A path may contain a percent sign, and %(name)s still expands."""

from __future__ import annotations

from eon.config import ConfigClass


def _load(tmp_path, main_directory):
    path = tmp_path / "config.ini"
    path.write_text(
        "[Main]\njob = akmc\n[Paths]\nmain_directory = %s\n" % main_directory
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def test_percent_in_main_directory(tmp_path):
    cfg = _load(tmp_path, "/tmp/a%b")
    assert cfg.path_root == "/tmp/a%b"
    assert cfg.path_jobs_out == "/tmp/a%b/jobs/out/"


def test_doubled_percent_is_one_percent(tmp_path):
    cfg = _load(tmp_path, "/tmp/a%%b")
    assert cfg.path_root == "/tmp/a%b"
    assert cfg.path_jobs_out == "/tmp/a%b/jobs/out/"
