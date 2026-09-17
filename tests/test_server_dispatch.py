"""Unit tests for pure CPython server entry / job dispatch."""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_select_job_runner_known_jobs():
    from eon.server import select_job_runner

    assert select_job_runner("akmc") is not None
    assert select_job_runner("AKMC") is not None
    assert select_job_runner("basin_hopping") is not None
    assert select_job_runner("parallel_replica") is not None
    assert select_job_runner("escape_rate") is not None
    assert select_job_runner("not_a_job") is None


def test_server_module_exports_main_and_server():
    srv = importlib.import_module("eon.server")

    assert callable(srv.server)
    assert callable(srv.main)
    assert srv.main.__module__ == "eon.server"


def test_pyproject_console_script_points_at_server():
    text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    # Active console_scripts entry (ignore historical mentions in comments).
    assert 'eon-server = "eon.server:main"' in text
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("eon-server"):
            assert "eon.server:main" in stripped
            assert "eon.main:main" not in stripped


def test_import_eon_does_not_load_server():
    env = dict(**__import__("os").environ)
    env["PYTHONPATH"] = str(REPO) + (
        ":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys, eon; assert 'eon.server' not in sys.modules; "
            "_ = eon.server; assert 'eon.server' in sys.modules",
        ],
        cwd=str(REPO),
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_import_eon_server_from_source_tree():
    # Real import path used by CPython package surface (no eonclient required).
    mod = importlib.import_module("eon.server")
    assert hasattr(mod, "server")
    assert hasattr(mod, "main")


def test_python_m_eon_server_help_exits_clean(tmp_path, monkeypatch):
    """``python -m eon.server`` must not raise ImportError on entry."""
    # Avoid config.init reading cwd junk: invoke module load only.
    env = dict(**__import__("os").environ)
    env["PYTHONPATH"] = str(REPO) + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    # Compile/import check: run a one-liner that imports the entry module
    proc = subprocess.run(
        [sys.executable, "-c", "from eon.server import main, server; assert callable(main)"],
        cwd=str(REPO),
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "ModuleNotFoundError" not in proc.stderr
    assert "eon.main" not in proc.stderr


def test_get_communicator_local_interface(tmp_path, monkeypatch):
    """Local communicator is selected via real get_communicator()."""
    from eon import communicator
    from eon.config import ConfigClass

    communicator.reset_communicators()

    # Dummy absolute client path so Local accepts the binary location.
    client = tmp_path / "eonclient"
    client.write_text("#!/bin/sh\n")
    client.chmod(0o755)

    cfg = ConfigClass()
    cfg.comm_type = "local"
    cfg.path_scratch = str(tmp_path / "scratch")
    cfg.comm_local_client = str(client)
    cfg.comm_local_ncpus = 1
    cfg.comm_job_bundle_size = 1

    comm = communicator.get_communicator(cfg)
    assert comm is not None
    assert hasattr(comm, "submit_jobs")
    assert hasattr(comm, "get_results")
    assert type(comm).__name__ == "Local"
