"""Golden masters for CPython server packaging / job-dispatch surfaces."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

# Job keys registered by the shipped select_job_runner (see eon/server.py).
# This is the public dispatch surface from #368, not a reimplementation.
EXPECTED_JOB_KEYS = {
    "akmc",
    "parallel_replica",
    "unbiased_parallel_replica",
    "basin_hopping",
    "escape_rate",
}


def test_eon_exports_version():
    import eon
    from eon import version

    assert version == eon.__version__
    assert isinstance(version, str)
    assert version


def test_select_job_runner_golden_keys():
    from eon.server import select_job_runner

    for key in EXPECTED_JOB_KEYS:
        assert select_job_runner(key) is not None, key
        assert callable(select_job_runner(key))
    # Unknown → None (fallback single-job path)
    assert select_job_runner("not_registered_xyz") is None
    # Case fold
    assert select_job_runner("AKMC") is not None


def test_console_script_and_main_chain():
    text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    assert 'eon-server = "eon.server:main"' in text
    from eon.server import main, server

    assert main.__module__ == "eon.server"
    assert server.__module__ == "eon.server"
    # main is the thin console wrapper
    assert main.__doc__ is not None
