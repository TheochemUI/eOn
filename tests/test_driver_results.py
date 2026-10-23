"""The results directory exists before a driver opens its log or lock."""

import logging
import os
import sys
from pathlib import Path

import pytest

from eon import fileio as io
from eon.akmc import main as akmc_main
from eon.basinhopping import main as hop_main
from eon.escaperate import main as escape_main
from eon.parallelreplica import main as replica_main
from tests.test_library_bodies import _atoms


def _reset_logging():
    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
        handler.close()


def _prepare(root: Path, job: str, client: Path, results: Path):
    root.mkdir()
    io.savecon(str(root / "pos.con"), _atoms(2.5))
    ini = root / "config.ini"
    ini.write_text(
        "\n".join(
            [
                "[Main]",
                f"job = {job}",
                "temperature = 300",
                "[Paths]",
                f"main_directory = {root}",
                f"results = {results}",
                f"states = {root / 'states'}",
                "[Communicator]",
                "type = local",
                f"client_path = {client}",
                "jobs_per_bundle = 1",
                "num_jobs = 1",
                "",
            ]
        )
    )
    return ini


@pytest.mark.parametrize(
    ("job", "entry", "log_name"),
    [
        ("basin_hopping", hop_main, "bh.log"),
        ("parallel_replica", replica_main, "pr.log"),
        ("escape_rate", escape_main, "pr.log"),
        ("akmc", akmc_main, "akmc.log"),
    ],
)
def test_missing_results_directory_is_created(tmp_path, monkeypatch, job, entry, log_name):
    client = tmp_path / "client.py"
    client.write_text("#!/usr/bin/env python3\n")
    client.chmod(0o755)
    _reset_logging()
    root = tmp_path / job
    results = root / "output"
    ini = _prepare(root, job, client, results)
    assert not results.exists()
    monkeypatch.chdir(root)
    monkeypatch.setattr(sys, "argv", [job, str(ini), "-q", "-n"])
    entry()
    assert results.is_dir()
    lock = results / "lockfile"
    assert lock.is_file()
    pid = int(lock.read_text().strip())
    assert pid == os.getpid()
    assert (results / log_name).is_file()
