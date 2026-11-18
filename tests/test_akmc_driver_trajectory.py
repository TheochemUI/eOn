"""A whole AKMC trajectory through the server driver, on the Cu fixture.

The recorded trajectory was produced by the driver before it was made
faster: dynamics rows, state table, process tables and superbasin files,
as text. A driver that takes a different decision anywhere along 150
cycles changes one of them. Alongside, the run is counted: a state's info
file is parsed once when its ``ini`` is made and once more after each
``set``, not once per ``get``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

pytest.importorskip("linkcell")
pytest.importorskip("minimage")

import akmc_synthetic as synth  # noqa: E402
from driver_oracles import count_events  # noqa: E402

GOLDEN = _TESTS / "data" / "akmc_synthetic_trajectory.json"
CELLS, CYCLES, SEARCHES, TRANSITIONS = 3, 150, 12, 2


@pytest.fixture(scope="module")
def trajectory(tmp_path_factory):
    root = tmp_path_factory.mktemp("trajectory")
    cwd = Path.cwd()
    try:
        with count_events() as counts:
            cfg, _ = synth.run(root, CELLS, CYCLES, SEARCHES, TRANSITIONS)
        return synth.snapshot(cfg), counts
    finally:
        import os

        os.chdir(cwd)


def test_the_trajectory_is_the_recorded_one(trajectory):
    snapshot, _ = trajectory
    golden = json.loads(GOLDEN.read_text())
    assert snapshot["dynamics"] == golden["dynamics"]
    assert snapshot["state_table"] == golden["state_table"]
    assert snapshot["processtables"] == golden["processtables"]
    assert snapshot["superbasins"] == golden["superbasins"]


def test_the_recorded_trajectory_exercises_basins_and_revisits():
    golden = json.loads(GOLDEN.read_text())
    steps = [row.split() for row in golden["dynamics"][2:]]
    reactants = [int(f[1]) for f in steps]
    assert len(set(reactants)) < len(reactants)
    assert golden["superbasins"]
    assert len(golden["processtables"]) > 10


def test_an_info_file_is_parsed_once_per_write_not_once_per_read(trajectory):
    _, counts = trajectory
    parses = counts["py:RawConfigParser._read"]
    made = counts["py:ini.__init__"]
    writes = counts["py:ini.set"]
    assert made > 0 and writes > 0
    assert parses <= made + writes, (parses, made, writes)
