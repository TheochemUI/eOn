"""``StateList.get_product_state`` reads the state table once per call.

Finding or minting the product state of a process needs the number of
states. Reading ``state_table`` again for every use of that number costs a
pass over every state for every process of every state visited.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_TESTS = Path(__file__).resolve().parent
if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

pytest.importorskip("linkcell")
pytest.importorskip("minimage")

import akmc_synthetic as synth  # noqa: E402
from driver_oracles import count_opens, reads_of  # noqa: E402

from eon import akmc  # noqa: E402


@pytest.fixture()
def walked(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = synth.setup_run(tmp_path, 3, transitions=2)
    client = synth.SyntheticClient(3, 7)
    for _ in range(25):
        synth.cycle(cfg, client, 8)
    return cfg, client


def fresh_states(cfg):
    kT = cfg.main_temperature / 11604.5
    return akmc.get_statelist(kT, config=cfg)


def test_every_product_lookup_reads_the_state_table_at_most_once(walked):
    cfg, client = walked
    states = fresh_states(cfg)
    here = states.get_state(0)
    for result in client.results(here, 16):
        if result["results"]["termination_reason"] == 0:
            here.add_process(result)
    here.load_process_table()
    pending = [pid for pid, p in here.procs.items() if p["product"] == -1]
    assert pending
    minted = 0
    for pid in pending:
        before = states.get_num_states()
        with count_opens() as seen:
            states.get_product_state(0, pid)
        assert reads_of(seen, "state_table") <= 1, pid
        minted += states.get_num_states() - before
    assert minted > 0
    assert len(states.state_table) == states.get_num_states()
