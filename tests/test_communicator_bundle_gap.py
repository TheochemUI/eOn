"""A bundled results index past the file count is still harvested."""

from __future__ import annotations

from types import SimpleNamespace

from eon.communicator import Communicator


def test_bundled_results_with_a_gap_are_kept(tmp_path):
    job = tmp_path / "0_10"
    job.mkdir()
    (job / "results_0.dat").write_text("slot zero\n")
    (job / "results_2.dat").write_text("slot two\n")
    comm = Communicator(tmp_path, 1, config=SimpleNamespace(path_pot=""))
    kept = list(comm.unbundle(str(tmp_path), lambda name: True))
    assert len(kept) == 1
    slots = kept[0]
    assert [slot["name"] for slot in slots] == ["0_10", "0_12"]
    texts = []
    for slot in slots:
        buf = slot["results.dat"]
        buf.seek(0)
        texts.append(buf.read())
    assert texts == ["slot zero\n", "slot two\n"]
