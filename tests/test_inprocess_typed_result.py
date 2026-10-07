"""The in-process path keeps a typed job result and no results.dat buffer."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_submit_does_not_require_a_results_dat_buffer():
    text = (ROOT / "eon" / "communicator_inprocess.py").read_text(encoding="utf-8")
    start = text.find("job_result = _job_result(")
    assert start != -1
    window = text[start : start + 700]
    assert '"job_result": job_result' in window
    assert '"results.dat"' not in window
    assert "StringIO" not in window
    assert "_results_dat(" not in window
