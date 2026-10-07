"""The in-process path keeps a typed job result beside results.dat."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_submit_keeps_the_typed_result_beside_results_dat():
    text = (ROOT / "eon" / "communicator_inprocess.py").read_text(encoding="utf-8")
    start = text.find("job_result = _job_result(")
    assert start != -1
    window = text[start : start + 900]
    assert '"job_result": job_result' in window
    assert '"results.dat": results' in window
    assert "_results_dat(" in window
    helper = text[text.find("def _results_dat(") : text.find("class LocalInProcess")]
    assert "_job_result(" in helper
