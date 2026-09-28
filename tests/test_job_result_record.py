"""The in-process result is a dict. results.dat is derived from it."""

from eon.communicator_inprocess import _job_result, _results_dat


def test_results_dat_is_the_typed_record():
    record = _job_result(0, -1.25, 7, "minimization")
    text = _results_dat(0, -1.25, 7, "minimization")
    assert record["termination_reason"] == 0
    assert record["job_type"] == "minimization"
    assert "minimization job_type" in text
    assert "7 total_force_calls" in text
    assert "-1.250000000000e+00 potential_energy" in text
