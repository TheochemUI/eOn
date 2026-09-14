from eon_schema.jobs import (
    dict_to_results_dat,
    job_result_capnp_path,
    job_result_scalars_from_results_dat,
    results_dat_to_dict,
)


def test_job_result_capnp_exists():
    path = job_result_capnp_path()
    assert path.is_file()
    text = path.read_text(encoding="utf-8")
    assert "struct JobResult" in text
    assert "struct Geometry" in text
    assert "struct JobRequest" in text
    assert "positions" in text
    assert "statusCode" in text


def test_results_dat_roundtrip_scalars():
    sample = (
        "0 termination_reason\n"
        "good termination_reason_text\n"
        "process_search job_type\n"
        "lj potential_type\n"
        "42 total_force_calls\n"
        "1.234567890123e+00 potential_energy_saddle\n"
        "1.000000000000e+00 potential_energy_reactant\n"
        "1.100000000000e+00 potential_energy_product\n"
        "2.345678901234e-01 barrier_reactant_to_product\n"
    )
    d = results_dat_to_dict(sample)
    assert d["termination_reason"] == 0
    assert d["job_type"] == "process_search"
    assert abs(d["potential_energy_saddle"] - 1.234567890123) < 1e-12
    mapped = job_result_scalars_from_results_dat(sample)
    assert mapped["status_code"] == 0
    assert mapped["force_calls"]["total"] == 42
    assert abs(mapped["barrier_reactant_to_product"] - 0.2345678901234) < 1e-10
    # round-trip text
    again = dict_to_results_dat(
        {
            "termination_reason": 0,
            "termination_reason_text": "good",
            "total_force_calls": 42,
        }
    )
    assert "0 termination_reason" in again
    assert "42 total_force_calls" in again
