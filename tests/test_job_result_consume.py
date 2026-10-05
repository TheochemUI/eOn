"""In-process records expose JobResult; results.dat stays a file adapter."""

from eon.job_result import results_mapping
from eon_schema.jobs import job_result_to_results_dat, results_dat_to_dict


def test_consume_job_result_without_results_dat():
    record = {
        "job_result": {
            "job_type": "process_search",
            "status_code": 0,
            "status_text": "good",
            "potential_energy_saddle": 1.5,
            "potential_energy_reactant": 1.0,
            "potential_energy_product": 1.2,
            "barrier_reactant_to_product": 0.5,
            "force_calls": {
                "total": 9,
                "minimization": 4,
                "saddle": 5,
            },
            "saddle": object(),
            "product": object(),
        }
    }
    mapped = results_mapping(record)
    assert mapped["termination_reason"] == 0
    assert mapped["job_type"] == "process_search"
    assert mapped["total_force_calls"] == 9
    assert mapped["force_calls_minimization"] == 4
    assert mapped["potential_energy_saddle"] == 1.5
    assert "saddle" not in mapped
    text = job_result_to_results_dat(record["job_result"])
    again = results_dat_to_dict(text)
    assert again["termination_reason"] == 0
    assert again["barrier_reactant_to_product"] == 0.5
    assert "positions" not in text
