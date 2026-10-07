"""Above the crossover the rate job names the parabolic correction."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_the_rate_job_does_not_evaluate_the_parabolic_factor():
    job = (ROOT / "client" / "InstantonJob.cpp").read_text(encoding="utf-8")
    assert "parabolicFactor(" not in job
    assert "the parabolic barrier correction is not evaluated" in job
    case = (ROOT / "client" / "unit_tests" / "JobIntegrationTest.cpp").read_text(
        encoding="utf-8"
    )
    assert 'REQUIRE(results.at("termination_reason") == "2");' in case
    assert 'REQUIRE(results.count("parabolic_factor") == 0);' in case
    guide = (ROOT / "docs" / "source" / "user_guide" / "instanton.md").read_text(
        encoding="utf-8"
    )
    assert "does not evaluate it" in guide
