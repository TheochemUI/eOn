"""Above the crossover the rate job evaluates the parabolic factor."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_the_rate_job_evaluates_the_parabolic_factor():
    job = (ROOT / "client" / "InstantonJob.cpp").read_text(encoding="utf-8")
    assert "tunneling::parabolicFactor(temperature, tc)" in job
    assert 'extras.emplace_back("parabolic_factor", factor);' in job
    assert 'extras.emplace_back("rate_parabolic", kPar);' in job
    case = (ROOT / "client" / "unit_tests" / "JobIntegrationTest.cpp").read_text(
        encoding="utf-8"
    )
    assert 'REQUIRE(std::stod(results.at("parabolic_factor")) > 1.0);' in case
    guide = (ROOT / "docs" / "source" / "user_guide" / "instanton.md").read_text(
        encoding="utf-8"
    )
    assert "writes the parabolic barrier factor" in guide
