"""A rate instanton writes its steepest-descent path for initial_path."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_rate_job_writes_the_steepest_descent_path():
    source = (ROOT / "client" / "InstantonJob.cpp").read_text(encoding="utf-8")
    assert 'writeConFrames("instanton_sd_path.con"' in source
    assert '{"arc_length", arc[k]}' in source
    assert 'extras.emplace_back("sd_path_force_calls"' in source
    case = (ROOT / "client" / "unit_tests" / "JobIntegrationTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "InstantonJob reuses the LJ13 steepest-descent path" in case
    assert 'REQUIRE(std::stod(second.at("sd_path_force_calls")) == 0.0);' in case
    assert "1e-6" in case
    guide = (ROOT / "docs" / "source" / "user_guide" / "instanton.md").read_text(
        encoding="utf-8"
    )
    assert "instanton_sd_path.con" in guide
