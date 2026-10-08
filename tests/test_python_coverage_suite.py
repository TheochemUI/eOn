"""The Python coverage script measures every module under eon/."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_coverage_python.sh"


def test_coverage_script_runs_the_package_suite():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "--cov=eon" in text
    assert "\n  tests \\\n" in text or "\n  tests\n" in text
    assert "eon/tests" in text
    assert "test_config_metadata.py" not in text
    assert "test_job_runners.py" not in text
    assert "scripts/eon-coverage.cfg" in text
    assert 'version = "0.dev-coverage"' in text
    assert "__version__" not in text
    cfg = (ROOT / "scripts" / "eon-coverage.cfg").read_text(encoding="utf-8")
    assert "*/eon/tests/*" in cfg
