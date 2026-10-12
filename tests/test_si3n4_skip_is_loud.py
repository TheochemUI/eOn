"""The Si3N4 point skip has to show up in the job log."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_opencpmd_job_warns_and_skips_with_a_reason():
    workflow = (ROOT / ".github/workflows/ci_opencpmd_point.yml").read_text(
        encoding="utf-8"
    )
    skip_test = (ROOT / "scripts/ci/test_si3n4_point_skip.py").read_text(
        encoding="utf-8"
    )
    reason = (
        "TheochemUI/sige_repro is not readable by this token, "
        "so the Si3N4 point is skipped."
    )
    assert "::warning title=Si3N4 point skipped::" + reason in workflow
    assert "::notice title=Si3N4 point skipped::" not in workflow
    assert "python -m pip install -q 'pytest>=8'" in workflow
    assert "scripts/ci/test_si3n4_point_skip.py" in workflow
    assert "steps.sige.outputs.readable != 'true'" in workflow
    assert reason in skip_test
    assert "pytest.skip" in skip_test
