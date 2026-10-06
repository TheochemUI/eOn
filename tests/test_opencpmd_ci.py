"""The OpenCPMD workflow is a real engine, not the fake library."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "ci_opencpmd_point.yml"
FAKE = ROOT / ".github" / "workflows" / "ci_build_akmc.yml"
PIN = "062582b7cfd832d36f88f504cd08e4ead42eb404"


def test_opencpmd_point_job_is_pinned():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert PIN in text
    assert "with_cpmd" in (ROOT / "scripts" / "ci" / "opencpmd_point.sh").read_text(
        encoding="utf-8"
    )
    assert "libcpmdc.so" in text or "libcpmdc.so" in (
        ROOT / "scripts" / "ci" / "opencpmd_point.sh"
    ).read_text(encoding="utf-8")
    assert "-1396.269526" in (ROOT / "scripts" / "ci" / "opencpmd_point.sh").read_text(
        encoding="utf-8"
    )
    assert "cpmdc_fake_engine" not in text
    # The fast suite still uses the fake library. This job does not replace it.
    assert "cpmdc_fake_engine" in FAKE.read_text(encoding="utf-8")
