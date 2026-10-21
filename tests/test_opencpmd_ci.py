"""The OpenCPMD workflow is a real engine, not the fake library."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "ci_opencpmd_point.yml"
FAKE = ROOT / ".github" / "workflows" / "ci_build_akmc.yml"
PIN = "062582b7cfd832d36f88f504cd08e4ead42eb404"


def test_opencpmd_point_job_is_pinned():
    text = WORKFLOW.read_text(encoding="utf-8")
    script = (ROOT / "scripts" / "ci" / "opencpmd_point.sh").read_text(encoding="utf-8")
    assert PIN in text
    assert "with_cpmd" in script
    assert "libcpmdc.so" in text or "libcpmdc.so" in script
    assert "-1396.269526" in script
    assert "cpmdc_fake_engine" not in text
    # The fast suite still uses the fake library. This job does not replace it.
    assert "cpmdc_fake_engine" in FAKE.read_text(encoding="utf-8")


def test_si3n4_point_skips_without_the_private_token():
    text = WORKFLOW.read_text(encoding="utf-8")
    # Job-level if cannot read secrets. The gate publishes a boolean and
    # the Si3N4 job is skipped when that boolean is not true.
    assert "needs.sige_repro_token.outputs.present == 'true'" in text
    assert "secrets.SIGE_REPRO_TOKEN != ''" in text
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("if:") and "secrets." in stripped:
            raise AssertionError(stripped)


def test_h2_point_uses_committed_inputs():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "name: H2 point" in text
    assert "examples/opencpmd-h2" in text
    assert "TheochemUI/sige_repro" not in text.split("h2_point:", 1)[1]
    kit = ROOT / "examples" / "opencpmd-h2"
    assert (kit / "config.ini").is_file()
    assert (kit / "h2.con").is_file()
    assert (kit / "ext_pot").is_file()
    assert (kit / "pseudo" / "H_GTH_BLYP").is_file()
    # Pins from the local 100 Ry BLYP scan, not placeholders.
    assert 'ENERGY_EV: "-31.795980465"' in text
    assert 'BOND_REF_ANG: "0.748890"' in text
    assert 'ENERGY_TOL_EV: "0.01"' in text
    assert 'BOND_TOL_ANG: "0.005"' in text
    assert "PENDING" not in text
