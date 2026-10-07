"""The rate instanton carries an implicit and an explicit friction bath."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_rate_job_reads_the_friction_bath():
    header = (ROOT / "include" / "eon" / "Tunneling.h").read_text(encoding="utf-8")
    assert "void addFrictionBath" in header
    assert "frictionExplicit" in header
    model = (
        ROOT / "packages" / "eon-schema" / "src" / "eon_schema" / "config" / "models.py"
    ).read_text(encoding="utf-8")
    assert 'Literal["none", "implicit", "explicit"]' in model
    job = (ROOT / "client" / "InstantonJob.cpp").read_text(encoding="utf-8")
    assert "ro.frictionExplicit = o.friction ==" in job
    guide = (ROOT / "eon" / "config.yaml").read_text(encoding="utf-8")
    assert "friction_eta_beads:" in guide
