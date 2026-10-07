"""Cancel is polled inside force, relax, NEB, saddle search and dynamics."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_compiled_loops_poll_the_cancel_token():
    matter = (ROOT / "client" / "Matter.cpp").read_text(encoding="utf-8")
    relax = (ROOT / "client" / "HelperFunctions.cpp").read_text(encoding="utf-8")
    neb = (ROOT / "client" / "NudgedElasticBand.cpp").read_text(encoding="utf-8")
    saddle = (ROOT / "client" / "MinModeSaddleSearch.cpp").read_text(encoding="utf-8")
    dynamics = (ROOT / "client" / "Dynamics.cpp").read_text(encoding="utf-8")
    assert 'cancel_token_.poll("force")' in matter
    assert 'pollCancel("relax")' in relax
    assert 'pollCancel("neb")' in neb
    assert 'pollCancel("saddle")' in saddle
    assert 'pollCancel("dynamics")' in dynamics
    case = (ROOT / "client" / "unit_tests" / "MatterTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "REQUIRE_THROWS_AS(matter.getPotentialEnergy(), JobCancelled)" in case
    assert "REQUIRE_THROWS_AS(eonc::helpers::relaxMatter(matter, params, true)" in case
