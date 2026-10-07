"""A solid-state band uses only the stress the potential reports."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_solid_state_band_refuses_a_missing_stress():
    source = (ROOT / "client" / "NudgedElasticBand.cpp").read_text(
        encoding="utf-8"
    )
    start = source.find("void NudgedElasticBand::projectSolidState")
    assert start != -1
    end = source.find("\nvoid ", start + 1)
    body = source[start:end]
    assert "finiteDifferenceCauchyStresses" not in body
    assert (
        "solid_state NEB requires a potential that reports the Cauchy stress"
        in body
    )
    case = (
        ROOT / "client" / "unit_tests" / "SolidStateNEBTest.cpp"
    ).read_text(encoding="utf-8")
    assert "solid_state refuses a potential that does not report stress" in case
    assert (
        "REQUIRE_THROWS_WITH(\n      neb->updateForces(),\n"
        '      "solid_state NEB requires a potential that reports the Cauchy stress");'
        in case
    )
    guide = (ROOT / "docs" / "source" / "user_guide" / "neb.md").read_text(
        encoding="utf-8"
    )
    assert "does not report\nthe Cauchy stress is refused" in guide
