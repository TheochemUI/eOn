"""The EMT factory arm is the rgpot kernel."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_emt_is_the_rgpot_adapter():
    source = (ROOT / "client" / "Potential.cpp").read_text(encoding="utf-8")
    assert "rgpot::EMTPot" in source
    assert "EMTConfig" in source
    assert "EffectiveMediumTheory" not in source
    assert "potentials/EAM/EAM.h" not in source
    pin = (ROOT / "client" / "unit_tests" / "EMTCuTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "5.129167" in pin
    assert "1.914263" in pin
    wrap = (ROOT / "subprojects" / "rgpot.wrap").read_text(encoding="utf-8")
    assert "d475f890dd3e3ee6b58643da0d2b9c92608b6476" in wrap
