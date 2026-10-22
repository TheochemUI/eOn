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
    assert "e40fd175f5230098d56333bf2aff38cd839dc810" in wrap
    assert "rgpot-eam-al-r2549-literals.patch" in wrap
    patch = (
        ROOT / "subprojects" / "packagefiles" / "rgpot-eam-al-r2549-literals.patch"
    ).read_text(encoding="utf-8")
    assert "real(7.1, real32)" in patch
    assert "9e1abb3cd3eeb5cd276c3ad20a809d46fc395ec6" in wrap
