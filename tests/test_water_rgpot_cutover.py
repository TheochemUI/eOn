"""The water factory arms are the rgpot kernels."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_water_family_uses_the_rgpot_adapters():
    source = (ROOT / "client" / "Potential.cpp").read_text(encoding="utf-8")
    assert "rgpot::TIP4PPot" in source
    assert "rgpot::SPCEPot" in source
    assert "rgpot::TIP4PPtPot" in source
    assert "rgpot::fortranpots::WaterHPot" in source
    assert "std::make_unique<Tip4p>" not in source
    assert "std::make_unique<SpceCcl>" not in source
    assert "std::make_unique<Tip4p_Pt>" not in source
    assert "potentials/Water/Water.hpp" not in source
    assert not (ROOT / "client" / "potentials" / "Water" / "tip4p_ccl.cpp").exists()
    assert not (ROOT / "client" / "potentials" / "Water_Pt" / "Tip4p_Pt.cpp").exists()
    meson = (ROOT / "client" / "meson.build").read_text(encoding="utf-8")
    assert "'with_water': get_option('with_water')" in meson
    assert "subdir('potentials/Water')" not in meson
    wrap = (ROOT / "subprojects" / "rgpot.wrap").read_text(encoding="utf-8")
    assert "e40fd175f5230098d56333bf2aff38cd839dc810" in wrap
