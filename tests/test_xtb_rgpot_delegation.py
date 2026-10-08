"""XTB is the rgpot kernel, and the engine libraries use the rgpot soname."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_xtb_delegates_and_the_engines_use_the_rgpot_soname():
    source = (ROOT / "client" / "Potential.cpp").read_text(encoding="utf-8")
    assert "makeRgpot<rgpot::XTBPot>" in source
    assert "rgpot::XTBConfig" in source
    assert ".accuracy = o.acc" in source
    assert ".electronic_temperature = o.elec_temperature" in source
    assert ".max_iterations = static_cast<int>(o.maxiter)" in source
    assert ".charge = o.charge" in source
    assert ".uhf = o.uhf" in source
    assert not (ROOT / "client" / "potentials" / "XTBPot" / "XTBPot.cpp").exists()
    assert not (
        ROOT / "include" / "eon" / "potentials" / "XTBPot" / "XTBPot.h"
    ).exists()
    adapter = (
        ROOT / "include" / "eon" / "potentials" / "RgpotAdapter" / "RgpotAdapter.h"
    ).read_text(encoding="utf-8")
    assert "pot_.forceBatchImpl(" in adapter
    xtb_loader = (
        ROOT / "client" / "potentials" / "Rgpot" / "XTBEngineLoader.cpp"
    ).read_text(encoding="utf-8")
    assert "librgpot_xtb_engine.so" in xtb_loader
    mta_loader = (
        ROOT / "client" / "potentials" / "Rgpot" / "MetatomicEngineLoader.cpp"
    ).read_text(encoding="utf-8")
    assert "librgpot_metatomic_engine.so" in mta_loader
    meson = (ROOT / "client" / "meson.build").read_text(encoding="utf-8")
    assert "libmetatomic_engine.so" in meson
