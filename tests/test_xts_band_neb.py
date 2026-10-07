"""The default NEB step is the xtsci band session."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_neb_lbfgs_uses_xts_band_step():
    band = (ROOT / "client" / "XtsciBand.cpp").read_text(encoding="utf-8")
    neb = (ROOT / "client" / "NudgedElasticBand.cpp").read_text(encoding="utf-8")
    assert "xts_band_step(" in band
    assert "bandMethod == OptType::LBFGS" in neb
    assert "bandMethod == OptType::FIRE" in neb
