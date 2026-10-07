"""Dimer and Lanczos rotation go through the xtsci min-mode session."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dimer_and_lanczos_use_xts_minmode_estimate():
    factory = (ROOT / "client" / "EigenmodeStrategy.cpp").read_text(encoding="utf-8")
    mode = (ROOT / "client" / "XtsciMinMode.cpp").read_text(encoding="utf-8")
    assert factory.count("make_shared<XtsciMinMode>") >= 3
    assert "xts_minmode_estimate(" in mode
