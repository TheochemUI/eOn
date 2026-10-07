"""The rate instanton past the Newton limit is a dimer search on one ring."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _rate_body(source: str) -> str:
    key = "RateInstanton optimizeRateInstanton("
    start = source.index(key)
    rest = source[start:]
    end = rest.index("\nvoid instantonRate")
    return rest[:end]


def test_rate_instanton_calls_the_shared_dimer():
    tunneling = (ROOT / "client" / "Tunneling.cpp").read_text(encoding="utf-8")
    body = _rate_body(tunneling)
    assert "searchRingWithDimer(" in body
    assert "lowestMode(" not in body
    assert "lowestMode(" not in tunneling
    assert "MINMODE_DIMER" in tunneling
    header = (ROOT / "include" / "eon" / "RingPolymerPotential.h").read_text(
        encoding="utf-8"
    )
    assert "class RingPolymerPotential" in header
    assert "forceBatch(" in header
    guide = (ROOT / "docs" / "source" / "user_guide" / "instanton.md").read_text(
        encoding="utf-8"
    )
    assert "one force batch" in guide
