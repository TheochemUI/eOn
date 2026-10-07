"""Atoms that share a Hessian color are farther apart than the cutoff."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_shared_colors_are_farther_than_the_cutoff():
    text = (ROOT / "client" / "unit_tests" / "HessianTest.cpp").read_text(
        encoding="utf-8"
    )
    start = text.find("Colored FD Hessian matches serial central difference")
    assert start != -1
    window = text[start : start + 1800]
    assert "REQUIRE(matter->distance(i, j) > kCutoff);" in window
    assert "REQUIRE(shared);" in window
    assert "kCutoff = 1.5" in window
