"""Finite-difference Hessians color the squared cutoff graph."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_fd_hessian_colors_the_cutoff_graph():
    header = (ROOT / "include" / "eon" / "Hessian.h").read_text(encoding="utf-8")
    assert "colorMobileCutoffGraph" in header
    assert "greedyColorCutoffGraph" in header
    source = (ROOT / "client" / "Hessian.cpp").read_text(encoding="utf-8")
    assert "Square of the cutoff graph" in source
    assert "calculateColored" in source
    case = (ROOT / "client" / "unit_tests" / "HessianTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "Colored FD Hessian matches serial central difference" in case
    assert "REQUIRE(nColors == 3);" in case
    assert "REQUIRE(calls == static_cast<size_t>(1 + 2 * 3 * nColors));" in case
    assert "WithinAbs(serial(i, j), 1e-8)" in case
    guide = (ROOT / "docs" / "source" / "user_guide" / "hessian.md").read_text(
        encoding="utf-8"
    )
    assert "three colors" in guide
    assert "19" in guide
