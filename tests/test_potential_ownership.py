"""Point and minimization keep an exclusive Potential until a Matter is copied."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_point_and_minimization_keep_a_unique_potential():
    point = (ROOT / "include" / "eon" / "PointJob.h").read_text(encoding="utf-8")
    mini = (ROOT / "include" / "eon" / "MinimizationJob.h").read_text(encoding="utf-8")
    matter = (ROOT / "include" / "eon" / "Matter.h").read_text(encoding="utf-8")
    assert "ExclusivePotential{}" in point
    assert "ExclusivePotential{}" in mini
    assert "holdsExclusivePotential" in matter
    assert "ownedPotential" in matter
    case = (ROOT / "client" / "unit_tests" / "MatterTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "REQUIRE(a.holdsExclusivePotential())" in case
    assert "REQUIRE(std::isfinite(a.getPotentialEnergy()))" in case
