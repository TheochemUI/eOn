"""The dense-mesh NEB spring is the geometric spring."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_geometric_spring_is_selectable():
    header = (ROOT / "include" / "eon" / "NEBSpringForce.h").read_text(
        encoding="utf-8"
    )
    ini = (ROOT / "client" / "ParametersINI.cpp").read_text(encoding="utf-8")
    assert "struct GeometricSpring" in header
    assert "geometric_spring" in ini
