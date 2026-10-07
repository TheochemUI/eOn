"""A rate instanton can divide each spring by a positive link weight."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_link_weights_divide_one_spring():
    header = (ROOT / "include" / "eon" / "Tunneling.h").read_text(encoding="utf-8")
    assert "std::vector<double> discretization" in header
    assert "double closedRingPotential(" in header
    source = (ROOT / "client" / "Tunneling.cpp").read_text(encoding="utf-8")
    assert "const double cj = c / wj;" in source
    assert "ro.discretization = o.discretization;" in (
        ROOT / "client" / "InstantonJob.cpp"
    ).read_text(encoding="utf-8")
    ini = (ROOT / "client" / "ParametersINI.cpp").read_text(encoding="utf-8")
    assert 'ini.Get("Instanton", "discretization"' in ini
    case = (ROOT / "client" / "unit_tests" / "TunnelingTest.cpp").read_text(
        encoding="utf-8"
    )
    assert "A link weight divides one spring of a closed ring" in case
    assert "0.75 * c" in case
    guide = (ROOT / "docs" / "source" / "user_guide" / "instanton.md").read_text(
        encoding="utf-8"
    )
    assert "one positive weight per bead" in guide
