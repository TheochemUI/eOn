"""A rate instanton takes one positive time step per link: an adaptive grid."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_link_weights_are_time_steps():
    header = (ROOT / "include" / "eon" / "Tunneling.h").read_text(encoding="utf-8")
    assert "std::vector<double> discretization" in header
    assert "double closedRingPotential(" in header
    source = (ROOT / "client" / "Tunneling.cpp").read_text(encoding="utf-8")
    assert "const double cj = c / wj;" in source
    assert "const double share = 0.5 * (wj + wp);" in source
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
    assert "A ring of unequal time steps converges to the continuum instanton" in case
    guide = " ".join(
        (ROOT / "docs" / "source" / "user_guide" / "instanton.md")
        .read_text(encoding="utf-8")
        .split()
    )
    assert "one positive weight per bead" in guide
    assert "inst-rommelAdaptiveIntegrationGrids2011" in guide
