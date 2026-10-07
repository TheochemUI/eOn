"""The band carries a harmonic quantum free energy with the tangent removed."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_quantum_free_energy_is_wired_through_the_band():
    header = (ROOT / "include" / "eon" / "QuantumFreeEnergy.h").read_text(
        encoding="utf-8"
    )
    assert "perpendicularHarmonicFreeEnergy" in header
    assert "quantumFreeEnergies" in header
    source = (ROOT / "client" / "QuantumFreeEnergy.cpp").read_text(
        encoding="utf-8"
    )
    assert "kHbar" in source
    assert "std::sinh" in source
    case = (
        ROOT / "client" / "unit_tests" / "QuantumFreeEnergyTest.cpp"
    ).read_text(encoding="utf-8")
    assert "curved valley free energy matches a direct sum" in case
    assert "perpendicularHarmonicFreeEnergy(hessian, tangent, 0.0)" in case
    job = (ROOT / "client" / "NudgedElasticBandJob.cpp").read_text(
        encoding="utf-8"
    )
    assert "quantumFreeEnergies(" in job
    assert "quantum_temperature" in job
    frames = (ROOT / "client" / "NEBSplineExtrema.cpp").read_text(
        encoding="utf-8"
    )
    assert "quantum_free_energy" in frames
    assert "zpe_corrected_barrier" in frames
