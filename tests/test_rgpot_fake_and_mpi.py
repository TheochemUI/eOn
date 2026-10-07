"""The fake CPMD energy is pinned, and the two-rank case checks the share."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "client" / "unit_tests" / "run_rgpot_mpi.sh"
MESON = ROOT / "client" / "meson.build"
CASE = ROOT / "client" / "unit_tests" / "RgpotPotTest.cpp"


def test_fake_engine_case_requires_the_written_energy():
    text = CASE.read_text(encoding="utf-8")
    start = text.index("RgpotPot in-process cpmdc force")
    body = text[start : text.index("RgpotPot reads params_path", start)]
    assert "0.75 + 0.001 * cell_zz" in body
    assert "std::abs(energy) > 1e-6" not in body


def test_mpirun_checks_owner_energy_and_error():
    script = SCRIPT.read_text(encoding="utf-8")
    assert "mpirun -np 2" in script
    fault = script.split('if [ "$MODE" = "fault" ]', 1)[1].split(
        'if [ "$MODE" = "params" ]', 1
    )[0]
    assert "owner=1" in fault
    assert "energy=" in fault
    assert "engine-rank1" in fault
    assert 'e0=$(sed' in fault
    assert 'e1=$(sed' in fault
    assert '"$e0" != "$e1"' in fault
    meson = MESON.read_text(encoding="utf-8")
    assert "--allow-running-no-tests" in meson
