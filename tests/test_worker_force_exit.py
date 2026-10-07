"""A force error reaches rank 0. Workers still stop with status 0."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "client" / "unit_tests" / "run_rgpot_mpi.sh"
POT = ROOT / "client" / "potentials" / "Rgpot" / "RgpotPot.cpp"


def test_rank_0_must_print_the_engine_error_and_exit_nonzero():
    script = SCRIPT.read_text(encoding="utf-8")
    assert "Their status is 0" in script
    assert "rank=0 fault owner=1 energy=" in script
    assert "mpirun exited 0" in script


def test_worker_stop_is_still_exit_zero():
    text = POT.read_text(encoding="utf-8")
    start = text.index("if (hdr[0] == kStop)")
    window = text[start : start + 700]
    assert "std::exit(0)" in window
