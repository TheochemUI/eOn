"""One rank that cannot read params_path stops the others before the split."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "client" / "potentials" / "Rgpot" / "RGPotEngine.cpp"
SCRIPT = ROOT / "client" / "unit_tests" / "run_rgpot_mpi.sh"
GUIDE = ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md"


def test_agreement_runs_without_the_eon_mpi_compile_flag():
    text = ENGINE.read_text(encoding="utf-8")
    start = text.index("bool agree_construction")
    window = text[max(0, start - 180) : start]
    assert "EON_RGPOT_MPI" not in window
    agree = text.index("agree_construction(local_error)")
    split = text.index("bindCalculators(opt.ranks_per_image)")
    assert agree < split
    assert "mpi_world_hint() > 1 && !agree_construction(local_error)" in text


def test_a_rank_left_in_the_split_fails_the_check():
    script = SCRIPT.read_text(encoding="utf-8")
    assert "MPI_Comm_split" in script
    assert 'rank=0 params-agreed .*cannot open params_path' in script


def test_guide_separates_the_two_mpi_builds():
    guide = GUIDE.read_text(encoding="utf-8")
    assert "builds the client/server program, not" in guide
    assert "-Dwith_rgpot_mpi=enabled" in guide
    assert "mpirun -np N eonclient" in guide
