"""eOn's client MPI switch and rgpot's calculator-group switch have two names."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OPTIONS = ROOT / "meson_options.txt"
MESON = ROOT / "client" / "meson.build"
COMM = ROOT / "docs" / "source" / "user_guide" / "communicator.md"
RGPOT = ROOT / "docs" / "source" / "user_guide" / "rgpot_pot.md"


def test_meson_options_name_the_two_mpi_switches():
    text = OPTIONS.read_text(encoding="utf-8")
    assert "option('with_mpi'" in text
    assert "option('with_rgpot_mpi'" in text
    mpi = text[text.index("option('with_mpi'") : text.index("option('with_rgpot_mpi'")]
    rgpot = text[text.index("option('with_rgpot_mpi'") :]
    assert "client/server" in mpi
    assert "calculator group" in rgpot


def test_wrap_forwards_with_rgpot_mpi():
    text = MESON.read_text(encoding="utf-8")
    assert "get_option('with_rgpot_mpi').enabled()" in text
    assert "'with_mpi': 'enabled'" in text


def test_guides_tie_each_flag_to_one_launch():
    communicator = COMM.read_text(encoding="utf-8")
    rgpot = RGPOT.read_text(encoding="utf-8")
    for page in (communicator, rgpot):
        assert "-Dwith_mpi=enabled" in page
        assert "-Dwith_rgpot_mpi=enabled" in page
    start = communicator.index("mpirun -np 42")
    window = communicator[max(0, start - 500) : start + 80]
    assert "-Dwith_mpi=enabled" in window
    assert "-Dwith_rgpot_mpi=enabled" in window
    assert "calculator group" in rgpot
    assert "-Dwith_rgpot_mpi=enabled" in rgpot
