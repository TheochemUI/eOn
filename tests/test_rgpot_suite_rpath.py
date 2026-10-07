"""The rgpot suite loads the library this build just linked."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MESON = ROOT / "client" / "meson.build"


def test_build_tree_rpath_is_recorded_before_the_prefix():
    text = MESON.read_text(encoding="utf-8")
    assert "build-tree rpath entries before the prefix libdir" in text
    assert "'-Wl,-rpath,'" in text
    gate = text.split(
        "if _rgpot_mpi_dep.found() and _rgpot_mpirun.found()", 1
    )[1].split("endif", 1)[0]
    assert "test_rgpot_mpi_fault" in gate
    assert "test_rgpot_mpi_abort" in gate
    assert "test_rgpot_mpi_params" in gate
