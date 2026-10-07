"""The LAMMPS loader and the MPI tests share one EONMPI layout."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MESON = ROOT / "client" / "meson.build"


def test_eonmpi_define_precedes_the_lammps_library():
    text = MESON.read_text(encoding="utf-8")
    define = text.index("_args += ['-DEONMPI']")
    library = text.index("subdir('potentials/LAMMPS')")
    assert define < library
