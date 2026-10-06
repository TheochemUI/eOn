"""The guide matches the meson flags, the MPI client, and the install version."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OPTIONS = ROOT / "meson_options.txt"
POTENTIAL = ROOT / "docs" / "source" / "user_guide" / "potential.md"
MPI = ROOT / "docs" / "source" / "user_guide" / "mpi_potential.md"
HOME = ROOT / "docs" / "source" / "index.md"
INSTALL = ROOT / "docs" / "source" / "install" / "index.md"
PIXI = ROOT / "pixi.toml"
CLIENT = ROOT / "client" / "ClientEON.cpp"


def boolean_options(text: str) -> dict[str, str]:
    return dict(
        re.findall(
            r"(?m)^option\('([^']+)',\s*type\s*:\s*'boolean',\s*value\s*:\s*(true|false)",
            text,
        )
    )


def test_potential_flags_match_meson_options():
    opts = boolean_options(OPTIONS.read_text(encoding="utf-8"))
    for name in ("with_xtb", "with_water", "with_metatomic", "with_vasp"):
        assert opts[name] == "false"
    raw = OPTIONS.read_text(encoding="utf-8")
    assert "option('with_lammps'" not in raw
    assert "option('with_dftd3'" not in raw
    page = POTENTIAL.read_text(encoding="utf-8")
    assert "`-Dwith_xtb` defaults false" in page
    assert "`-Dwith_vasp` defaults false" in page
    assert "`-Dwith_metatomic` defaults false" in page
    assert "`with_water` defaults false" in page
    assert "The option `-Dwith_lammps` does not exist." in page
    assert "The option `-Dwith_dftd3` does not exist." in page


def test_mpi_client_names_match_the_binary():
    source = CLIENT.read_text(encoding="utf-8")
    assert 'getenv("EON_NUMBER_OF_CLIENTS")' in source
    assert 'getenv("EON_CLIENT_STANDALONE")' in source
    page = MPI.read_text(encoding="utf-8")
    assert "EON_NUMBER_OF_CLIENTS" in page
    assert "EON_CLIENT_STANDALONE" in page
    assert "eOn_NUMBER_OF_CLIENTS" not in page
    assert "eOn_CLIENT_STANDALONE" not in page
    assert "client_mpi" not in page
    assert "-n 1 eonclient" in page


def test_homepage_command_starts_from_pos_con():
    page = HOME.read_text(encoding="utf-8")
    assert "eonclient -s" not in page
    assert "The working directory needs pos.con." in page
    assert "eonclient # reads config.ini and runs" in page
    assert "python -m eon.server" in page
    assert "rgpot>=3.2.0" in page


def test_install_page_records_the_pixi_version():
    version = re.search(r'(?m)^version = "([^"]+)"', PIXI.read_text(encoding="utf-8"))
    assert version is not None
    page = INSTALL.read_text(encoding="utf-8")
    assert f"records version {version.group(1)}" in page
