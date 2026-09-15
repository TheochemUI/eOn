"""Integration tests and sample configs must not embed workstation home paths."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TRACKED = [
    ROOT / "tests" / "test_pyeonclient_metatomic_cookbook.py",
    ROOT / "tests" / "test_pyeonclient_neb.py",
    ROOT / "tests" / "test_pyeonclient_rgpot_metatomic.py",
    ROOT / "client" / "potentials" / "ASE_NWCHEM" / "client.py",
    ROOT / "client" / "potentials" / "ASE_NWCHEM" / "blah.py",
    ROOT / "client" / "potentials" / "ASE_NWCHEM" / "design.org",
    ROOT / "client" / "potentials" / "ASE_ORCA" / "design.org",
    ROOT / "client" / "potentials" / "Metatomic" / "readme.org",
    ROOT / "client" / "potentials" / "SocketNWChem" / "design" / "config.ini",
    ROOT / "client" / "tests" / "bgsd" / "config.ini",
    ROOT / "client" / "tests" / "lj_cluster" / "config.ini",
    ROOT / "tools" / "clusters" / "bgp" / "surveyor.sh",
]


def test_tracked_files_have_no_home_user_paths():
    offenders: list[str] = []
    for path in TRACKED:
        text = path.read_text(encoding="utf-8")
        for i, line in enumerate(text.splitlines(), 1):
            if "/home/" in line or "/Users/" in line:
                offenders.append(f"{path.relative_to(ROOT)}:{i}:{line.strip()}")
    assert not offenders, "workstation paths remain:\n" + "\n".join(offenders)


def test_cookbook_tests_name_env_skip():
    neb = (ROOT / "tests" / "test_pyeonclient_neb.py").read_text(encoding="utf-8")
    cook = (ROOT / "tests" / "test_pyeonclient_metatomic_cookbook.py").read_text(
        encoding="utf-8"
    )
    rg = (ROOT / "tests" / "test_pyeonclient_rgpot_metatomic.py").read_text(
        encoding="utf-8"
    )
    assert "EON_PET_NEB_ROOT" in neb
    assert "EON_PET_NEB_ROOT" in rg
    assert "EON_PET_MAD_POS" in cook or "EON_PET_NEB_ROOT" in cook
    nw = (ROOT / "client" / "potentials" / "ASE_NWCHEM" / "client.py").read_text(
        encoding="utf-8"
    )
    assert "NWCHEM_COMMAND" in nw
