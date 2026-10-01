"""CLI config init takes the positional path, not a rewritten sys.argv."""

from pathlib import Path

import pytest

from eon.config import ConfigClass, canonical_listed_value


def test_init_from_cli_uses_positional_path(tmp_path, monkeypatch):
    cfgfile = tmp_path / "config.ini"
    cfgfile.write_text("[Main]\njob = process_search\n")
    monkeypatch.chdir(tmp_path)
    cfg = ConfigClass()
    cfg.init_from_cli([str(cfgfile)])
    assert Path(cfg.config_path).resolve() == cfgfile.resolve()


def test_canonical_listed_value_keeps_exact_spelling():
    values = ["socketnwchem", "SocketNWChem", "rgpot"]
    assert canonical_listed_value("RGPOT", values) == "rgpot"
    assert canonical_listed_value("SocketNWChem", values) == "SocketNWChem"
    assert canonical_listed_value("socketnwchem", values) == "socketnwchem"
    assert canonical_listed_value("SOCKETNWCHEM", values) == "socketnwchem"
    assert canonical_listed_value("nope", values) is None


def test_server_accepts_potential_spellings(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for raw in ("RGPOT", "rgpot", "lenosky_si", "lenosky_Si"):
        cfgfile = tmp_path / f"{raw}.ini"
        cfgfile.write_text(
            "[Main]\njob = point\n[Potential]\npotential = %s\n" % raw
        )
        cfg = ConfigClass()
        cfg.init(str(cfgfile))
        assert cfg.init_done


def test_server_rejects_unknown_potential(tmp_path, monkeypatch):
    cfgfile = tmp_path / "config.ini"
    cfgfile.write_text(
        "[Main]\njob = point\n[Potential]\npotential = not-a-potential\n"
    )
    monkeypatch.chdir(tmp_path)
    cfg = ConfigClass()
    with pytest.raises(SystemExit):
        cfg.init(str(cfgfile))


def test_server_accepts_cpmd_cutoff_and_functional_aliases(tmp_path, monkeypatch):
    """Every cutoff and functional spelling the client reads loads here."""
    monkeypatch.chdir(tmp_path)
    for section in ("cpmd", "RgpotPot"):
        for key in ("cutOffRy", "cutoff_ry", "cpmd_cut_off_ry"):
            cfgfile = tmp_path / f"{section}-{key}.ini"
            aliases = f"{key} = 60.0\ncpmd_functional = PBE\n"
            rgpot = "backend = cpmdc\n" + (aliases if section == "RgpotPot" else "")
            cpmd = f"[cpmd]\n{aliases}" if section == "cpmd" else ""
            cfgfile.write_text(
                "[Main]\njob = point\n[Potential]\npotential = rgpot\n"
                f"[RgpotPot]\n{rgpot}{cpmd}"
            )
            cfg = ConfigClass()
            cfg.init(str(cfgfile))
            assert cfg.init_done
