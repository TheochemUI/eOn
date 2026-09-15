"""CLI config init takes the positional path, not a rewritten sys.argv."""

from pathlib import Path

from eon.config import ConfigClass


def test_init_from_cli_uses_positional_path(tmp_path, monkeypatch):
    cfgfile = tmp_path / "config.ini"
    cfgfile.write_text("[Main]\njob = process_search\n")
    monkeypatch.chdir(tmp_path)
    cfg = ConfigClass()
    cfg.init_from_cli([str(cfgfile)])
    assert Path(cfg.config_path).resolve() == cfgfile.resolve()
