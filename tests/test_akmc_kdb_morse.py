"""Morse aKMC stores the first good process and suggests it on the next search."""

from pathlib import Path

import pytest


def test_morse_akmc_inserts_then_suggests(tmp_path, monkeypatch, eon):
    pytest.importorskip("amsel")
    pytest.importorskip("readcon_db")
    ini = Path("tests/test_akmc_pt/morse_dimer.ini").read_text()
    ini += (
        "\n[KDB]\n"
        "use_kdb = true\n"
        "kdb_nf = 0.2\n"
        "kdb_dc = 0.3\n"
        "kdb_mac = 0.7\n"
    )
    pos = Path("tests/data/server/Pt_Heptamer_oneLayer/pos.con")
    (tmp_path / "config.ini").write_text(ini)
    (tmp_path / "pos.con").write_bytes(pos.read_bytes())
    monkeypatch.chdir(tmp_path)
    log_path = tmp_path / "akmc.log"
    for _ in range(6):
        eon()
        if log_path.is_file() and "Made a KDB suggestion" in log_path.read_text():
            break
    text = log_path.read_text()
    assert "Python module kdb not found" not in text
    assert "Made a KDB suggestion" in text
    assert "readcon.db" in text
    assert (tmp_path / "kdb" / "data.mdb").is_file()
    assert (tmp_path / "readcon.db" / "data.mdb").is_file()
