"""A missing config does not start the client logger."""

import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_missing_config_is_rejected_before_the_logger():
    source = (ROOT / "client" / "ClientEON.cpp").read_text(encoding="utf-8")
    guard = source.find("problem loading parameter file, stopping")
    logger = source.find("eonc::log::init_client()")
    assert guard != -1 and logger != -1
    assert guard < logger

    binary = os.environ.get("EONCLIENT")
    if not binary:
        return
    work = Path(os.environ.get("TMPDIR", "/tmp")) / "eon-startup-missing"
    work.mkdir(exist_ok=True)
    for name in ("client_quill.log", "client_traceback.log", "config.ini"):
        path = work / name
        if path.exists():
            path.unlink()
    proc = subprocess.run(
        [binary],
        cwd=work,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert "Can't load INI file: config.ini" in proc.stderr
    assert not (work / "client_quill.log").exists()
    assert not (work / "client_traceback.log").exists()
