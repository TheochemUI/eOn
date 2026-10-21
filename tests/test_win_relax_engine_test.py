"""The relax-engine launcher puts linked bins ahead of PATH."""

import importlib.util
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    path = ROOT / "scripts" / "ci" / "win_relax_engine_test.py"
    spec = importlib.util.spec_from_file_location("win_relax_engine_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_front_dirs_lead_path_and_missing_dirs_drop(tmp_path: Path):
    mod = _mod()
    present = tmp_path / "torch"
    present.mkdir()
    missing = tmp_path / "gone"
    got = mod.path_with_front([str(present), str(missing), ""], "/usr/bin")
    assert got.split(os.pathsep) == [str(present), "/usr/bin"]


def test_relax_engine_path_is_wrapped():
    text = (ROOT / "client" / "meson.build").read_text()
    assert "win_relax_engine_test.py" in text
    assert "_relax_front" in text
