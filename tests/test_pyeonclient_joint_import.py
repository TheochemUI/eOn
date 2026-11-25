"""The joint-import check needs the bundled libeon_rgpot.so.N itself mapped."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "pyeonclient_joint_import.py"


@pytest.fixture
def joint(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    spec = importlib.util.spec_from_file_location("joint_import", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    libs = tmp_path / "site" / ".pyeonclient.mesonpy.libs"
    package = tmp_path / "site" / "pyeonclient"
    rg = tmp_path / "site" / "rgpot"
    for d in (libs, package, rg):
        d.mkdir(parents=True)
    core = package / "_core.abi3.so"
    eonclib = libs / "libeonclib.so"
    private = libs / "libeon_rgpot.so.6"
    rg_core = rg / "_rgpot.so"
    for so in (core, eonclib, private, rg_core):
        so.write_bytes(b"")

    needs = {
        core: [eonclib.name],
        eonclib: [private.name, "libc.so.6"],
        private: ["libc.so.6"],
        rg_core: ["librgpot.so.6"],
    }
    wheels = {"pyeonclient": [core, eonclib, private], "rgpot": [rg_core]}
    monkeypatch.setattr(module, "needed", lambda path: needs[path])
    monkeypatch.setattr(module, "wheel_shared_objects", lambda dist: wheels[dist])
    for name in ("pyeonclient", "rgpot"):
        fake = types.ModuleType(name)
        fake.__version__ = "0"
        fake.__file__ = str(tmp_path / "site" / name / "__init__.py")
        monkeypatch.setitem(sys.modules, name, fake)

    def run(mapped: list[Path]) -> int:
        monkeypatch.setattr(
            module, "mapped_libraries", lambda: sorted(str(p) for p in mapped)
        )
        return module.main()

    run.eonclib = eonclib
    run.private = private
    run.rg_core = rg_core
    return run


def test_mapped_libeon_rgpot_passes(joint, capsys) -> None:
    assert joint([joint.eonclib, joint.private, joint.rg_core]) == 0
    assert "JOINT_IMPORT_OK" in capsys.readouterr().out


def test_libeonclib_alone_does_not_count_as_libeon_rgpot(joint, capsys) -> None:
    assert joint([joint.eonclib, joint.rg_core]) == 1
    assert "the bundled libeon_rgpot is not loaded" in capsys.readouterr().err
