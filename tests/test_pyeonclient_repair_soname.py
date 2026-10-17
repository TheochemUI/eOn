"""The wheel repair must vendor librgpot when the filename is not the SONAME.

The wrap writes librgpot.so.<project version> (librgpot.so.3.5.1) and the
loader NEEDED entry is librgpot.so.<kCxxSoversion> (librgpot.so.4). The
pyeonclient custom target used to install that file as librgpot.so.3.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPAIR = ROOT / "scripts" / "pyeonclient_repair_wheel.sh"


def _have_tools() -> bool:
    return all(shutil.which(tool) for tool in ("gcc", "readelf", "patchelf"))


def _so_needed(path: Path) -> str:
    return subprocess.check_output(["readelf", "-d", str(path)], text=True)


@pytest.mark.skipif(not _have_tools(), reason="gcc, readelf, and patchelf are required")
def test_repair_vendors_soname_when_the_member_name_is_old(tmp_path: Path) -> None:
    src = tmp_path / "src"
    src.mkdir()
    (src / "rgpot.c").write_text("void rgpot_anchor(void) {}\n")
    (src / "core.c").write_text(
        "void rgpot_anchor(void);\nvoid core_anchor(void) { rgpot_anchor(); }\n"
    )
    real = src / "librgpot.so.3.5.1"
    subprocess.check_call(
        [
            "gcc",
            "-shared",
            "-fPIC",
            "-Wl,-soname,librgpot.so.4",
            "-o",
            str(real),
            str(src / "rgpot.c"),
        ]
    )
    core = src / "_core.abi3.so"
    subprocess.check_call(
        [
            "gcc",
            "-shared",
            "-fPIC",
            "-Wl,--no-as-needed",
            str(real),
            "-Wl,--as-needed",
            "-Wl,-soname,_core.abi3.so",
            "-o",
            str(core),
            str(src / "core.c"),
        ]
    )
    assert "librgpot.so.4" in _so_needed(core)
    assert "librgpot.so.4" in _so_needed(real)

    pkg = tmp_path / "pkg"
    (pkg / "pyeonclient").mkdir(parents=True)
    shutil.copy2(core, pkg / "pyeonclient" / "_core.abi3.so")
    # The CI wheel member was the real library under the old soversion name.
    shutil.copy2(real, pkg / "pyeonclient" / "librgpot.so.3")
    dist = pkg / "pyeonclient-0.0.0.dist-info"
    dist.mkdir()
    (dist / "METADATA").write_text("Name: pyeonclient\nVersion: 0.0.0\n")
    (dist / "WHEEL").write_text("Wheel-Version: 1.0\nRoot-Is-Purelib: false\n")
    (dist / "RECORD").write_text("pyeonclient/_core.abi3.so,,\n")

    wheel = tmp_path / "pyeonclient-0.0.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        for path in pkg.rglob("*"):
            if path.is_file():
                archive.write(path, path.relative_to(pkg).as_posix())

    env = os.environ.copy()
    env["PYEONCLIENT_BUILD_ROOT"] = str(tmp_path / "missing-build")
    proc = subprocess.run(
        ["bash", str(REPAIR), str(wheel)],
        cwd=str(tmp_path),
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.fail(proc.stdout + "\n" + proc.stderr)

    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        assert any(name.endswith("libeon_rgpot.so.4") for name in names), names
        assert not any(Path(name).name.startswith("librgpot.so") for name in names), names
        extracted = tmp_path / "out"
        archive.extractall(extracted)
    bundled = next(extracted.rglob("_core.abi3.so"))
    needed = _so_needed(bundled)
    assert "libeon_rgpot.so.4" in needed
    assert "librgpot.so" not in needed
