"""Wheel repair renames librgpot.so.4 even when the file is misnamed."""

from __future__ import annotations

import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPAIR = ROOT / "scripts" / "pyeonclient_repair_wheel.sh"


def _tool(name: str) -> bool:
    return shutil.which(name) is not None


@pytest.mark.skipif(
    not (_tool("gcc") and _tool("patchelf") and _tool("readelf")),
    reason="gcc, patchelf and readelf build the synthetic wheel",
)
def test_repair_renames_a_misnamed_so4(tmp_path: Path):
    """rgpot 3.5.1's file is librgpot.so.3.5.1 and its SONAME is .so.4.

    A wheel that stored that file as librgpot.so.3 left the NEEDED
    entry unresolved. Repair vendors it and renames it to
    libeon_rgpot.so.4.
    """
    src = tmp_path / "src"
    src.mkdir()
    (src / "wrap.c").write_text("void rgpot_marker(void) {}\n")
    (src / "core.c").write_text("void rgpot_marker(void); void core_marker(void) { rgpot_marker(); }\n")
    subprocess.check_call(
        [
            "gcc",
            "-shared",
            "-fPIC",
            "-o",
            "librgpot.so.3.5.1",
            "wrap.c",
            "-Wl,-soname,librgpot.so.4",
        ],
        cwd=src,
    )
    subprocess.check_call(
        ["ln", "-s", "librgpot.so.3.5.1", "librgpot.so.4"],
        cwd=src,
    )
    subprocess.check_call(
        [
            "gcc",
            "-shared",
            "-fPIC",
            "-o",
            "_core.so",
            "core.c",
            "-L.",
            "-l:librgpot.so.4",
            "-Wl,-rpath,$ORIGIN",
        ],
        cwd=src,
    )
    pkg = tmp_path / "pkg"
    (pkg / "pyeonclient").mkdir(parents=True)
    (pkg / "pyeonclient-0.0.0.dist-info").mkdir()
    shutil.copy(src / "librgpot.so.3.5.1", pkg / "pyeonclient" / "librgpot.so.3")
    shutil.copy(src / "_core.so", pkg / "pyeonclient" / "_core.so")
    (pkg / "pyeonclient-0.0.0.dist-info" / "RECORD").write_text("")
    whl = tmp_path / "pyeonclient-0.0.0-py3-none-any.whl"
    with zipfile.ZipFile(whl, "w") as z:
        for path in pkg.rglob("*"):
            if path.is_file():
                z.write(path, path.relative_to(pkg).as_posix())

    env = os.environ.copy()
    env["PYEONCLIENT_BUILD_ROOT"] = str(tmp_path)
    proc = subprocess.run(
        ["bash", str(REPAIR), str(whl)],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    out = tmp_path / "out"
    with zipfile.ZipFile(whl) as z:
        z.extractall(out)
    libs = out / ".pyeonclient.mesonpy.libs"
    private = libs / "libeon_rgpot.so.4"
    assert private.is_file()
    assert not list((out / "pyeonclient").glob("librgpot.so*"))
    core_need = subprocess.check_output(
        ["readelf", "-d", str(out / "pyeonclient" / "_core.so")],
        text=True,
    )
    soname = subprocess.check_output(
        ["readelf", "-d", str(private)],
        text=True,
    )
    assert "[libeon_rgpot.so.4]" in core_need
    assert "librgpot.so" not in core_need
    assert "Library soname: [libeon_rgpot.so.4]" in soname
