"""install_wrap_rgpot links the DT_SONAME, including over a versioned filename."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "install_wrap_rgpot.py"


def _tool(name: str) -> bool:
    return shutil.which(name) is not None


def _elf(directory: Path, filename: str, soname: str) -> None:
    source = directory / "wrap.c"
    if not source.exists():
        source.write_text("void rgpot_marker(void) {}\n")
    subprocess.check_call(
        [
            "gcc",
            "-shared",
            "-fPIC",
            "-o",
            filename,
            "wrap.c",
            f"-Wl,-soname,{soname}",
        ],
        cwd=directory,
    )


def _meson_links(directory: Path, filename: str, soname: str) -> None:
    """Real file plus the two links meson writes for a versioned library."""
    _elf(directory, filename, soname)
    subprocess.check_call(["ln", "-s", filename, soname], cwd=directory)
    unversioned = filename.split(".so", 1)[0] + ".so"
    subprocess.check_call(["ln", "-s", soname, unversioned], cwd=directory)


def _install(build: Path, prefix: Path) -> subprocess.CompletedProcess[str]:
    real = next(path for path in build.iterdir() if path.is_file() and path.name.startswith("librgpot.so"))
    env = os.environ.copy()
    env["MESON_INSTALL_DESTDIR_PREFIX"] = str(prefix)
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(real), "lib"],
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def _assert_resolves(path: Path, link_target: str) -> None:
    assert path.is_symlink()
    assert os.readlink(path) == link_target
    resolved = path.resolve()
    assert resolved.is_file()


@pytest.mark.skipif(not (_tool("gcc") and _tool("readelf")), reason="gcc and readelf build the library")
def test_install_keeps_soname_links_for_versioned_filename(tmp_path: Path):
    """The build dir holds librgpot.so.3.5.1 with SONAME librgpot.so.4.

    Meson also writes librgpot.so.4 and librgpot.so. The install keeps
    those links pointing at the regular file, and the same for
    librgpot_mpi.so.3.
    """
    build = tmp_path / "build"
    build.mkdir()
    _meson_links(build, "librgpot.so.3.5.1", "librgpot.so.4")
    _meson_links(build, "librgpot_mpi.so.3.5.1", "librgpot_mpi.so.3")
    prefix = tmp_path / "prefix"
    proc = _install(build, prefix)
    assert proc.returncode == 0, proc.stderr
    lib = prefix / "lib"
    _assert_resolves(lib / "librgpot.so.4", "librgpot.so.3.5.1")
    _assert_resolves(lib / "librgpot.so", "librgpot.so.4")
    _assert_resolves(lib / "librgpot_mpi.so.3", "librgpot_mpi.so.3.5.1")
    _assert_resolves(lib / "librgpot_mpi.so", "librgpot_mpi.so.3")
    assert not (lib / "librgpot.so.3").exists()
    soname = subprocess.check_output(["readelf", "-d", str(lib / "librgpot.so.3.5.1")], text=True)
    assert "Library soname: [librgpot.so.4]" in soname
    again = _install(build, prefix)
    assert again.returncode == 0, again.stderr
    _assert_resolves(lib / "librgpot.so.4", "librgpot.so.3.5.1")
    _assert_resolves(lib / "librgpot_mpi.so.3", "librgpot_mpi.so.3.5.1")


@pytest.mark.skipif(not (_tool("gcc") and _tool("readelf")), reason="gcc and readelf build the library")
def test_install_creates_soname_link_when_the_build_dir_has_none(tmp_path: Path):
    build = tmp_path / "build"
    build.mkdir()
    _elf(build, "librgpot.so.3.5.1", "librgpot.so.4")
    prefix = tmp_path / "prefix"
    proc = _install(build, prefix)
    assert proc.returncode == 0, proc.stderr
    _assert_resolves(prefix / "lib" / "librgpot.so.4", "librgpot.so.3.5.1")
