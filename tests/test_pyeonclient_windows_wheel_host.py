"""Windows pyeonclient wheels need a host that CMake can probe."""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOST = ROOT / "scripts" / "pyeonclient_cibw_windows_host.py"
BASH = ROOT / "scripts" / "ci" / "win_bash_path.sh"
WORKFLOW = ROOT / ".github" / "workflows" / "pyeonclient-wheels.yml"
MESON = ROOT / "client" / "meson.build"
INIT = ROOT / "client" / "python" / "pyeonclient" / "__init__.py"
PYPROJECT = ROOT / "pyproject-pyeonclient.toml"
WRAP = ROOT / "subprojects" / "rgpot.wrap"


def _host():
    spec = importlib.util.spec_from_file_location("wheel_host", HOST)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_cleaned_path_keeps_the_prefix_and_bash_after_msvc(tmp_path: Path):
    host = _host()
    prefix = tmp_path / "prefix"
    (prefix / "Library" / "bin").mkdir(parents=True)
    (prefix / "bin").mkdir()
    path = ":".join(
        [
            "/usr/bin",
            "/c/mingw64/bin",
            "/c/Strawberry/perl/bin",
            "/c/Program Files/Microsoft Visual Studio/2022/Enterprise/VC/Tools/MSVC/14.44/bin/Hostx64/x64",
            "/c/Windows/System32",
        ]
    )
    cleaned = host.cleaned_path(path, prefix)
    low = cleaned.lower()
    assert "mingw" not in low
    assert "strawberry" not in low
    assert low.index("hostx64/x64") < low.index("/usr/bin")
    assert str(prefix / "Library" / "bin") in cleaned


def test_host_environment_points_at_flang_and_eigen(tmp_path: Path):
    host = _host()
    prefix = tmp_path / "prefix"
    eigen = prefix / "Library" / "include" / "eigen3" / "Eigen"
    eigen.mkdir(parents=True)
    (eigen / "Core").write_text("/* eigen */\n", encoding="utf-8")
    bindir = prefix / "Library" / "bin"
    bindir.mkdir()
    (bindir / "flang-new.exe").write_text("", encoding="utf-8")
    (bindir / "pkg-config.exe").write_text("", encoding="utf-8")
    pc = prefix / "Library" / "lib" / "pkgconfig"
    pc.mkdir(parents=True)
    (pc / "eigen3.pc").write_text("Name: eigen3\n", encoding="utf-8")
    cmake = prefix / "Library" / "share" / "eigen3" / "cmake"
    cmake.mkdir(parents=True)
    (cmake / "Eigen3Config.cmake").write_text("# config\n", encoding="utf-8")
    values = host.host_environment(
        prefix,
        {"PATH": "/c/mingw64/bin:/usr/bin", "LIB": "C:\\flang\\lib"},
    )
    assert values["FC"].endswith("flang-new.exe")
    assert values["CC"] == "cl"
    assert values["CXX"] == "cl"
    assert "mingw" not in values["PATH"].lower()
    assert values["Eigen3_DIR"] == str(cmake)
    assert str(pc) in values["PKG_CONFIG_PATH"]
    assert "C:\\flang\\lib" in values["LIB"]
    block = host.github_env_block(values)
    assert "\nFC=" in "\n" + block
    assert "\n" not in values["PATH"]


def test_bash_poison_filter_matches_the_python_helper():
    script = f"""
      set -euo pipefail
      source {BASH}
      eon_drop_cmake_poison "/usr/bin:/c/mingw64/bin:/c/Program Files/Microsoft Visual Studio/2022/Enterprise/VC/Tools/MSVC/14.44/bin/Hostx64/x64:/c/Windows/System32"
    """
    out = subprocess.check_output(["bash", "-c", script], text=True)
    low = out.lower()
    assert "mingw" not in low
    assert low.index("hostx64/x64") < low.index("/usr/bin")


def test_workflow_installs_the_wheel_host_before_cibuildwheel():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "scripts/ci/wheel-host/pixi.toml" in text
    assert "pyeonclient_cibw_windows_host.sh" in text
    assert "pyeonclient_cibw_before_windows.py" in text
    assert text.index("Install Windows wheel host tools") < text.index(
        "pypa/cibuildwheel@"
    )
    pixi = (ROOT / "scripts" / "ci" / "wheel-host" / "pixi.toml").read_text(
        encoding="utf-8"
    )
    assert 'platforms = ["win-64"]' in pixi
    assert 'flang = ">=19.1.7,<20"' in pixi
    assert "eigen" in pixi
    assert "capnproto" in pixi
    assert 'rust = ">=1.88"' in pixi


def test_windows_client_stays_an_archive_and_ships_wheel_dlls():
    meson = MESON.read_text(encoding="utf-8")
    assert "eclib = static_library(" in meson
    ext = (ROOT / "client" / "potentials" / "ExtPot" / "meson.build").read_text(
        encoding="utf-8"
    )
    lammps = (ROOT / "client" / "potentials" / "LAMMPS" / "meson.build").read_text(
        encoding="utf-8"
    )
    assert "extpot = static_library(" in ext
    assert "lammps_pot = static_library(" in lammps
    py = (ROOT / "client" / "python" / "meson.build").read_text(encoding="utf-8")
    assert "_rgpot_wheel_name = 'rgpot.dll'" in py
    assert "output: 'readcon_db.dll'" in py
    init = INIT.read_text(encoding="utf-8")
    assert init.index("_extend_windows_dll_search()") < init.index(
        "from pyeonclient._core import"
    )
    assert "os.add_dll_directory" in init
    project = PYPROJECT.read_text(encoding="utf-8")
    assert "allow-windows-internal-shared-libs = true" in project
    py_meson = (ROOT / "client" / "python" / "meson.build").read_text(
        encoding="utf-8"
    )
    assert "cppc.get_supported_arguments(['-Wno-unused-parameter'])" in py_meson
    wrap = WRAP.read_text(encoding="utf-8")
    assert "rgpot-windows-link-whole.patch" in wrap
    assert "rgpot-windows-flang-rt.patch" in wrap
    whole = (
        ROOT / "subprojects" / "packagefiles" / "rgpot-windows-link-whole.patch"
    ).read_text(encoding="utf-8")
    assert "link_whole: [rgpot_core]" in whole
    flang = (
        ROOT / "subprojects" / "packagefiles" / "rgpot-windows-flang-rt.patch"
    ).read_text(encoding="utf-8")
    assert "FortranRuntime.dynamic.lib" in flang
