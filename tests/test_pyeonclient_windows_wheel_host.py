"""The Windows wheel host drops the tools that break the MSVC CMake probe."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "pyeonclient_cibw_windows_host",
    ROOT / "scripts" / "pyeonclient_cibw_windows_host.py",
)
hostmod = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(hostmod)

WORKFLOW = ROOT / ".github" / "workflows" / "pyeonclient-wheels.yml"
HOST_SH = ROOT / "scripts" / "pyeonclient_cibw_windows_host.sh"


def test_msys_prefix_becomes_a_windows_path():
    assert (
        hostmod.msys_to_windows("/d/a/eOn/eOn/scripts/ci/wheel-host", enable=True)
        == "D:\\a\\eOn\\eOn\\scripts\\ci\\wheel-host"
    )
    assert hostmod.msys_to_windows("/usr/bin", enable=False) == "/usr/bin"
    assert hostmod.msys_to_windows("D:\\already", enable=True) == "D:\\already"


def test_poison_paths_are_the_cmake_probe_breakers():
    assert hostmod.is_poison_path(r"C:\mingw64\bin")
    assert hostmod.is_poison_path(r"C:\Strawberry\perl\bin")
    assert hostmod.is_poison_path(r"C:\Program Files\Git\ucrt64\bin")
    assert not hostmod.is_poison_path(r"C:\Program Files\Git\usr\bin")
    assert not hostmod.is_poison_path(r"C:\Program Files\Git\bin")
    assert not hostmod.is_poison_path(
        r"C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.44\bin\HostX64\x64"
    )
    assert hostmod.split_path(r"C:\Windows Kits\Include") == [
        r"C:\Windows Kits\Include"
    ]


def test_cleaned_path_keeps_msvc_and_the_host_bins(tmp_path: Path):
    prefix = tmp_path / "env"
    raw = ";".join(
        [
            r"C:\mingw64\bin",
            r"C:\Program Files\Git\usr\bin",
            r"C:\Strawberry\c\bin",
            r"C:\Program Files\Git\bin",
            r"C:\Windows\System32",
        ]
    )
    cleaned = hostmod.cleaned_path(raw, prefix)
    assert str(prefix / "Library" / "bin") in cleaned
    assert r"C:\Program Files\Git\bin" in cleaned
    assert r"C:\Windows\System32" in cleaned
    assert "mingw" not in cleaned.lower()
    assert "strawberry" not in cleaned.lower()
    assert r"C:\Program Files\Git\usr\bin" in cleaned


def test_cleaned_path_puts_msvc_link_ahead_of_git_sh(tmp_path: Path):
    prefix = tmp_path / "env"
    msvc = (
        r"C:\Program Files\Microsoft Visual Studio\2022\Enterprise"
        r"\VC\Tools\MSVC\14.44.35207\bin\HostX64\x64"
    )
    git_sh = r"C:\Program Files\Git\usr\bin"
    raw = ";".join(
        [
            r"C:\Program Files\Git\ucrt64\bin",
            git_sh,
            r"C:\Users\runneradmin\bin",
            msvc,
            r"C:\Windows\System32",
        ]
    )
    parts = hostmod.cleaned_path(raw, prefix).split(";")
    assert "ucrt64" not in hostmod.cleaned_path(raw, prefix).lower()
    assert parts.index(msvc) < parts.index(git_sh)
    assert parts.index(str(prefix / "Library" / "bin")) < parts.index(msvc)


def test_host_environment_points_at_eigen_and_flang(tmp_path: Path):
    prefix = tmp_path / "env"
    (prefix / "Library" / "include" / "eigen3" / "Eigen").mkdir(parents=True)
    (prefix / "Library" / "include" / "eigen3" / "Eigen" / "Core").write_text(
        "", encoding="utf-8"
    )
    (prefix / "Library" / "share" / "eigen3" / "cmake").mkdir(parents=True)
    (prefix / "Library" / "share" / "eigen3" / "cmake" / "Eigen3Config.cmake").write_text(
        "", encoding="utf-8"
    )
    (prefix / "Library" / "lib" / "pkgconfig").mkdir(parents=True)
    share_pc = prefix / "Library" / "share" / "pkgconfig"
    share_pc.mkdir(parents=True)
    (share_pc / "eigen3.pc").write_text("Name: Eigen3\n", encoding="utf-8")
    (prefix / "Library" / "bin").mkdir(parents=True)
    (prefix / "Library" / "bin" / "pkg-config.exe").write_text("", encoding="utf-8")
    (prefix / "Library" / "bin" / "flang-new.exe").write_text("", encoding="utf-8")

    assert hostmod.find_eigen_header(prefix) is not None
    assert hostmod.find_eigen_cmake_dir(prefix) is not None
    assert hostmod.find_tool(prefix, ("flang-new",)) is not None

    values = hostmod.host_environment(
        prefix,
        {
            "PATH": r"C:\mingw64\bin;C:\Windows\System32",
            "LIB": r"C:\flang\lib;C:\Windows Kits\Lib",
            "INCLUDE": r"C:\Windows Kits\Include",
        },
    )
    assert values["FC"].endswith("flang-new.exe")
    assert values["CC"] == "cl"
    assert values["CXX"] == "cl"
    assert values["CMAKE_GENERATOR"] == "Ninja"
    assert values["MSYS2_ARG_CONV_EXCL"] == "*"
    assert "mingw" not in values["PATH"].lower()
    pc_path = values["PKG_CONFIG_PATH"].replace("\\", "/")
    assert "Library/lib/pkgconfig" in pc_path
    assert pc_path.endswith("Library/share/pkgconfig") or (
        "Library/share/pkgconfig" in pc_path
    )
    assert values["Eigen3_DIR"].endswith("cmake")
    assert values["LIB"].startswith(r"C:\flang\lib")
    assert str(prefix / "Library" / "lib") in values["LIB"]
    assert values["INCLUDE"].startswith(
        str(prefix / "Library" / "include" / "eigen3")
    )
    block = hostmod.github_env_block(values)
    assert "\nFC=" in "\n" + block
    assert "PKG_CONFIG=" in block


def test_workflow_installs_the_host_prefix_before_cibuildwheel():
    text = WORKFLOW.read_text(encoding="utf-8")
    shell = HOST_SH.read_text(encoding="utf-8")
    assert "pixi-version: v0.82.0" in text
    assert "scripts/ci/wheel-host/pixi.toml" in text
    assert "pyeonclient_cibw_windows_host.sh" in text
    assert "activate-environment: false" in text
    host = text.split("pyeonclient_cibw_windows_host.sh", 1)[0]
    assert "ilammy/msvc-dev-cmd" in host
    assert shell.index("eon_drop_cmake_poison") < shell.index("eon_export_flang_rt")
    assert "win_bash_path.sh" in shell
    assert "win_flang_rt.sh" in shell
    assert "flang-new" in shell
    assert "cargo install --root" in shell
    host_toml = (ROOT / "scripts" / "ci" / "wheel-host" / "pixi.toml").read_text(
        encoding="utf-8"
    )
    assert 'rust = ">=1.88"' in host_toml
    assert 'flang = ">=19.1.7,<20"' in host_toml
