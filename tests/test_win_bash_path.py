"""Git bash PATH helpers keep sh and drop the CMake-probe breakers."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "ci" / "win_bash_path.sh"


def _bash(func: str, arg: str) -> str:
    return subprocess.check_output(
        [
            "bash",
            "-c",
            f'source "{SCRIPT}"; {func} "$1"',
            "bash",
            arg,
        ],
        text=True,
    ).strip()


def test_posix_path_is_unchanged_without_cygpath():
    assert _bash("eon_posix_path", "/c/hostedtoolcache/Python") == (
        "/c/hostedtoolcache/Python"
    )


def test_filter_drops_mingw_and_strawberry_and_keeps_git_sh():
    raw = ":".join(
        [
            "/c/mingw64/bin",
            "/c/Program Files/Git/usr/bin",
            "/c/Strawberry/perl/bin",
            "/c/Windows/system32",
        ]
    )
    cleaned = _bash("eon_drop_cmake_poison", raw)
    assert "/c/Program Files/Git/usr/bin" in cleaned
    assert "/c/Windows/system32" in cleaned
    assert "mingw" not in cleaned.lower()
    assert "strawberry" not in cleaned.lower()


def test_filter_puts_msvc_link_ahead_of_git_sh():
    msvc = (
        "/c/Program Files/Microsoft Visual Studio/2022/Enterprise"
        "/VC/Tools/MSVC/14.44.35207/bin/HostX64/x64"
    )
    git_sh = "/c/Program Files/Git/usr/bin"
    raw = ":".join(
        [
            "/c/Program Files/Git/ucrt64/bin",
            git_sh,
            "/c/Users/runneradmin/bin",
            msvc,
            "/c/Windows/system32",
        ]
    )
    parts = _bash("eon_drop_cmake_poison", raw).split(":")
    assert "ucrt64" not in ":".join(parts).lower()
    assert parts.index(msvc) < parts.index(git_sh)
    assert parts[parts.index(msvc) + 1] == git_sh
