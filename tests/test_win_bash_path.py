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
