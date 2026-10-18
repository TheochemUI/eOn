"""The Windows wheel host drops the tools that break the CMake probe."""

from __future__ import annotations

import subprocess
from pathlib import Path

BASH = Path(__file__).resolve().parents[1] / "scripts" / "ci" / "win_bash_path.sh"


def _drop(path: str) -> str:
    script = f"""
      set -euo pipefail
      source {BASH}
      eon_drop_cmake_poison {path!r}
    """
    return subprocess.check_output(["bash", "-c", script], text=True)


def test_mingw_and_strawberry_are_removed_and_bash_follows_msvc():
    out = _drop(
        ":".join(
            [
                "/usr/bin",
                "/c/mingw64/bin",
                "/c/Strawberry/c/bin",
                "/c/msys64/ucrt64/bin",
                "/c/Program Files/Microsoft Visual Studio/2022/Enterprise/VC/Tools/MSVC/14.44/bin/Hostx64/x64",
                "/c/Windows/System32",
            ]
        )
    )
    low = out.lower()
    assert "mingw" not in low
    assert "strawberry" not in low
    assert "ucrt64" not in low
    assert "/usr/bin" in low
    assert low.index("hostx64/x64") < low.index("/usr/bin")
    assert low.endswith("/c/Windows/System32") or "/c/windows/system32" in low
