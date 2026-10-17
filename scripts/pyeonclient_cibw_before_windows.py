#!/usr/bin/env python3
"""Install pyeonclient's Windows wheel build dependencies.

cibuildwheel runs CIBW_BEFORE_BUILD_WINDOWS through cmd.exe. A multiline
value executes only the first line, so numpy (a pyproject build-system
requirement under --no-isolation) was never installed.

The variant is an argument because the build environment does not inherit
every job variable. It is ``base`` or ``metatomic``.
"""

from __future__ import annotations

import subprocess
import sys


def run(args: list[str]) -> None:
    print("+", " ".join(args), flush=True)
    subprocess.check_call(args)


def main(argv: list[str]) -> int:
    variant = argv[1] if len(argv) > 1 else ""
    py = sys.executable
    run(
        [
            py,
            "-m",
            "pip",
            "install",
            "-U",
            "pip",
            "nanobind>=2.4,<3",
            # 1.12 adds -Werror=implicit-function-declaration to every flang
            # probe. flang 19 rejects that flag, so FortranRuntime looks missing.
            "meson>=1.10.1,<1.11",
            "ninja",
            "meson-python",
        ]
    )
    run([py, "-m", "pip", "install", "--only-binary=:all:", "numpy"])
    if variant == "metatomic":
        run(
            [
                py,
                "-m",
                "pip",
                "install",
                "--only-binary=:all:",
                "torch",
                "--index-url",
                "https://download.pytorch.org/whl/cpu",
            ]
        )
        run(
            [
                py,
                "-m",
                "pip",
                "install",
                "metatensor",
                "metatensor-torch",
                "metatomic-torch",
                "vesin",
            ]
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
