#!/usr/bin/env python3
"""Install pyeonclient's Windows wheel build dependencies.

cibuildwheel runs CIBW_BEFORE_BUILD_WINDOWS through cmd.exe. A multiline
value executes only the first line, so numpy (a pyproject build-system
requirement under --no-isolation) was never installed.

The variant is an argument because the build environment does not inherit
every job variable. It is ``base`` or ``metatomic``. The metatomic variant
takes the pinned requirements after it; the test env installs the same pins.
"""

from __future__ import annotations

import subprocess
import sys


def run(args: list[str]) -> None:
    print("+", " ".join(args), flush=True)
    subprocess.check_call(args)


def main(argv: list[str]) -> int:
    variant = argv[1] if len(argv) > 1 else ""
    pins = argv[2:]
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
            "meson>=1.3",
            "ninja",
            "meson-python",
        ]
    )
    run([py, "-m", "pip", "install", "--only-binary=:all:", "numpy"])
    if variant == "metatomic":
        torch_pins = [p for p in pins if p.startswith("torch=")]
        if not torch_pins:
            raise SystemExit("metatomic variant needs a torch== pin")
        run(
            [
                py,
                "-m",
                "pip",
                "install",
                "--only-binary=:all:",
                *torch_pins,
                "--index-url",
                "https://download.pytorch.org/whl/cpu",
            ]
        )
        run([py, "-m", "pip", "install", *pins])
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
