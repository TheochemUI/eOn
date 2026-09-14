#!/usr/bin/env python3
"""Ensure _core.__version__ is derived from pyproject-pyeonclient.toml."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib  # type: ignore

    meta = tomllib.loads((ROOT / "pyproject-pyeonclient.toml").read_text())
    v_py = meta["project"]["version"]

    mod = (ROOT / "client/python/bind/module.cpp").read_text()
    if re.search(r'__version__"\)\s*=\s*"', mod):
        print(
            "ERROR: module.cpp hardcodes __version__; it must publish "
            "PYEONCLIENT_VERSION",
            file=sys.stderr,
        )
        return 1
    if not re.search(r'__version__"\)\s*=\s*PYEONCLIENT_VERSION', mod):
        print("ERROR: no __version__ = PYEONCLIENT_VERSION in module.cpp", file=sys.stderr)
        return 1

    build = (ROOT / "client/python/meson.build").read_text()
    if "pyproject-pyeonclient.toml" not in build or "PYEONCLIENT_VERSION" not in build:
        print(
            "ERROR: client/python/meson.build does not define "
            "PYEONCLIENT_VERSION from pyproject-pyeonclient.toml",
            file=sys.stderr,
        )
        return 1

    print(f"pyeonclient version OK: {v_py} (generated into _core)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
