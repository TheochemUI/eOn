#!/usr/bin/env python3
"""Run test_relax_engine with the linked torch bins ahead of PATH.

Meson prepends metatomic wheel bins. A c10.dll from that directory is
missing an entry point, and the process exits 0xC0000139. Directories
after the executable are placed at the front of PATH. They are not
copied beside the other tests.
"""
from __future__ import annotations

import os
import sys


def path_with_front(dirs: list[str], old: str) -> str:
    """Directories that exist, then the previous PATH."""
    kept = [d for d in dirs if d and os.path.isdir(d)]
    parts = kept + ([old] if old else [])
    return os.pathsep.join(parts)


def main() -> int:
    if len(sys.argv) < 2:
        print(
            "usage: win_relax_engine_test.py EXE [DIR ...]",
            file=sys.stderr,
        )
        return 2
    exe = sys.argv[1]
    old = os.environ.get("PATH", "")
    os.environ["PATH"] = path_with_front(sys.argv[2:], old)
    os.execv(exe, [exe])
    return 127


if __name__ == "__main__":
    raise SystemExit(main())
