"""Launch test_relax_engine with the DLL search order the link used.

Meson prepends the directories of every linked DLL to PATH. That puts a
wheel bin ahead of the prefix, and the loader then binds a same-named
DLL that does not export the entry point (0xC0000139). This process
rebuilds PATH as the executable directory, the staged search list, then
whatever was inherited, and reports a missing import before the test runs.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from win_stage_imports import audit_missing


def _unique(parts: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for part in parts:
        key = part.rstrip("\\/").lower()
        if not part or key in seen:
            continue
        seen.add(key)
        ordered.append(part)
    return ordered


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: win_relax_launch.py EXE SEARCH_DIR...", file=sys.stderr)
        return 2
    exe = Path(argv[1])
    search = [Path(p) for p in argv[2:] if Path(p).is_dir()]
    front = _unique([str(exe.parent)] + [str(p) for p in search])
    inherited = os.environ.get("PATH", "").split(os.pathsep)
    os.environ["PATH"] = os.pathsep.join(_unique(front + inherited))
    problems = audit_missing(exe, [exe.parent] + search)
    if problems:
        print(f"win_relax_launch: {len(problems)} missing imports", file=sys.stderr)
        for line in problems[:40]:
            print(f"win_relax_launch: {line}", file=sys.stderr)
    completed = subprocess.run([str(exe)], check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
