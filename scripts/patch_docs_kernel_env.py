#!/usr/bin/env python3
"""Patch the active Jupyter kernel environment for docs builds."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

def main() -> int:
    kernel_spec = Path(sys.prefix) / "share/jupyter/kernels/python3/kernel.json"
    spec = json.loads(kernel_spec.read_text())
    env = spec.setdefault("env", {})

    path_entries = []
    existing_path = env.get("PATH") or os.environ.get("PATH", "")
    if existing_path:
        path_entries.extend(p for p in existing_path.split(":") if p)

    def _dedupe(entries: list[str]) -> str:
        seen: set[str] = set()
        ordered: list[str] = []
        for entry in entries:
            if not entry or entry in seen:
                continue
            seen.add(entry)
            ordered.append(entry)
        return ":".join(ordered)

    env["PATH"] = _dedupe(path_entries)

    kernel_spec.write_text(json.dumps(spec, indent=1))

    print(f"Patched {kernel_spec}")
    print(f"  PATH entries: {len(env['PATH'].split(':'))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
