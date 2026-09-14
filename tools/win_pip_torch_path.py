"""Emit torch/metatensor native lib dirs for Windows LoadLibrary.

Git Bash PATH must use colon-separated /c/... paths. Also call
os.add_dll_directory for the current process when --register is passed.
"""
from __future__ import annotations

import argparse
import os
import pathlib
import site
import sys


def native_lib_dirs() -> list[pathlib.Path]:
    roots: list[pathlib.Path] = []
    for sp in site.getsitepackages():
        base = pathlib.Path(sp)
        candidates = [
            base / "torch" / "lib",
            base / "torch" / "bin",
            base / "vesin" / "lib",
            base / "metatensor" / "lib",
        ]
        for pat in (
            "metatensor_torch/torch-*/lib",
            "metatomic/torch/torch-*/lib",
        ):
            candidates.extend(base.glob(pat))
        for d in candidates:
            if d.is_dir():
                roots.append(d.resolve())
    # stable unique
    out: list[pathlib.Path] = []
    seen: set[str] = set()
    for d in roots:
        key = str(d).lower()
        if key not in seen:
            seen.add(key)
            out.append(d)
    return out


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--register",
        action="store_true",
        help="os.add_dll_directory for this process and prepend PATH",
    )
    p.add_argument(
        "--bash-path",
        action="store_true",
        help="print colon-separated Git-Bash PATH prefix (cygpath if available)",
    )
    p.add_argument(
        "--import-check",
        action="store_true",
        help="register dirs then import vesin/metatensor/metatomic torch stack",
    )
    args = p.parse_args()
    dirs = native_lib_dirs()
    if not dirs:
        print("no native lib dirs found under site-packages", file=sys.stderr)
        return 1

    if args.register or args.import_check:
        for d in dirs:
            if hasattr(os, "add_dll_directory"):
                os.add_dll_directory(str(d))
            os.environ["PATH"] = str(d) + os.pathsep + os.environ.get("PATH", "")

    if args.bash_path:
        parts: list[str] = []
        for d in dirs:
            s = str(d)
            # Prefer cygpath for Git Bash; fall back to /d/ form
            try:
                import subprocess

                u = subprocess.check_output(
                    ["cygpath", "-u", s], text=True
                ).strip()
                parts.append(u)
            except Exception:
                # D:\foo\bar -> /d/foo/bar
                if len(s) >= 2 and s[1] == ":":
                    parts.append("/" + s[0].lower() + s[2:].replace("\\", "/"))
                else:
                    parts.append(s.replace("\\", "/"))
        print(":".join(parts))
        return 0

    if args.import_check:
        import metatensor  # noqa: F401
        import metatensor.torch  # noqa: F401
        import metatomic.torch  # noqa: F401
        import vesin  # noqa: F401

        print("pip_metatomic imports ok")
        return 0

    for d in dirs:
        print(d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
