#!/usr/bin/env python3
"""Fail if a Highway kernel in libeonclib contains a fused multiply-add.

Scans objdump of the built client library. A function is in scope when its
demangled name contains a kernel (MaxFreeAtomForceNorm, wrapMinimumImage,
wrapLegacyUnit) and an ``N_`` Highway namespace. x86 mnemonics vfmadd,
vfmsub, vfnmadd and vfnmsub fail the check, as do aarch64 fmla, fmls,
fmadd and fmsub.

Exit 77 (meson skip) when objdump is missing, the platform is Windows, or
the library has no Highway kernel symbols. The allow list starts empty.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

# Substrings of demangled names. Fusion in a listed function is ignored.
# Empty until a kernel proves the fused value cannot reach an output.
ALLOW_FUSION_SUBSTRINGS: tuple[str, ...] = ()

KERNELS = (
    "MaxFreeAtomForceNorm",
    "wrapMinimumImage",
    "wrapLegacyUnit",
)

# Meson treats 77 as a skipped test.
SKIP = 77


def skip(message: str) -> None:
    print(f"SKIP: {message}")
    raise SystemExit(SKIP)


def in_scope(name: str) -> bool:
    if "N_" not in name:
        return False
    if any(allow and allow in name for allow in ALLOW_FUSION_SUBSTRINGS):
        return False
    return any(kernel in name for kernel in KERNELS)


def function_name(line: str) -> str | None:
    gnu = re.match(r"^[0-9a-fA-F]+\s+<(.+)>:\s*$", line)
    if gnu:
        return gnu.group(1)
    if line[:1].isspace():
        return None
    llvm = re.match(r"^(.+\(.*\)):\s*$", line)
    if llvm:
        return llvm.group(1)
    return None


def mnemonic(line: str) -> str | None:
    match = re.match(r"^\s*[0-9a-fA-F]+:\s+(\S+)", line)
    if not match:
        return None
    return match.group(1).lower()


def is_fused(op: str) -> bool:
    if op.startswith(("vfmadd", "vfmsub", "vfnmadd", "vfnmsub")):
        return True
    base = op.split(".", 1)[0]
    return base in {"fmla", "fmls", "fmadd", "fmsub"}


def main() -> None:
    if sys.platform == "win32":
        skip("objdump kernel check is not run on Windows")
    if len(sys.argv) != 2:
        print("usage: check_fp_contract.py <libeonclib>", file=sys.stderr)
        raise SystemExit(2)
    library = Path(sys.argv[1])
    if not library.is_file():
        print(f"library not found: {library}", file=sys.stderr)
        raise SystemExit(1)
    objdump = shutil.which("objdump")
    if objdump is None:
        skip("objdump is not on PATH")

    # --no-show-raw-insn is GNU and LLVM. Retry without it if this objdump
    # rejects the flag (older Xcode objdump).
    disassembly = ""
    last_err = ""
    for extra in (["--no-show-raw-insn"], []):
        proc = subprocess.run(
            [objdump, "-d", *extra, "-C", str(library)],
            check=False,
            capture_output=True,
            text=True,
        )
        if proc.returncode == 0:
            disassembly = proc.stdout
            break
        last_err = proc.stderr.strip() or proc.stdout.strip()
    if not disassembly:
        print(f"objdump failed: {last_err}", file=sys.stderr)
        raise SystemExit(1)

    current = ""
    scoped = False
    scanned = 0
    failures: list[str] = []
    for line in disassembly.splitlines():
        name = function_name(line)
        if name is not None:
            current = name
            scoped = in_scope(name)
            if scoped:
                scanned += 1
            continue
        if not scoped:
            continue
        op = mnemonic(line)
        if op is not None and is_fused(op):
            failures.append(f"{current}\n  {line.strip()}")

    if scanned == 0:
        skip(
            "no Highway kernel symbols in "
            f"{library.name} (Highway disabled, or names were not demangled)"
        )
    if failures:
        print(
            f"fused multiply-add in {len(failures)} instruction(s) "
            f"across Highway kernels of {library}:"
        )
        print("\n".join(failures))
        raise SystemExit(1)
    print(
        f"fp-contract: no fused multiply-add in {scanned} Highway kernel "
        f"function(s) of {library.name}"
    )


if __name__ == "__main__":
    main()
