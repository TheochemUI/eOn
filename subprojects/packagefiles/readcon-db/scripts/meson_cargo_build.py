#!/usr/bin/env python3
"""Build the readcon-db C ABI via cargo rustc and copy artifacts.

Used by meson.build. Headers are shipped; this script never invokes cbindgen.
On MSVC the cdylib's import library is ``<dll>.dll.lib`` (often only under
``deps/``). The static archive is ``readcon_db.lib`` and is not that file.
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


def locate_implib(built: Path, shared_name: str) -> Path | None:
    """Return rustc's MSVC import library, or None when it was not emitted."""
    exact = f"{shared_name}.lib"
    for cand in (built / exact, built / "deps" / exact):
        if cand.is_file() and cand.stat().st_size > 0:
            return cand
    matches = [p for p in built.rglob(exact) if p.is_file() and p.stat().st_size > 0]
    if len(matches) == 1:
        return matches[0]
    for dll in built.rglob("*.dll"):
        if "readcon_db" not in dll.name or not dll.is_file():
            continue
        sibling = dll.parent / f"{dll.name}.lib"
        if sibling.is_file() and sibling.stat().st_size > 0:
            return sibling
    return None


def export_names(text: str) -> list[str]:
    """Names from ``dumpbin /exports``. The table ends at the Summary."""
    names: list[str] = []
    in_table = False
    for line in text.splitlines():
        lower = line.lower()
        if "ordinal" in lower and "hint" in lower and "name" in lower:
            in_table = True
            continue
        if not in_table:
            continue
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.lower().startswith("summary"):
            break
        parts = stripped.split()
        if len(parts) < 4 or not parts[0].isdigit():
            continue
        name = parts[-1]
        if name.lower() in {"[noname]", "noname"}:
            continue
        names.append(name.split("=", 1)[0])
    return names


def _msvc_machine() -> str:
    arch = os.environ.get("VSCMD_ARG_TGT_ARCH", "").lower()
    named = {"x64": "x64", "x86": "X86", "arm64": "ARM64"}
    if arch in named:
        return named[arch]
    machine = platform.machine().lower()
    if machine in {"amd64", "x86_64"}:
        return "x64"
    if machine in {"x86", "i386"}:
        return "X86"
    if machine in {"arm64", "aarch64"}:
        return "ARM64"
    return "x64"


def write_implib(dll: Path, out_lib: Path) -> None:
    """Write an MSVC import library from the DLL export table."""
    dumpbin = shutil.which("dumpbin") or shutil.which("dumpbin.exe")
    lib_exe = shutil.which("lib") or shutil.which("lib.exe")
    if dumpbin is None or lib_exe is None:
        sys.exit(
            "MSVC import library was not produced and dumpbin/lib are not on PATH"
        )
    text = subprocess.check_output(
        [dumpbin, "/nologo", "/exports", str(dll)],
        text=True,
        errors="replace",
    )
    names = export_names(text)
    if not names:
        sys.exit(f"no exports in {dll}\n{text[:2000]}")
    out_lib.parent.mkdir(parents=True, exist_ok=True)
    out_def = out_lib.with_suffix(".def")
    out_def.write_text(
        f"LIBRARY {dll.name}\nEXPORTS\n" + "\n".join(names) + "\n",
        encoding="ascii",
    )
    subprocess.check_call(
        [
            lib_exe,
            "/nologo",
            f"/machine:{_msvc_machine()}",
            f"/def:{out_def}",
            f"/out:{out_lib}",
        ]
    )
    if not out_lib.is_file() or out_lib.stat().st_size == 0:
        sys.exit(f"lib.exe did not write {out_lib}")


def copy_or_make_implib(
    built: Path, shared_name: str, out_shared: Path, out_implib: Path
) -> None:
    if out_implib.is_file() and out_implib.stat().st_size > 0:
        return
    found = locate_implib(built, shared_name)
    if found is not None:
        if found.resolve() != out_implib.resolve():
            shutil.copy2(found, out_implib)
        return
    if not out_shared.is_file():
        sys.exit(f"cannot make an import library; DLL is missing: {out_shared}")
    write_implib(out_shared, out_implib)


def _parse(argv: list[str]) -> tuple[list[str], str, str]:
    if len(argv) not in (9, 10, 11):
        print(
            "usage: meson_cargo_build.py CARGO SRC_ROOT TARGET_DIR PROFILE "
            "SHARED_NAME STATIC_NAME OUT_SHARED OUT_STATIC [OUT_IMPLIB] [FEATURES]",
            file=sys.stderr,
        )
        sys.exit(2)
    base = argv[1:9]
    rest = argv[9:]
    out_implib = ""
    features = ""
    if len(rest) == 2:
        out_implib, features = rest
    elif len(rest) == 1:
        if rest[0].endswith(".lib"):
            out_implib = rest[0]
        else:
            features = rest[0]
    return base, out_implib, features


def main(argv: list[str]) -> int:
    (
        cargo,
        src_root,
        target_dir,
        profile,
        shared_name,
        static_name,
        out_shared,
        out_static,
    ), out_implib, features = _parse(argv)
    src_root_p = Path(src_root)
    cmd = [
        cargo,
        "rustc",
        "--lib",
        "--manifest-path",
        str(src_root_p / "Cargo.toml"),
        "--target-dir",
        target_dir,
    ]
    if (src_root_p / "Cargo.lock").is_file():
        cmd.append("--locked")
    if (src_root_p / "vendor").is_dir():
        cmd.append("--offline")
    if profile == "release":
        cmd.append("--release")
    if features.strip():
        cmd.extend(["--features", features.strip()])
    rustc_args = ["--crate-type=cdylib", "--crate-type=staticlib"]
    if sys.platform.startswith("linux"):
        rustc_args.append(f"-Clink-arg=-Wl,-soname,{shared_name}")
    elif sys.platform == "darwin":
        rustc_args.append(f"-Clink-arg=-Wl,-install_name,@rpath/{shared_name}")
    cmd.extend(["--", *rustc_args])
    env = os.environ.copy()
    subprocess.check_call(cmd, cwd=src_root, env=env)
    built = Path(target_dir) / profile
    shutil.copy2(built / shared_name, out_shared)
    shutil.copy2(built / static_name, out_static)
    if out_implib:
        copy_or_make_implib(built, shared_name, Path(out_shared), Path(out_implib))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
