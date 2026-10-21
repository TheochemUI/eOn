#!/usr/bin/env python3
"""Build the eon-rust-abi staticlib and a C order stamp for Meson.

The linker variable is set only when Meson passes one. MSVC must not: rustc
drives link.exe, and pointing CARGO_TARGET_*_LINKER at cl fails the build
script (LNK2001 memcpy, mainCRTStartup, __CxxFrameHandler3). That is the
same override lode-org/minimage's meson.build sets unconditionally.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 9:
        print(
            "usage: cargo_staticlib.py CARGO MANIFEST TARGET PROFILE "
            "OUT_LIB OUT_STAMP LINKER_KEY LINKER",
            file=sys.stderr,
        )
        return 2
    cargo, manifest, target, profile, out_lib, out_stamp, linker_key, linker = argv[1:9]
    manifest_path = Path(manifest)
    cmd = [
        cargo,
        "rustc",
        "--lib",
        "--manifest-path",
        str(manifest_path),
        "--target-dir",
        target,
    ]
    if (manifest_path.parent / "Cargo.lock").is_file():
        cmd.append("--locked")
    if profile == "release":
        cmd.append("--release")
    rustc_args = ["--crate-type=staticlib", "-C", "link-dead-code"]
    if sys.platform != "win32":
        rustc_args.extend(["-C", "relocation-model=pic"])
    cmd.extend(["--", *rustc_args])
    env = os.environ.copy()
    # Meson drops empty command arguments, so an unset linker is "-".
    if linker not in ("", "-") and linker_key not in ("", "-"):
        env[linker_key] = linker
    subprocess.check_call(cmd, cwd=manifest_path.parent, env=env)
    built = Path(target) / profile
    artifact = "eon_rust_abi.lib" if sys.platform == "win32" else "libeon_rust_abi.a"
    shutil.copy2(built / artifact, out_lib)
    Path(out_stamp).write_text("typedef int eon_rust_abi_order;\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
