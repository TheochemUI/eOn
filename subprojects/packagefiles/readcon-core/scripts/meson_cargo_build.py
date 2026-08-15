#!/usr/bin/env python3
"""Build only the readcon-core static C ABI and copy it to the meson output.

The wrap must expose one custom_target output. A dual shared+static
custom_target indexed in declare_dependency.sources does not create a
ninja edge, so eonclib links before the archive exists.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) not in (7, 8):
        print(
            "usage: meson_cargo_build.py CARGO SRC_ROOT TARGET_DIR PROFILE "
            "STATIC_NAME OUT_STATIC [FEATURES]",
            file=sys.stderr,
        )
        return 2
    cargo, src_root, target_dir, profile, static_name, out_static = argv[1:7]
    features = argv[7] if len(argv) == 8 else ""
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
    cmd.extend(["--", "--crate-type=staticlib"])
    subprocess.check_call(cmd, cwd=src_root, env=os.environ.copy())
    built = Path(target_dir) / profile
    shutil.copy2(built / static_name, out_static)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
