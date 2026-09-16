#!/usr/bin/env python3

import datetime
import subprocess
from pathlib import Path


def get_semantic_version():
    """Read version from pyproject.toml (single source of truth)."""
    toml_path = Path(__file__).resolve().parent.parent / "pyproject.toml"
    with toml_path.open() as fid:
        for line in fid:
            if line.startswith("version ="):
                return line.strip().split(" = ")[1].strip('"').strip("'")
    return "unknown"


def get_git_hash():
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
            )
            .decode()
            .strip()
        )
    except Exception:
        return "unknown"


def get_build_date():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%a %b %d %I:%M:%S %p GMT %Y"
    )


def main():
    version = get_semantic_version()
    build_date = get_build_date()
    git_hash = get_git_hash()
    print(f"{version},{build_date},{git_hash}")


if __name__ == "__main__":
    main()
