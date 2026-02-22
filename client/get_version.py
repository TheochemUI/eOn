#!/usr/bin/env python3

import subprocess
import datetime
import os


def get_semantic_version():
    """Read version from pyproject.toml (single source of truth)."""
    toml_path = os.path.join(os.path.dirname(__file__), "..", "pyproject.toml")
    with open(toml_path) as fid:
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
