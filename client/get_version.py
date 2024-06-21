#!/usr/bin/env python3

import subprocess
import platform
import datetime
import os


def get_git_version():
    try:
        version = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL
            )
            .decode()
            .strip()
        )
        return version
    except Exception:
        return None


def get_build_date():
    return datetime.datetime.now(datetime.UTC).strftime("%a %b %d %I:%M:%S %p GMT %Y")


def main():
    version = get_git_version() or "unknown"
    build_date = get_build_date()
    print(f"{version},{build_date}")


if __name__ == "__main__":
    main()
