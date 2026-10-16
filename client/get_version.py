#!/usr/bin/env python3

import datetime
import argparse
import re
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


def get_git_revision(source_dir=None):
    source_dir = Path(source_dir) if source_dir else Path(__file__).resolve().parent
    try:
        revision = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=source_dir,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
        return revision if re.fullmatch(r"[0-9a-fA-F]{40}", revision) else "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def get_git_hash():
    revision = get_git_revision()
    return revision[:8] if revision != "unknown" else revision


def get_build_date():
    return datetime.datetime.now(datetime.timezone.utc).strftime(
        "%a %b %d %I:%M:%S %p GMT %Y"
    )


def write_header(template_path, output_path, source_dir=None):
    """Refresh revision metadata while preserving unchanged output timestamps."""
    template = Path(template_path).read_text()
    output = Path(output_path)
    revision = get_git_revision(source_dir)
    short_revision = revision[:8] if revision != "unknown" else revision
    rendered = template.replace("@GIT_HASH_FULL@", revision).replace(
        "@GIT_HASH@", short_revision)
    existing = output.read_text() if output.exists() else ""
    date = re.search(r'const std::string BUILD_DATE = "([^"\n]*)";', existing)
    if date and rendered.replace("@BUILD_DATE@", date[1]) == existing:
        return False
    rendered = rendered.replace("@BUILD_DATE@", get_build_date())
    if rendered == existing:
        return False
    output.write_text(rendered)
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--template")
    parser.add_argument("--output")
    args = parser.parse_args()
    if args.template or args.output:
        if not (args.template and args.output):
            parser.error("--template and --output are required together")
        write_header(args.template, args.output)
        return
    version = get_semantic_version()
    build_date = get_build_date()
    git_hash = get_git_hash()
    print(f"{version},{build_date},{git_hash}")


if __name__ == "__main__":
    main()
