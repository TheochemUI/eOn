#!/usr/bin/env python3
"""Assert eOn release surfaces agree and optional CHANGELOG section exists.

Mirrors rgpot's potctl ``release assert`` for eOn's real surfaces:
``pyproject.toml``, ``pixi.toml``, optional ``CHANGELOG.md`` section.

Exit codes:
  0 — all checks pass
  1 — surfaces disagree or required changelog section missing
  2 — usage / IO error
"""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

SEMVER_RE = re.compile(
    r"^(?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)"
    r"(?:-(?P<pre>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)"
    r"(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?"
    r"(?:\+(?P<build>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
)


# PyPI distribution name must not collide with occupied global names (see ci/gha/pypi.ncl).
# `eon` on PyPI is EoN (epidemics); our wheel/sdist is eon-akmc (import package eon/).
BLOCKED_PYPI_NAMES = frozenset({"eon"})


def read_pyproject_dist_name(path: Path) -> str:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    proj = data.get("project")
    if not isinstance(proj, dict) or not isinstance(proj.get("name"), str):
        raise KeyError(f"{path}: project.name missing")
    return proj["name"]

CHANGELOG_SECTION_RE = re.compile(
    r"^## \[(?P<ver>[^\]]+)\]\([^\)]*\) - \d{4}-\d{2}-\d{2}\s*$",
    re.MULTILINE,
)


def load_toml_version(path: Path, key_path: tuple[str, ...] = ("project", "version")) -> str:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    cur: object = data
    for key in key_path:
        if not isinstance(cur, dict) or key not in cur:
            # pixi.toml uses top-level version
            if key_path == ("project", "version") and path.name == "pixi.toml":
                ver = data.get("version")
                if isinstance(ver, str):
                    return ver
            raise KeyError(f"{path}: missing {'/'.join(key_path)}")
        cur = cur[key]
    if not isinstance(cur, str):
        raise TypeError(f"{path}: version is not a string")
    return cur


def read_pixi_version(path: Path) -> str:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    ver = data.get("version")
    if not isinstance(ver, str):
        # workspace / project block fallbacks
        for block in ("workspace", "project"):
            blk = data.get(block)
            if isinstance(blk, dict) and isinstance(blk.get("version"), str):
                return blk["version"]
        raise KeyError(f"{path}: no top-level or workspace/project version")
    return ver


def read_pyproject_version(path: Path) -> str:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    proj = data.get("project")
    if not isinstance(proj, dict) or not isinstance(proj.get("version"), str):
        raise KeyError(f"{path}: project.version missing")
    return proj["version"]


def changelog_has_section(changelog: Path, version: str) -> bool:
    if not changelog.is_file():
        return False
    text = changelog.read_text(encoding="utf-8")
    for m in CHANGELOG_SECTION_RE.finditer(text):
        if m.group("ver") == version:
            return True
    # Also accept a simpler ## [version] heading without date (manual edits)
    simple = re.compile(rf"^## \[{re.escape(version)}\]", re.MULTILINE)
    return bool(simple.search(text))


def collect_versions(root: Path) -> dict[str, str]:
    pyproject = root / "pyproject.toml"
    pixi = root / "pixi.toml"
    versions: dict[str, str] = {}
    versions["pyproject.toml"] = read_pyproject_version(pyproject)
    versions["pixi.toml"] = read_pixi_version(pixi)
    return versions


def assert_lockstep(
    root: Path,
    expected: str | None = None,
    *,
    require_changelog: bool = False,
) -> list[str]:
    """Return list of error strings (empty if ok)."""
    errors: list[str] = []
    try:
        versions = collect_versions(root)
    except (OSError, KeyError, TypeError, tomllib.TOMLDecodeError) as exc:
        return [f"failed to read version surfaces: {exc}"]

    vals = set(versions.values())
    if len(vals) != 1:
        detail = ", ".join(f"{k}={v!r}" for k, v in sorted(versions.items()))
        errors.append(f"version surfaces disagree: {detail}")
    else:
        current = next(iter(vals))
        if not SEMVER_RE.match(current):
            errors.append(f"version {current!r} is not valid semver")
        if expected is not None and current != expected:
            errors.append(
                f"surfaces at {current!r} but expected {expected!r} "
                f"({', '.join(f'{k}={v}' for k, v in sorted(versions.items()))})"
            )

    try:
        dist_name = read_pyproject_dist_name(root / "pyproject.toml")
        if dist_name.lower() in BLOCKED_PYPI_NAMES:
            errors.append(
                f"pyproject.toml project.name {dist_name!r} is blocked on PyPI "
                f"(use eon-akmc; see ci/gha/pypi.ncl — {dist_name!r} is EoN epidemics)"
            )
    except (OSError, KeyError, TypeError, tomllib.TOMLDecodeError) as exc:
        errors.append(f"failed to read project.name: {exc}")

    target = expected or (next(iter(vals)) if len(vals) == 1 else None)
    if require_changelog and target:
        cl = root / "CHANGELOG.md"
        if not changelog_has_section(cl, target):
            errors.append(
                f"CHANGELOG.md missing section for {target!r} "
                f"(towncrier should have written ## [{target}](...))"
            )

    return errors


def emit_release_notes(root: Path, version: str, out: Path | None = None) -> str:
    """Extract CHANGELOG section for *version* as release notes body."""
    cl = root / "CHANGELOG.md"
    if not cl.is_file():
        raise FileNotFoundError(cl)
    text = cl.read_text(encoding="utf-8")
    # Match ## [version] ... until next ## [ or EOF
    pat = re.compile(
        rf"(^## \[{re.escape(version)}\][^\n]*\n)(.*?)(?=^## \[|\Z)",
        re.MULTILINE | re.DOTALL,
    )
    m = pat.search(text)
    if not m:
        # Fallback: curated release-notes.md
        curated = root / "docs" / "source" / "releases" / f"v{version}" / "release-notes.md"
        if curated.is_file():
            body = curated.read_text(encoding="utf-8")
        else:
            body = f"Release v{version}\n\nSee CHANGELOG.md for details.\n"
    else:
        body = m.group(1) + m.group(2)
        body = body.strip() + "\n"

    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(body, encoding="utf-8")
    return body


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assert eOn pyproject/pixi versions lockstep; optional CHANGELOG gate.",
    )
    parser.add_argument(
        "version",
        nargs="?",
        default=None,
        help="Expected semver (default: any, as long as surfaces agree)",
    )
    parser.add_argument(
        "--require-changelog",
        action="store_true",
        help="Require CHANGELOG.md section for the target version",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=REPO_ROOT,
        help="Repository root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--notes",
        type=Path,
        default=None,
        help="If set with version, write release notes body to this path",
    )
    parser.add_argument(
        "--print-version",
        action="store_true",
        help="Print agreed version to stdout on success",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()

    errors = assert_lockstep(
        root,
        args.version,
        require_changelog=args.require_changelog,
    )
    if errors:
        for err in errors:
            print(f"error: {err}", file=sys.stderr)
        return 1

    versions = collect_versions(root)
    ver = args.version or versions["pyproject.toml"]
    if args.notes is not None:
        if not args.version:
            print("error: --notes requires an explicit version argument", file=sys.stderr)
            return 2
        try:
            emit_release_notes(root, args.version, args.notes)
        except OSError as exc:
            print(f"error: could not write notes: {exc}", file=sys.stderr)
            return 2

    if args.print_version:
        print(ver)
    else:
        print(f"ok: lockstep version {ver} (pyproject.toml, pixi.toml)")
        if args.require_changelog:
            print(f"ok: CHANGELOG.md has section for {ver}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
