"""``import eon`` must work in a checkout that has not been built.

``eon/version.py`` is written by ``tools/gitversion.py`` at build time (see
``eon/meson.build``) and is gitignored, so it is absent in a fresh clone and
in CI before a meson build. ``eon/__init__`` falls back when it is missing
and re-exports the name, which is what the job modules read. A module that
imports the generated file directly bypasses both and fails there.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

EON_PKG = Path(__file__).resolve().parents[1] / "eon"


def _modules_imported_by(path: Path) -> list[str]:
    names: list[str] = []
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            names.append(node.module)
        elif isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
    return names


def test_only_package_init_touches_the_generated_version_module():
    """Everything else reads `version` off the package, which has a fallback."""
    offenders = [
        str(path.relative_to(EON_PKG.parent))
        for path in sorted(EON_PKG.rglob("*.py"))
        if path.name != "__init__.py"
        and "eon.version" in _modules_imported_by(path)
    ]
    assert not offenders, (
        "these import the generated eon.version directly, so they raise "
        f"ModuleNotFoundError in an unbuilt checkout: {offenders}"
    )


def test_package_exposes_both_version_names():
    eon = pytest.importorskip("eon")
    assert isinstance(eon.__version__, str) and eon.__version__
    # Job modules do `from eon import version`; it is the same string.
    assert eon.version == eon.__version__
