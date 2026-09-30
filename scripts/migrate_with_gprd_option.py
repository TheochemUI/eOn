#!/usr/bin/env python3
"""Map a stored boolean ``with_gprd`` onto the feature option.

Meson applies the new default ``auto`` to the boolean stored in
``meson-private/coredata.dat`` and rejects it (``value auto is not
boolean``) before ``meson.build`` runs. ``true`` becomes ``enabled``
and ``false`` becomes ``disabled``. A feature value is left as it is.

``python scripts/migrate_with_gprd_option.py <builddir>`` rewrites the
build directory. Pass ``--reconfigure`` to run ``meson setup
--reconfigure <builddir>`` after the rewrite.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def _meson_python() -> str:
    try:
        import mesonbuild  # noqa: F401

        return sys.executable
    except ImportError:
        pass
    meson = shutil.which("meson")
    if meson is None:
        raise SystemExit("meson is not on PATH and mesonbuild is not importable")
    first = Path(meson).read_text(encoding="utf-8", errors="replace").splitlines()[0]
    if first.startswith("#!"):
        return first[2:].strip().split()[0]
    return sys.executable


def _reexec_under_meson() -> None:
    try:
        import mesonbuild  # noqa: F401

        return
    except ImportError:
        py = _meson_python()
        if os.path.realpath(py) == os.path.realpath(sys.executable):
            raise SystemExit("mesonbuild is not importable")
        os.execv(py, [py, *sys.argv])


def migrate(builddir: Path) -> str:
    from mesonbuild import coredata
    from mesonbuild.options import UserBooleanOption, UserFeatureOption

    core = coredata.load(str(builddir))
    store = core.optstore
    key = None
    for candidate in store.options:
        if candidate.name == "with_gprd" and candidate.subproject == "":
            key = candidate
            break
    if key is None:
        raise SystemExit(f"no project option with_gprd in {builddir}")
    old = store.options[key]
    if isinstance(old, UserFeatureOption):
        return "feature"
    if not isinstance(old, UserBooleanOption):
        raise SystemExit(
            f"with_gprd is {type(old).__name__}, not a boolean or a feature"
        )
    mapped = "enabled" if old.value else "disabled"
    updated = UserFeatureOption(
        old.name,
        old.description,
        mapped,
        yielding=old.yielding,
        deprecated=old.deprecated,
        readonly=old.readonly,
    )
    store.options[key] = updated
    coredata.save(core, str(builddir))
    return mapped


def main(argv: list[str]) -> int:
    _reexec_under_meson()
    reconfigure = False
    args = []
    for arg in argv:
        if arg == "--reconfigure":
            reconfigure = True
        else:
            args.append(arg)
    if len(args) != 1:
        raise SystemExit(
            "usage: migrate_with_gprd_option.py [--reconfigure] <builddir>"
        )
    builddir = Path(args[0]).resolve()
    if not (builddir / "meson-private" / "coredata.dat").is_file():
        raise SystemExit(f"no meson build directory at {builddir}")
    mapped = migrate(builddir)
    if mapped == "feature":
        print(f"{builddir}: with_gprd is already a feature")
    else:
        print(f"{builddir}: with_gprd boolean mapped to {mapped}")
    if reconfigure:
        meson = shutil.which("meson")
        if meson is None:
            raise SystemExit("meson is not on PATH")
        subprocess.check_call(
            [meson, "setup", "--reconfigure", str(builddir)],
            cwd=builddir.parent,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
