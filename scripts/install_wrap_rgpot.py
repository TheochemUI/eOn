#!/usr/bin/env python3
"""Copy wrap-built librgpot into the install prefix.

rgpot sets install: not meson.is_subproject(). eOn CI uses
meson install --skip-subprojects. eonclient then looks for
librgpot.<kCxxSoversion> (currently 4) next to the prefix.
The SONAME link is taken from the regular file. Following the
unversioned symlink wrote librgpot.so.4 back at librgpot.so and
meson install stopped on that cycle.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def elf_soname(path: Path) -> str | None:
    """DT_SONAME of an ELF shared object.

    rgpot 3.5.1 writes librgpot.so.3.5.1 and records librgpot.so.4.
    The filename major is the project version, not the SONAME.
    """
    readelf = shutil.which("readelf")
    if readelf is None:
        return None
    try:
        out = subprocess.check_output(
            [readelf, "-d", str(path)],
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    for line in out.splitlines():
        if "SONAME" in line and "[" in line:
            return line.split("[", 1)[1].split("]", 1)[0]
    return None


def filename_soname(name: str) -> str | None:
    """librgpot.so.3.5.1 -> librgpot.so.3 when readelf is unavailable."""
    if ".so." not in name:
        return None
    head, rest = name.split(".so.", 1)
    major = rest.split(".", 1)[0]
    if not major:
        return None
    return head + ".so." + major


def place_soname_alias(dest_dir: Path, soname: str, filename: str, src: Path) -> None:
    """Link ``soname`` to ``filename``.

    The alias comes from the regular file's DT_SONAME. A link already at
    that name stays when it names ``filename``, and is replaced when it
    names another path.
    """
    alias = dest_dir / soname
    if alias.name == filename:
        return
    if alias.is_symlink() and os.readlink(alias) == filename:
        return
    if os.path.lexists(alias) or alias.is_symlink():
        alias.unlink()
    try:
        alias.symlink_to(filename)
    except OSError:
        shutil.copy2(src, alias)


def main() -> int:
    src = Path(sys.argv[1])
    prefix = Path(os.environ["MESON_INSTALL_DESTDIR_PREFIX"])
    libdir_arg = Path(sys.argv[2])
    dest_dirs = [libdir_arg if libdir_arg.is_absolute() else prefix / libdir_arg]
    # meson-python packs the extension dir, not prefix/lib. Drop a copy
    # next to _core*.so so the wheel carries the SONAME.
    dest_dirs.extend(
        core.parent
        for core in prefix.rglob("_core*.so")
        if core.is_file()
    )
    dest_dirs.extend(
        core.parent
        for core in prefix.rglob("_core*.pyd")
        if core.is_file()
    )
    # meson library.full_path() can be the .dylib.p object dir.
    parent = src.parent if src.is_dir() else src.parent
    copied = 0
    for dest_dir in dest_dirs:
        dest_dir.mkdir(parents=True, exist_ok=True)
        for path in parent.iterdir():
            name = path.name
            if path.is_dir():
                continue
            if not (
                name.startswith("librgpot")
                or name.startswith("rgpot.")
                or name == "rgpot.dll"
                or name == "rgpot.lib"
            ):
                continue
            dest = dest_dir / name
            if path.is_symlink():
                # os.path.lexists sees a dangling link and a link cycle.
                # Path.exists follows the cycle and reports the path missing.
                if os.path.lexists(dest) or dest.is_symlink():
                    dest.unlink()
                dest.symlink_to(os.readlink(path))
            else:
                shutil.copy2(path, dest)
                soname = elf_soname(dest) or filename_soname(name)
                if soname:
                    place_soname_alias(dest_dir, soname, name, path)
            copied += 1
    if copied == 0:
        print(f"install_wrap_rgpot: no librgpot next to {src}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
