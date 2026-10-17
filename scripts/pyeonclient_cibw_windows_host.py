#!/usr/bin/env python3
"""Emit the Windows cibuildwheel host environment.

The wheel job links with MSVC. MinGW's gfortran, ld.exe, and windres.exe,
and Strawberry's pkg-config, make Meson's CMake compiler probe fail before
Eigen is found. A small pixi prefix supplies Eigen, pkg-config, and
flang-new. This script turns that prefix plus the process environment into
``KEY=VALUE`` lines for ``GITHUB_ENV``.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path


def msys_to_windows(path: str, *, enable: bool) -> str:
    """Turn a Git bash ``/d/a/...`` prefix into ``D:\\a\\...``.

    ``enable`` is true on the Windows runner. A normal POSIX path is left
    alone when conversion is off, so ``/usr/bin`` is not read as drive U.
    """

    if not enable:
        return path
    norm = path.replace("\\", "/")
    match = re.match(r"^/([A-Za-z])/(.*)$", norm)
    if match is None:
        return path
    rest = match.group(2).replace("/", "\\")
    return f"{match.group(1).upper()}:\\{rest}"


def is_poison_path(entry: str) -> bool:
    """True for PATH entries that break the MSVC CMake probe.

    MinGW supplies ``ld.exe``, ``windres.exe``, and ``gfortran``. Git's
    ``ucrt64`` directory is that same toolchain under another name.
    Strawberry supplies a Perl pkg-config Meson refuses and an MSYS
    ``gmake`` that rewrites ``/nologo``. Git ``usr\\bin`` is not poison:
    ``sh.exe`` is there. Its GNU ``link.exe`` is moved behind MSVC instead
    of deleting the directory.
    """

    norm = entry.replace("\\", "/").rstrip("/").lower()
    return "mingw" in norm or "strawberry" in norm or "ucrt64" in norm


def _norm_path(entry: str) -> str:
    return entry.replace("\\", "/").rstrip("/").lower()


def is_git_usr_bin(entry: str) -> bool:
    """True for Git's ``usr\\bin``, which holds both ``sh.exe`` and GNU link."""

    return _norm_path(entry).endswith("/git/usr/bin")


def is_msvc_link_dir(entry: str) -> bool:
    """True for the MSVC ``HostX64\\x64`` directory that contains ``link.exe``."""

    norm = _norm_path(entry)
    return "vc/tools/msvc" in norm and norm.endswith("/bin/hostx64/x64")


def place_git_sh_after_msvc(entries: list[str]) -> list[str]:
    """Keep Git ``usr\\bin`` on PATH, immediately after MSVC's linker.

    Git bash prepends ``usr\\bin``. Meson uses the first ``link.exe`` it
    finds and rejects the GNU one. ``sh.exe`` still has to stay reachable
    for the rgpot configure.
    """

    git_at = [i for i, entry in enumerate(entries) if is_git_usr_bin(entry)]
    msvc_at = [i for i, entry in enumerate(entries) if is_msvc_link_dir(entry)]
    if not git_at or not msvc_at:
        return entries
    first_msvc = msvc_at[0]
    if all(i > first_msvc for i in git_at):
        return entries
    git_entries = [entries[i] for i in git_at]
    rest = [entry for i, entry in enumerate(entries) if i not in set(git_at)]
    removed_before = sum(1 for i in git_at if i < first_msvc)
    insert_at = first_msvc - removed_before + 1
    return rest[:insert_at] + git_entries + rest[insert_at:]


def split_path(value: str) -> list[str]:
    """Split a PATH-like string without breaking a drive letter.

    Windows lists use ``;``. A lone ``C:\\...`` path has no separator, and
    splitting it on ``:`` would leave the drive letter behind.
    """

    if not value:
        return []
    if ";" in value:
        parts = value.split(";")
    elif "\\" in value or re.match(r"^[A-Za-z]:", value):
        return [value]
    else:
        parts = value.split(":")
    return [part for part in parts if part]


def cleaned_path(value: str, prefix: Path) -> str:
    """Drop poison entries and keep the pixi host bins in front."""

    extras = [str(prefix / "Library" / "bin"), str(prefix / "bin")]
    entries: list[str] = []
    seen: set[str] = set()
    for entry in extras + split_path(value):
        if is_poison_path(entry):
            continue
        key = entry.replace("\\", "/").rstrip("/").lower()
        if key in seen:
            continue
        seen.add(key)
        entries.append(entry)
    return ";".join(place_git_sh_after_msvc(entries))


def _first_file(candidates: list[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def find_eigen_header(prefix: Path) -> Path | None:
    return _first_file(
        [
            prefix / "Library" / "include" / "eigen3" / "Eigen" / "Core",
            prefix / "Library" / "include" / "Eigen" / "Core",
            prefix / "include" / "eigen3" / "Eigen" / "Core",
            prefix / "include" / "Eigen" / "Core",
        ]
    )


def find_eigen_cmake_dir(prefix: Path) -> Path | None:
    relative = (
        "Library/share/eigen3/cmake",
        "Library/share/eigen3",
        "Library/lib/cmake/Eigen3",
        "lib/cmake/Eigen3",
        "share/eigen3/cmake",
    )
    for rel in relative:
        directory = prefix / rel
        if (directory / "Eigen3Config.cmake").is_file() or (
            directory / "eigen3-config.cmake"
        ).is_file():
            return directory
    return None


def find_tool(prefix: Path, names: tuple[str, ...]) -> Path | None:
    bases = (prefix / "Library" / "bin", prefix / "bin", prefix / "Scripts")
    candidates: list[Path] = []
    for base in bases:
        for name in names:
            candidates.append(base / name)
            candidates.append(base / f"{name}.exe")
    return _first_file(candidates)


def pkgconfig_dirs(prefix: Path) -> list[Path]:
    relative = (
        "Library/lib/pkgconfig",
        "lib/pkgconfig",
        "Library/share/pkgconfig",
        "share/pkgconfig",
    )
    return [prefix / rel for rel in relative if (prefix / rel).is_dir()]


def prepend_list(existing: str, extra: list[str]) -> str:
    entries: list[str] = []
    seen: set[str] = set()
    for entry in extra + split_path(existing):
        key = entry.replace("\\", "/").rstrip("/").lower()
        if not key or key in seen:
            continue
        seen.add(key)
        entries.append(entry)
    return ";".join(entries)


def host_environment(prefix: Path, environ: dict[str, str]) -> dict[str, str]:
    """Variables the Windows wheel build needs from ``prefix``."""

    header = find_eigen_header(prefix)
    if header is None:
        raise SystemExit(f"eigen headers missing under {prefix}")
    pkg = find_tool(prefix, ("pkg-config", "pkgconf"))
    if pkg is None:
        raise SystemExit(f"pkg-config missing under {prefix}")
    flang = find_tool(prefix, ("flang-new", "flang"))
    if flang is None:
        raise SystemExit(f"flang-new missing under {prefix}")
    pc_dirs = pkgconfig_dirs(prefix)
    if not pc_dirs:
        raise SystemExit(f"pkgconfig directory missing under {prefix}")

    values: dict[str, str] = {
        "PATH": cleaned_path(environ.get("PATH", ""), prefix),
        "FC": str(flang),
        "CC": "cl",
        "CXX": "cl",
        "AR": "lib",
        "ARFLAGS": "",
        "CMAKE_GENERATOR": "Ninja",
        "CMAKE_C_COMPILER": "cl",
        "CMAKE_CXX_COMPILER": "cl",
        "CMAKE_Fortran_COMPILER": str(flang),
        "PKG_CONFIG": str(pkg),
        "PKG_CONFIG_PATH": ";".join(str(path) for path in pc_dirs),
        "CMAKE_PREFIX_PATH": prepend_list(
            environ.get("CMAKE_PREFIX_PATH", ""),
            [str(prefix / "Library"), str(prefix)],
        ),
        "INCLUDE": prepend_list(
            environ.get("INCLUDE", ""),
            # parent of Eigen/ so includes of <Eigen/Core> resolve
            [str(header.parents[1]), str(prefix / "Library" / "include")],
        ),
        # The flang helper's LIB is already in environ and must stay ahead
        # of the prefix lib dir.
        "LIB": prepend_list(
            str(prefix / "Library" / "lib"),
            split_path(environ.get("LIB", "")),
        ),
        "MSYS_NO_PATHCONV": "1",
        "MSYS2_ARG_CONV_EXCL": "*",
    }
    eigen_dir = find_eigen_cmake_dir(prefix)
    if eigen_dir is not None:
        values["Eigen3_DIR"] = str(eigen_dir)
    for key in ("FLANG_RT_LIB_DIR", "FLANG_RT_DLL_DIR", "FLANG_RT_DIR"):
        if environ.get(key):
            values[key] = environ[key]
    return values


def github_env_block(values: dict[str, str]) -> str:
    lines: list[str] = []
    for key, value in values.items():
        if "\n" in key or "\r" in key or "\n" in value or "\r" in value:
            raise SystemExit(f"refusing to emit {key}: value contains a newline")
        lines.append(f"{key}={value}")
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, help="pixi env prefix")
    parser.add_argument(
        "--emit",
        action="store_true",
        help="write KEY=VALUE lines to stdout",
    )
    args = parser.parse_args(argv)
    prefix_text = msys_to_windows(args.prefix, enable=os.name == "nt")
    prefix = Path(prefix_text)
    values = host_environment(prefix, dict(os.environ))
    if args.emit:
        sys.stdout.write(github_env_block(values))
    else:
        for key, value in values.items():
            print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
