"""Pack and inspect eOn LAMMPS run-input bundles (.eonlpb).

The bundle is a single self-contained file that holds in.lammps plus
every file the run needs (pair_coeff data, custom pair_style .so
plugins, KIM tables, read_data inputs, shell helpers, ...). The
client-side reader in client/potentials/LAMMPS/LammpsBundle.{h,cpp}
extracts the bundle to a per-instance scratch dir and pins liblammps's
working directory there, so the eonclient CWD no longer matters for
LAMMPS file lookups.

Format (all little-endian):
    [0..7]   magic    : b"EONLPB1\\0"  (8 bytes)
    [8..15]  m_len    : uint64          (manifest length in bytes)
    [16..]   manifest : JSON UTF-8      (m_len bytes)
    [...]    bodies   : concatenated file contents in manifest order

Manifest schema:
    {"files": [{"name": "in.lammps", "size": 1234},
               {"name": "Pd.eam.alloy", "size": 56789}]}

Usage::

    python -m eon.lammps_bundle pack potfiles/ bundle.eonlpb
    python -m eon.lammps_bundle list bundle.eonlpb
"""

from __future__ import annotations

import argparse
import json
import os
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

MAGIC = b"EONLPB1\x00"
HEADER_LEN = len(MAGIC) + 8  # magic + uint64 manifest length


@dataclass(frozen=True)
class BundleEntry:
    name: str
    size: int


def _iter_files(root: Path, exclude: Iterable[str] = ()) -> list[tuple[Path, str]]:
    """Walk ``root`` and yield ``(absolute_path, bundle_relative_name)`` pairs.

    Excludes anything in ``exclude`` (matched on the bundle-relative
    name) plus the bundle output file itself when it lives inside ``root``.
    """
    skip = {Path(p).as_posix() for p in exclude}
    out: list[tuple[Path, str]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root).as_posix()
        if rel in skip:
            continue
        out.append((path, rel))
    return out


def pack(source_dir: str | os.PathLike, bundle_path: str | os.PathLike) -> Path:
    """Pack every regular file under ``source_dir`` into ``bundle_path``.

    Refuses to pack a directory that does not contain ``in.lammps`` --
    the client treats a bundle without it as a configuration error.
    Returns the absolute path of the written bundle.
    """
    src = Path(source_dir).resolve()
    bundle = Path(bundle_path).resolve()
    if not src.is_dir():
        raise FileNotFoundError(f"{src} is not a directory")
    if not (src / "in.lammps").is_file():
        raise FileNotFoundError(
            f"{src}/in.lammps is missing; refusing to pack a bundle without it"
        )

    # Skip the bundle output itself if it lives inside the source dir.
    exclude: list[str] = []
    try:
        rel_self = bundle.relative_to(src)
    except ValueError:
        pass
    else:
        exclude.append(rel_self.as_posix())

    entries: list[tuple[Path, str]] = _iter_files(src, exclude=exclude)
    manifest = {
        "files": [
            {"name": rel, "size": path.stat().st_size} for path, rel in entries
        ]
    }
    manifest_bytes = json.dumps(
        manifest, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")

    bundle.parent.mkdir(parents=True, exist_ok=True)
    with bundle.open("wb") as out:
        out.write(MAGIC)
        out.write(struct.pack("<Q", len(manifest_bytes)))
        out.write(manifest_bytes)
        for path, _rel in entries:
            with path.open("rb") as src_f:
                while True:
                    chunk = src_f.read(64 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
    return bundle


def read_manifest(bundle_path: str | os.PathLike) -> list[BundleEntry]:
    """Parse the magic + manifest from a bundle and return its entry list."""
    with Path(bundle_path).open("rb") as f:
        header = f.read(HEADER_LEN)
        if len(header) != HEADER_LEN or not header.startswith(MAGIC):
            raise ValueError(
                f"{bundle_path}: bad magic (expected EONLPB1, got {header[:8]!r})"
            )
        (m_len,) = struct.unpack("<Q", header[len(MAGIC) :])
        manifest = json.loads(f.read(m_len))
    return [BundleEntry(e["name"], int(e["size"])) for e in manifest["files"]]


def _cmd_pack(args: argparse.Namespace) -> int:
    out = pack(args.source_dir, args.bundle)
    entries = read_manifest(out)
    total = sum(e.size for e in entries)
    print(f"wrote {out} ({len(entries)} files, {total} bytes payload)")
    return 0


def _cmd_list(args: argparse.Namespace) -> int:
    entries = read_manifest(args.bundle)
    width = max((len(e.name) for e in entries), default=0)
    for e in entries:
        print(f"{e.name:<{width}}  {e.size:>12} bytes")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_pack = sub.add_parser(
        "pack", help="pack a directory tree into a .eonlpb bundle"
    )
    p_pack.add_argument("source_dir", help="directory containing in.lammps + deps")
    p_pack.add_argument("bundle", help="output bundle path (e.g. bundle.eonlpb)")
    p_pack.set_defaults(func=_cmd_pack)

    p_list = sub.add_parser("list", help="list entries in a bundle")
    p_list.add_argument("bundle", help="input bundle path")
    p_list.set_defaults(func=_cmd_list)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
