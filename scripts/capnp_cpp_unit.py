"""Republish a Cap'n Proto C++ unit so MSVC cl will compile it.

``capnp compile -oc++`` writes ``<schema>.c++``. cl warning D9024 treats
that extension as an object, D9027 ignores the file, and the librarian
then fails with LNK1181 because the object was never produced.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def compile_unit(
    capnp: str,
    import_path: str,
    src_prefix: str,
    schema: str,
    header: str,
    source: str,
) -> None:
    header_path = Path(header)
    source_path = Path(source)
    if header_path.parent != source_path.parent:
        raise SystemExit("header and source must share a directory")
    header_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            capnp,
            "compile",
            f"-oc++:{header_path.parent}",
            f"--import-path={import_path}",
            f"--src-prefix={src_prefix}",
            schema,
        ],
        check=True,
    )
    base = Path(schema).name
    produced_h = header_path.parent / f"{base}.h"
    produced_c = header_path.parent / f"{base}.c++"
    if not produced_h.is_file():
        raise SystemExit(f"capnp did not write {produced_h}")
    if not produced_c.is_file():
        raise SystemExit(f"capnp did not write {produced_c}")
    if produced_h.resolve() != header_path.resolve():
        os.replace(produced_h, header_path)
    os.replace(produced_c, source_path)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capnp", required=True)
    parser.add_argument("--import-path", required=True)
    parser.add_argument("--src-prefix", required=True)
    parser.add_argument("--schema", required=True)
    parser.add_argument("--header", required=True)
    parser.add_argument("--source", required=True)
    args = parser.parse_args(argv)
    compile_unit(
        args.capnp,
        args.import_path,
        args.src_prefix,
        args.schema,
        args.header,
        args.source,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
