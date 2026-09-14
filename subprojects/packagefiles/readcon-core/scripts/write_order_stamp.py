#!/usr/bin/env python3
"""Write a compileable empty translation unit used as the cargo order edge."""

from __future__ import annotations

import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print("usage: write_order_stamp.py OUT.c", file=sys.stderr)
        return 2
    Path(argv[1]).write_text("typedef int readcon_order_stamp;\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
