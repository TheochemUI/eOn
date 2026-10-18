#!/usr/bin/env python3
"""Write a C translation unit that only orders a link after a cargo target.

usage: meson_order_stamp.py OUT.c NAME
"""

import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: meson_order_stamp.py OUT.c NAME", file=sys.stderr)
        return 2
    Path(argv[1]).write_text(f"typedef int {argv[2]};\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
