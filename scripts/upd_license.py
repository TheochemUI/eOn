from pathlib import Path
import re

directory = Path(".")

file_extensions = {".h", ".hpp", ".c", ".cpp"}

new_license = """/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
"""

def find_license(top_lines):
    _eonpat = re.compile(r"(.*eOn is free)")
    _gnupat = re.compile(r"(GNU General Public)")
    return _eonpat.search(top_lines) and _gnupat.search(top_lines)

def has_any_license(top_lines):
    license_patterns = [
        re.compile(r"(SPDX-License-Identifier:)"),
        re.compile(r"(GNU General Public License)"),
        re.compile(r"(BSD 3-Clause License)"),
    ]
    return any(pat.search(top_lines) for pat in license_patterns)

def replace_license(file_path):
    with file_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()

    top_lines = '\n'.join(lines[:10])

    if find_license(top_lines):
        # Replace the first 10 lines with the new license
        new_content = new_license + "".join(lines[10:])
        print(f"Replaced license in {file_path}")
    elif not has_any_license(top_lines):
        new_content = new_license + "".join(lines)
        print(f"Prepended license in {file_path}")
    else:
        new_content = "".join(lines)
        print(f"No match found in {file_path}")

    with file_path.open("w", encoding="utf-8") as f:
        f.write(new_content)

ignore_if = ["Asap", "mcamc", ".venv", "thirdparty", "tools"]
for file_path in directory.rglob("*"):
    if file_path.suffix in file_extensions:
        if not any(x in str(file_path) for x in ignore_if):
            replace_license(file_path)

print("License replacement completed.")
