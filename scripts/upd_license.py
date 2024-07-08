from pathlib import Path
import re

directory = Path(".")

file_extensions = {".h", ".hpp", ".c", ".cpp"}

new_license = """/*
** This file is part of eON.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eON Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eON
*/
"""


def find_license(top_lines):
    _eonpat = re.compile(r"(.*eOn is free)")
    _gnupat = re.compile(r"(GNU General Public)")
    if _eonpat.search(top_lines) and _gnupat.search(top_lines):
        return True
    else:
        return False


def replace_license(file_path):
    with file_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()

    if find_license('\n'.join(lines[:10])):
        # Replace the first 10 lines with the new license
        new_content = new_license + "".join(lines[9:])

        # Write the new content back to the file
        with file_path.open("w") as f:
            f.write(new_content)
        print(f"Replaced license in {file_path}")
    else:
        print(f"No match found in {file_path}")


ignore_if = ["Asap", "mcamc", ".venv"]
for file_path in directory.rglob("*"):
    if file_path.suffix in file_extensions:
        if not any(x in str(file_path) for x in ignore_if):
            replace_license(file_path)

print("License replacement completed.")
