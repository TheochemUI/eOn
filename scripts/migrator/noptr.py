#!/usr/bin/env python3

import subprocess
from pathlib import Path

def replace_make_shared():
    FROM_MAKE = 'std::shared_ptr<Matter>'
    TO_MAKE = 'Matter*'

    result = subprocess.run(['ag', '-l', FROM_MAKE],
                            capture_output=True, text=True)
    files = result.stdout.splitlines()
    actionable_files = [x for x in files if 'migrator' not in x]

    for f in actionable_files:
        file_path = Path(f)
        if file_path.exists():
            subprocess.run(['sd', FROM_MAKE, TO_MAKE,
                            str(file_path)], check=True)
        else:
            print(f"File not found: {file_path}")

replace_make_shared()
