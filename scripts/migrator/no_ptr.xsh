#!/usr/bin/env xonsh

from pathlib import Path

def replace_make_shared():
    FROM_MAKE = "std::shared_ptr<Matter>"
    TO_MAKE = 'Matter*'
    
    files = $(ag -l @(FROM_MAKE)).splitlines()
    actionable_files = [x for x in files if 'migrator' not in x]

    for f in actionable_files:
        if Path(f).exists():
            sd @(FROM_MAKE) @(TO_MAKE) @(f)
            echo @(f)
        else:
            echo "File not found: @(f)"

replace_make_shared()
