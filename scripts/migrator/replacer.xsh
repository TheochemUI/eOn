#!/usr/bin/env xonsh

import sys
from pathlib import Path

import tomli
from loguru import logger

logger.remove(0)
fmt = "{message}"
logger.add(sys.stdout, format=fmt)

def replacer(tomlname=Path(__file__).resolve().parent / "substitutions.toml"):
    logger.info(f"Loading substitutions from {tomlname}")

    substitutions = tomli.load(tomlname.open('rb'))['substitutions']

    for FROM, TO in substitutions.items():
        logger.info(f"Processing substitution: {FROM} -> {TO}")
        files = $(ag -l @(FROM)).splitlines()
        actionable_files = [x for x in files if 'migrator' not in x]

        ncounter = 0
        for f in actionable_files:
            ncounter +=1
            if Path(f).exists():
                logger.info(f"Replacing in file: {f}")
                sd @(FROM) @(TO) @(f)
            else:
                logger.warning(f"File not found: {f}")

        logger.info(f"Done processing {ncounter} files for {FROM} -> {TO}")


replacer()
