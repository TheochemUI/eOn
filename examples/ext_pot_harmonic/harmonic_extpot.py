#!/usr/bin/env python3
"""Harmonic energy in eV and its negative derivative in eV/Angstrom."""
import math
from pathlib import Path

rows = [line.split() for line in
        Path("from_eon_to_extpot").read_text().splitlines() if line.strip()]
if len(rows) < 4 or any(len(row) != 3 for row in rows[:3]):
    raise ValueError("Expected three box rows followed by at least one atom")
positions = []
for row in rows[3:]:
    if len(row) != 4:
        raise ValueError("Expected atomic_number x y z")
    int(row[0])
    xyz = [float(value) for value in row[1:]]
    if not all(math.isfinite(value) for value in xyz):
        raise ValueError("Non-finite coordinate")
    positions.append(xyz)
k = 1.0  # eV/Angstrom squared
energy = 0.5 * k * sum(x*x for xyz in positions for x in xyz)
forces = [[-k*x for x in xyz] for xyz in positions]
with Path("from_extpot_to_eon").open("w") as output:
    output.write(f"{energy:.15g}\n")
    for force in forces:
        output.write(" ".join(f"{value:.15g}" for value in force) + "\n")
