"""Check the harmonic example's exact energy and maximum force."""
import math
from pathlib import Path
import sys

path = Path(sys.argv[1])
values = {}
for line in path.read_text().splitlines():
    row = line.split(maxsplit=1)
    if len(row) == 2:
        values[row[1]] = row[0]
for key, expected in [("potential_energy", 0.5), ("Max_Force", 1.0)]:
    value = float(values[key])
    if not math.isfinite(value) or abs(value - expected) > 1e-12:
        raise SystemExit(f"{key}: got {value}, expected {expected}; inspect {path.parent}")
print("Harmonic point verified: energy 0.5 eV; maximum force 1 eV/Angstrom")
