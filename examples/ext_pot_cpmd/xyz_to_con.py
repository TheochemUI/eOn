#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["readcon-chemfiles==0.17.1"]
# ///
"""Read an XYZ trajectory with readcon and record the selected CON geometry."""

import json
import math
import os
from pathlib import Path
import sys

import readcon

BOHR_ANGSTROM = 0.529177210903


def positive_triplet(values, name):
    values = [float(value) for value in values]
    if len(values) != 3 or not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError(f"{name} needs three positive finite lengths")
    return values


def convert(source, output):
    source, output = Path(source).resolve(), Path(output)
    frames = readcon.read_chemfiles(str(source))
    if not frames:
        raise ValueError("XYZ contains no frames")
    requested = int(os.environ.get("XYZ_FRAME", "-1"))
    selected = requested if requested >= 0 else len(frames) + requested
    if not 0 <= selected < len(frames):
        raise ValueError(f"XYZ_FRAME={requested} is outside {len(frames)} frames")
    frame = frames[selected]
    units = os.environ.get("XYZ_UNITS", "angstrom").lower()
    if units not in {"angstrom", "bohr"}:
        raise ValueError("XYZ_UNITS must be angstrom or bohr")
    scale = BOHR_ANGSTROM if units == "bohr" else 1.0
    atoms = frame.atoms
    if not atoms:
        raise ValueError("Selected XYZ frame contains no atoms")
    for atom in atoms:
        atom.x, atom.y, atom.z = atom.x * scale, atom.y * scale, atom.z * scale
        if not all(math.isfinite(v) for v in (atom.x, atom.y, atom.z)):
            raise ValueError("XYZ coordinates must be finite")
    translation = [0.0, 0.0, 0.0]
    explicit = os.environ.get("CELL_ANGSTROM", "").split()
    angles = list(frame.angles)
    if explicit:
        lengths = positive_triplet(explicit, "CELL_ANGSTROM")
        angles, cell_source = [90.0] * 3, "CELL_ANGSTROM"
    elif any(v != 0 for v in frame.cell):
        lengths = positive_triplet([v * scale for v in frame.cell], "XYZ cell")
        cell_source = "XYZ metadata"
    else:
        padding = float(os.environ.get("VACUUM_ANGSTROM", "10"))
        if not math.isfinite(padding) or padding <= 0:
            raise ValueError("VACUUM_ANGSTROM must be positive and finite")
        coordinates = [(a.x, a.y, a.z) for a in atoms]
        low = [min(row[k] for row in coordinates) for k in range(3)]
        high = [max(row[k] for row in coordinates) for k in range(3)]
        lengths = [high[k] - low[k] + 2 * padding for k in range(3)]
        translation = [padding - low[k] for k in range(3)]
        for atom in atoms:
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
        angles, cell_source = [90.0] * 3, "inferred bounding box"
        print("\n*** WARNING: XYZ HAS NO CELL; USING AN ARTIFICIAL VACUUM BOX. ***", file=sys.stderr)
        print(f"Cell lengths [Angstrom]: {lengths}; padding per side: {padding}.", file=sys.stderr)
        print("Coordinates are translated. This is NOT the original bulk cell.", file=sys.stderr)
        print("Set CELL_ANGSTROM to the original orthorhombic cell for bulk calculations.\n", file=sys.stderr)
    if len(angles) != 3 or any(not math.isfinite(a) or abs(a - 90) > 1e-8 for a in angles):
        raise ValueError("This CPMD example requires an orthorhombic cell; use a validated custom deck for other cells")
    converted = readcon.ConFrame(cell=lengths, angles=[90.0] * 3, atoms=atoms)
    readcon.write_con(str(output), [converted], precision=17)
    written = readcon.read_con(str(output))[0]
    input_indices = {atom.atom_id: i for i, atom in enumerate(atoms)}
    record = {
        "source": str(source), "frame": selected, "requested_frame": requested,
        "frames": len(frames), "input_units": units,
        "reader": "readcon-chemfiles==0.17.1", "cell_source": cell_source,
        "lengths_angstrom": lengths, "angles_degree": [90.0] * 3,
        "translation_angstrom": translation,
        "xyz_atom_ids": [atom.atom_id for atom in atoms],
        "con_atom_ids": [atom.atom_id for atom in written.atoms],
        "con_order_to_xyz_index": [input_indices[atom.atom_id] for atom in written.atoms],
    }
    output.with_suffix(".source.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: uv run --script xyz_to_con.py INPUT.xyz OUTPUT.con")
    convert(*sys.argv[1:])
