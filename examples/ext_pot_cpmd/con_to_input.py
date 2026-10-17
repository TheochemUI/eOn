#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["readcon-chemfiles==0.17.1"]
# ///
"""Prepare a CON geometry for an independent potential evaluation."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import readcon


def prepare(source, output, requested=0):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite input: {output}")
    frames = readcon.read_con(str(source))
    selected = requested if requested >= 0 else len(frames)+requested
    if not 0 <= selected < len(frames):
        raise ValueError(f"CON_FRAME={requested} is outside {len(frames)} frames")
    frame = frames[selected]
    if not frame.atoms:
        raise ValueError("CON contains no atoms")
    if not all(math.isfinite(v) and v > 0 for v in frame.cell):
        raise ValueError("CON needs positive finite cell lengths")
    if any(not math.isfinite(a) or abs(a-90) > 1e-8 for a in frame.angles):
        raise ValueError("This CPMD example requires an orthorhombic cell")
    atoms = []
    for atom in frame.atoms:
        if not all(math.isfinite(v) for v in (atom.x, atom.y, atom.z)):
            raise ValueError("CON contains a non-finite coordinate")
        atoms.append(readcon.Atom(atom.symbol, atom.x, atom.y, atom.z,
                                  fixed=atom.fixed, atom_id=atom.atom_id, mass=atom.mass))
    clean = readcon.ConFrame(cell=frame.cell, angles=frame.angles, atoms=atoms)
    readcon.write_con(str(output), [clean], precision=17)
    written = readcon.read_first_frame(str(output))
    assert written.energy is None and not written.has_forces
    assert all(not atom.has_energy for atom in written.atoms)
    assert written.cell == frame.cell and written.angles == frame.angles
    fields = lambda a: (a.symbol, a.atom_id, a.mass, tuple(a.fixed), a.x, a.y, a.z)
    assert [fields(a) for a in written.atoms] == [fields(a) for a in frame.atoms]
    record = dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  frame=selected, requested_frame=requested, frames=len(frames),
                  reader="readcon-chemfiles==0.17.1", input_units="angstrom",
                  input_had_energy=frame.energy is not None, input_had_forces=frame.has_forces,
                  cached_energy_and_forces_removed=True, geometry_mass_constraints_preserved=True,
                  output_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
    output.with_suffix(".source.json").write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(record, indent=2))
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--frame", type=int, default=int(os.environ.get("CON_FRAME", "0")))
    args = parser.parse_args()
    prepare(args.source, args.output, args.frame)
