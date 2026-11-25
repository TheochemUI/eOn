#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["readcon-chemfiles==0.17.1"]
# ///
"""Prepare two NEB endpoints with a shared cell and coordinate origin."""

import json
import math
import os
from pathlib import Path
import sys

import readcon

BOHR_ANGSTROM = 0.529177210903


def read_endpoint(path):
    path = Path(path).resolve()
    is_con = path.suffix.lower() == ".con"
    if not is_con and path.suffix.lower() != ".xyz":
        raise ValueError(f"Expected an XYZ or CON endpoint: {path}")
    frames = readcon.read_con(str(path)) if is_con else readcon.read_chemfiles(str(path))
    requested = int(os.environ.get("XYZ_FRAME", "-1"))
    selected = requested if requested >= 0 else len(frames) + requested
    if not 0 <= selected < len(frames):
        raise ValueError(f"XYZ_FRAME={requested} is outside {len(frames)} frames in {path}")
    frame = frames[selected]
    units = "angstrom" if is_con else os.environ.get("XYZ_UNITS", "angstrom").lower()
    if units not in {"angstrom", "bohr"}:
        raise ValueError("XYZ_UNITS must be angstrom or bohr; CON coordinates are always Angstrom")
    scale = BOHR_ANGSTROM if units == "bohr" else 1.0
    atoms = frame.atoms
    if not atoms:
        raise ValueError(f"Endpoint contains no atoms: {path}")
    for atom in atoms:
        atom.x, atom.y, atom.z = atom.x * scale, atom.y * scale, atom.z * scale
        if not all(math.isfinite(v) for v in (atom.x, atom.y, atom.z)):
            raise ValueError(f"Non-finite endpoint coordinate: {path}")
    cell = [v * scale for v in frame.cell]
    if any(v != 0 for v in cell):
        if not all(math.isfinite(v) and v > 0 for v in cell):
            raise ValueError(f"Cell needs positive finite lengths: {path}")
        if any(not math.isfinite(a) or abs(a - 90) > 1e-8 for a in frame.angles):
            raise ValueError(f"This CPMD example requires an orthorhombic cell: {path}")
    else:
        cell = None
    ids = [a.atom_id for a in atoms]
    if len(set(ids)) != len(ids):
        raise ValueError(f"Endpoint atom IDs must be unique: {path}")
    return atoms, cell, {"source": str(path), "frame": selected,
                         "frames": len(frames), "input_units": units,
                         "input_cell_angstrom": cell, "atom_ids": ids}


def prepare(reactant, product, outdir):
    left, left_cell, left_record = read_endpoint(reactant)
    right, right_cell, right_record = read_endpoint(product)
    if len(left) != len(right):
        raise ValueError("NEB endpoints have different atom counts")
    if [a.symbol for a in left] != [a.symbol for a in right]:
        raise ValueError("NEB endpoint element order differs; provide matching atom identities/order")
    if [a.atom_id for a in left] != [a.atom_id for a in right]:
        raise ValueError("NEB endpoint atom ID order differs; provide matching atom identities/order")
    if [a.fixed for a in left] != [a.fixed for a in right]:
        raise ValueError("NEB endpoint fixed-atom flags differ")
    translation = [0.0] * 3
    explicit = os.environ.get("CELL_ANGSTROM", "").split()
    if explicit:
        cell = [float(v) for v in explicit]
        if len(cell) != 3 or not all(math.isfinite(v) and v > 0 for v in cell):
            raise ValueError("CELL_ANGSTROM needs three positive finite lengths")
        cell_source = "CELL_ANGSTROM"
    elif left_cell is not None or right_cell is not None:
        if left_cell is not None and right_cell is not None:
            if not all(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-8) for a, b in zip(left_cell, right_cell)):
                raise ValueError("NEB endpoint cells differ; provide compatible geometries in one shared cell")
        cell = left_cell if left_cell is not None else right_cell
        cell_source = "endpoint metadata"
        if left_cell is None or right_cell is None:
            print("WARNING: one endpoint has no cell; using the other endpoint's cell for both. Coordinates retain their common input origin.", file=sys.stderr)
    else:
        padding = float(os.environ.get("VACUUM_ANGSTROM", "10"))
        if not math.isfinite(padding) or padding <= 0:
            raise ValueError("VACUUM_ANGSTROM must be positive and finite")
        coordinates = [(a.x, a.y, a.z) for a in left + right]
        low = [min(row[k] for row in coordinates) for k in range(3)]
        high = [max(row[k] for row in coordinates) for k in range(3)]
        cell = [high[k] - low[k] + 2 * padding for k in range(3)]
        translation = [padding - low[k] for k in range(3)]
        cell_source = "shared inferred bounding box"
        for atom in left + right:
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
        print("\n*** WARNING: NEB ENDPOINTS HAVE NO CELL; USING ONE ARTIFICIAL VACUUM BOX. ***", file=sys.stderr)
        print(f"Shared cell lengths [Angstrom]: {cell}; padding per side: {padding}.", file=sys.stderr)
        print(f"Both endpoints receive the SAME translation [Angstrom]: {translation}.", file=sys.stderr)
        print("This is NOT the original bulk cell. Set CELL_ANGSTROM for a periodic material.\n", file=sys.stderr)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for filename, atoms in (("reactant.con", left), ("product.con", right)):
        output = outdir / filename
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite endpoint: {output}")
    for filename, atoms in (("reactant.con", left), ("product.con", right)):
        frame = readcon.ConFrame(cell=cell, angles=[90.0] * 3, atoms=atoms)
        readcon.write_con(str(outdir / filename), [frame], precision=17)
    input_indices = {a.atom_id: i for i, a in enumerate(left)}
    written = readcon.read_con(str(outdir / "reactant.con"))[0]
    record = {"reader": "readcon-chemfiles==0.17.1", "reactant": left_record,
              "product": right_record, "cell_source": cell_source,
              "lengths_angstrom": cell, "angles_degree": [90.0] * 3,
              "translation_angstrom": translation,
              "con_atom_ids": [a.atom_id for a in written.atoms],
              "con_order_to_input_index": [input_indices[a.atom_id] for a in written.atoms]}
    (outdir / "endpoints.source.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("usage: uv run --script prepare_neb.py REACTANT.xyz|con PRODUCT.xyz|con OUTDIR")
    prepare(*sys.argv[1:])
