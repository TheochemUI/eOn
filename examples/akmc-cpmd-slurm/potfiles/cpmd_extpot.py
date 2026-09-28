#!/usr/bin/env python3
"""eOn ExtPot wrapper around a multi-rank cpmd.x single point.

eOn writes from_eon_to_extpot (3 cell rows, then "Z x y z" per atom, in
Angstrom) in a per-potential exchange directory and runs this script there;
the script answers with from_extpot_to_eon (energy in eV, then forces in
eV/Angstrom, in eOn's atom order).

Settings come from the environment of the client job:
  CPMD_BIN        cpmd.x
  CPMD_LAUNCH     MPI launcher, default "mpirun -np 4"
  PP_LIBRARY_PATH pseudopotential directory (CPMD reads it too)
  CPMD_PP         "Z:file[:options]" entries separated by commas,
                  e.g. "14:Si_MT_BLYP.psp:LMAX=P"
  CPMD_CUTOFF     plane-wave cutoff in Ry, default 30
  CPMD_FUNCTIONAL default BLYP

The wavefunction of the previous call stays in RESTART.1 in the exchange
directory, so every call after the first starts from it.
"""
from __future__ import annotations

import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903


def read_input(path: Path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    cell = [[float(v) for v in rows[i][:3]] for i in range(3)]
    atoms = [(int(r[0]), float(r[1]), float(r[2]), float(r[3])) for r in rows[3:]]
    return cell, atoms


def pp_table(spec: str):
    table = {}
    for entry in filter(None, (s.strip() for s in spec.split(","))):
        z, fname, *opts = entry.split(":")
        table[int(z)] = (fname, opts[0] if opts else "LMAX=P")
    return table


def species_order(atoms):
    """CPMD lists atoms species by species; return eOn indices in that order."""
    order = []
    for z in dict.fromkeys(a[0] for a in atoms):
        order.extend(i for i, a in enumerate(atoms) if a[0] == z)
    return order


def write_input(path: Path, cell, atoms, pps, restart: bool):
    order = species_order(atoms)
    lines = ["&CPMD", " OPTIMIZE WAVEFUNCTION", " CONVERGENCE ORBITALS", "  1.0D-6",
             " MAXITER", "  300", " PRINT FORCES ON", " STORE", "  1"]
    if restart:
        lines.append(" RESTART WAVEFUNCTION LATEST")
    # CELL VECTORS takes the full cell; CPMD refuses SYMMETRY beside it.
    lines += ["&END", "&SYSTEM", " ANGSTROM", " CELL VECTORS"]
    lines += ["  " + " ".join(f"{v:.10f}" for v in row) for row in cell]
    lines += [" CUTOFF", f"  {float(os.environ.get('CPMD_CUTOFF', '30'))}", "&END",
              "&DFT", f" FUNCTIONAL {os.environ.get('CPMD_FUNCTIONAL', 'BLYP')}", "&END",
              "&ATOMS"]
    for z in dict.fromkeys(a[0] for a in atoms):
        fname, lmax = pps[z]
        idx = [i for i in order if atoms[i][0] == z]
        lines += [f"*{fname} KLEINMAN-BYLANDER", f" {lmax}", f" {len(idx)}"]
        lines += [" " + " ".join(f"{atoms[i][k]:.10f}" for k in (1, 2, 3)) for i in idx]
    lines.append("&END")
    path.write_text("\n".join(lines) + "\n")
    return order


def table_forces(text: str, n: int):
    """The last force table in cpmd.out.

    wrgeo prints fion, the ionic force, under the header
    "GRADIENTS (-FORCES)", in 1PE11.3 (three significant digits). The first
    table, printed before the SCF, is all zeros; the last one is the result.
    """
    block = text.rsplit("GRADIENTS (-FORCES)", 1)
    if len(block) != 2:
        raise RuntimeError("cpmd.x output has no force table (PRINT FORCES ON)")
    rows = []
    for line in block[1].splitlines()[1:]:
        cols = line.split()
        if len(cols) < 8 or not cols[0].isdigit():
            if rows:
                break
            continue
        rows.append([float(c) for c in cols[-3:]])
    if len(rows) != n:
        raise RuntimeError(f"cpmd.x printed {len(rows)} forces for {n} atoms")
    return rows


def geometry_forces(path: Path, n: int):
    """Forces from GEOMETRY: after a wavefunction optimisation, columns 4-6
    hold fion in Ha/Bohr at full precision, in species order."""
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    if len(rows) != n:
        raise RuntimeError(f"{path.name} holds {len(rows)} atoms, expected {n}")
    return [[float(c) for c in r[3:6]] for r in rows]


def parse_output(text: str, work: Path, n: int, order):
    energies = re.findall(r"TOTAL ENERGY =\s+(-?\d+\.\d+)\s+A\.U\.", text)
    if not energies:
        raise RuntimeError("cpmd.x output has no TOTAL ENERGY line")
    if "BUT NO CONVERGENCE" in text:
        raise RuntimeError("cpmd.x wavefunction optimisation did not converge")
    table = table_forces(text, n)
    geo = geometry_forces(work / "GEOMETRY", n)
    for t_row, g_row in zip(table, geo):
        for t, g in zip(t_row, g_row):
            # The table rounds to three significant digits.
            if abs(t - g) > 1e-3 * abs(g) + 1e-10:
                raise RuntimeError("GEOMETRY and the cpmd.out force table disagree")
    forces = [None] * n
    scale = HARTREE_EV / BOHR_ANG
    for pos, i in enumerate(order):
        forces[i] = [g * scale for g in geo[pos]]
    return float(energies[-1]) * HARTREE_EV, forces


def main() -> int:
    work = Path.cwd()
    cell, atoms = read_input(work / "from_eon_to_extpot")
    pps = pp_table(os.environ["CPMD_PP"])
    restart = (work / "RESTART.1").exists() or (work / "LATEST").exists()
    order = write_input(work / "cpmd.inp", cell, atoms, pps, restart)
    # A GEOMETRY left by the previous call must not pass for this one.
    (work / "GEOMETRY").unlink(missing_ok=True)
    launch = shlex.split(os.environ.get("CPMD_LAUNCH", "mpirun -np 4"))
    cmd = launch + [os.environ["CPMD_BIN"], "cpmd.inp", os.environ["PP_LIBRARY_PATH"]]
    with open(work / "cpmd.out", "w") as out:
        rc = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, cwd=work).returncode
    text = (work / "cpmd.out").read_text(errors="replace")
    if rc != 0:
        sys.stderr.write(text[-4000:])
        return rc
    energy, forces = parse_output(text, work, len(atoms), order)
    with open(work / "from_extpot_to_eon", "w") as f:
        f.write(f"{energy:.12f}\n")
        for fx, fy, fz in forces:
            f.write(f"{fx:.12f} {fy:.12f} {fz:.12f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
