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
  CPMD_CONVERGENCE orbital convergence (CONVERGENCE ORBITALS), default 1e-6
  CPMD_MAXITER    wavefunction optimisation steps, default 300

Units: eOn passes Angstrom and the deck says ANGSTROM, so CPMD converts with
its own Bohr (cnst.mod: 0.529177210859 A). CPMD answers in Hartree and
Hartree/Bohr; the script converts with CODATA 2018 (27.211386245988 eV,
0.529177210903 A), as rgpot and cpmdc do for the in-process route. The two
Bohr values differ by 8.3e-11 relative.

The wavefunction of the previous call stays in RESTART.1 in the exchange
directory, so every call after the first starts from it. With
CPMD_RESTART_POOL set to a shared directory, each converged call also
publishes its RESTART.1 there, keyed by species order, cell, cutoff,
functional and pseudopotentials, and a job's first call starts from the
newest compatible one instead of an atomic guess. CPMDC_RESTART names the
cpmdc-restart tool (OmniPotentRPC/cpmdc); when set, a pooled file is used
only if `cpmdc-restart info` reports the same atom count, species counts
and cutoff.
"""
from __future__ import annotations

import os
import re
import hashlib
import json
import shlex
import shutil
import tempfile
import subprocess
import sys
from pathlib import Path

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903


def read_input(path: Path):
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    cell = [[float(v) for v in rows[i][:3]] for i in range(3)]
    atoms = [(int(r[0]), float(r[1]), float(r[2]), float(r[3])) for r in rows[3:]]
    # eOn sends a zero box for a system that is not periodic. CPMD needs a
    # cell, and an invented one would change the energy without a word.
    if all(abs(v) < 1e-12 for row in cell for v in row):
        raise RuntimeError("from_eon_to_extpot has a zero cell; CPMD needs a "
                           "periodic box (give the system a cell in pos.con)")
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


def write_input(path: Path, cell, atoms, pps, restart: bool, pcg: bool = False):
    order = species_order(atoms)
    conv = float(os.environ.get("CPMD_CONVERGENCE", "1e-6"))
    lines = ["&CPMD", " OPTIMIZE WAVEFUNCTION", " CONVERGENCE ORBITALS", f"  {conv:.3E}",
             " MAXITER", f"  {int(os.environ.get('CPMD_MAXITER', '300'))}",
             " PRINT FORCES ON", " STORE", "  1"]
    if pcg:
        lines.append(" PCG MINIMIZE")
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


def pool_dir(cell, atoms, pps):
    root = os.environ.get("CPMD_RESTART_POOL")
    if not root:
        return None
    species = list(dict.fromkeys(a[0] for a in atoms))
    key = json.dumps({
        "species": [(z, sum(1 for a in atoms if a[0] == z), pps[z]) for z in species],
        "cell": [[round(v, 6) for v in row] for row in cell],
        "cutoff": float(os.environ.get("CPMD_CUTOFF", "30")),
        "functional": os.environ.get("CPMD_FUNCTIONAL", "BLYP"),
    }, sort_keys=True)
    return Path(root) / hashlib.sha256(key.encode()).hexdigest()[:16]


def restart_matches(path: Path, atoms) -> bool:
    """Check a pooled RESTART.1 with cpmdc-restart info, when available."""
    tool = os.environ.get("CPMDC_RESTART")
    if not tool:
        return True
    out = subprocess.run([tool, "info", str(path)], capture_output=True, text=True)
    if out.returncode != 0:
        return False
    info = out.stdout
    m = re.search(r"coordinates:\s+(\d+)", info)
    counts = [int(c) for c in re.findall(r"na\[\d+\]=(\d+)", info)]
    want = [sum(1 for a in atoms if a[0] == z) for z in dict.fromkeys(a[0] for a in atoms)]
    cut = re.search(r"ecut=([0-9.eE+-]+)", info)
    return (m is not None and int(m.group(1)) == len(atoms) and counts == want
            and cut is not None
            and abs(float(cut.group(1)) - float(os.environ.get("CPMD_CUTOFF", "30"))) < 1e-9)


def seed_from_pool(work: Path, pool, atoms) -> bool:
    if pool is None or not (pool / "RESTART.1").exists():
        return False
    if not restart_matches(pool / "RESTART.1", atoms):
        return False
    shutil.copyfile(pool / "RESTART.1", work / "RESTART.1")
    (work / "LATEST").write_text("./RESTART.1\n           1\n")
    return True


def publish_to_pool(work: Path, pool) -> None:
    if pool is None or not (work / "RESTART.1").exists():
        return
    pool.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=pool, prefix=".RESTART.")
    os.close(fd)
    shutil.copyfile(work / "RESTART.1", tmp)
    os.replace(tmp, pool / "RESTART.1")


def main() -> int:
    work = Path.cwd()
    cell, atoms = read_input(work / "from_eon_to_extpot")
    pps = pp_table(os.environ["CPMD_PP"])
    pool = pool_dir(cell, atoms, pps)
    restart = (work / "RESTART.1").exists() and (work / "LATEST").exists()
    seeded = False
    if not restart:
        seeded = restart = seed_from_pool(work, pool, atoms)
    launch = shlex.split(os.environ.get("CPMD_LAUNCH", "mpirun -np 4"))
    cmd = launch + [os.environ["CPMD_BIN"], "cpmd.inp", os.environ["PP_LIBRARY_PATH"]]

    def attempt(restart: bool, pcg: bool):
        order = write_input(work / "cpmd.inp", cell, atoms, pps, restart, pcg)
        # A GEOMETRY left by the previous call must not pass for this one.
        (work / "GEOMETRY").unlink(missing_ok=True)
        with open(work / "cpmd.out", "w") as out:
            rc = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, cwd=work).returncode
        text = (work / "cpmd.out").read_text(errors="replace")
        if rc != 0:
            sys.stderr.write(text[-4000:])
            raise RuntimeError(f"cpmd.x exited with status {rc}")
        return parse_output(text, work, len(atoms), order)

    try:
        energy, forces = attempt(restart, pcg=False)
    except RuntimeError as err:
        if "did not converge" not in str(err):
            raise
        # A wavefunction carried over from a distant geometry can stall the
        # SCF. Retry once from an atomic guess with PCG before giving up.
        sys.stderr.write(f"{err}; retrying cold with PCG MINIMIZE\n")
        (work / "RESTART.1").unlink(missing_ok=True)
        (work / "LATEST").unlink(missing_ok=True)
        seeded = False
        energy, forces = attempt(False, pcg=True)
        (work / "retried_cold").write_text("1\n")
    publish_to_pool(work, pool)
    if seeded:
        (work / "seeded_from_pool").write_text(str(pool) + "\n")
    with open(work / "from_extpot_to_eon", "w") as f:
        f.write(f"{energy:.12f}\n")
        for fx, fy, fz in forces:
            f.write(f"{fx:.12f} {fy:.12f} {fz:.12f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
