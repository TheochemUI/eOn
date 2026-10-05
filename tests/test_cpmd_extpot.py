"""The CPMD ExtPot wrapper against a stand-in cpmd.x.

examples/akmc-cpmd-slurm/potfiles/cpmd_extpot.py turns eOn's exchange file
into a CPMD deck, runs cpmd.x, and converts what CPMD prints back into eOn's
units and atom order. The stand-in below reads the deck the wrapper wrote,
evaluates a harmonic well around the cell centre in CPMD's own units
(Bohr from cnst.mod, Hartree), and writes cpmd.out and GEOMETRY the way
cpmd.x lays them out: species by species, forces in Hartree/Bohr. The
tests check the energy, the force sign, the units, the atom order for
interleaved species, the cell, and the refusals.
"""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "examples" / "akmc-cpmd-slurm" / "potfiles" / "cpmd_extpot.py"
HA_EV = 27.211386245988
A0_2018 = 0.529177210903
A0_CPMD = 0.529177210859
K_HA = 0.01  # Hartree per Bohr^2 per atom
# cpmd.x prints TOTAL ENERGY with 8 decimals (F20.8): the file route's
# energy carries half of 1e-8 Hartree, 1.4e-7 eV, of rounding.
E_ROUND_EV = 0.5e-8 * HA_EV

FAKE_CPMD = textwrap.dedent(
    f"""\
    #!{sys.executable}
    import os, sys
    from pathlib import Path
    A0 = {A0_CPMD!r}
    K = {K_HA!r}
    lines = Path(sys.argv[1]).read_text().splitlines()
    Path("argv").write_text(" ".join(sys.argv[1:]) + "\\n")
    i = lines.index(" CELL VECTORS")
    cell = [[float(v) for v in lines[i + 1 + r].split()] for r in range(3)]
    centre = [sum(cell[r][c] for r in range(3)) / 2 / A0 for c in range(3)]
    atoms, species = [], []
    j = lines.index("&ATOMS") + 1
    while lines[j] != "&END":
        species.append(lines[j][1:].split()[0])
        n = int(lines[j + 2])
        for a in range(n):
            atoms.append((species[-1], [float(v) / A0 for v in lines[j + 3 + a].split()]))
        j += 3 + n
    energy = sum(K * sum((x - c) ** 2 for x, c in zip(r, centre)) for _, r in atoms)
    fion = [[-2 * K * (x - c) for x, c in zip(r, centre)] for _, r in atoms]
    if Path("no_convergence_once").exists():
        Path("no_convergence_once").unlink()
        print(" *** BUT NO CONVERGENCE ***")
    print(" ATOM          COORDINATES            GRADIENTS (-FORCES)")
    for k, (s, r) in enumerate(atoms, 1):
        print(f"{{k:4d}}  {{s[:2]:>2s}}" + "".join(f"{{x:12.6f}}" for x in r) + "   "
              + "".join(f"{{0.0:11.3E}}" for _ in r))
    print(f" (K+E1+L+N+X)           TOTAL ENERGY = {{energy:20.8f}} A.U.")
    print(" ATOM          COORDINATES            GRADIENTS (-FORCES)")
    for k, ((s, r), f) in enumerate(zip(atoms, fion), 1):
        print(f"{{k:4d}}  {{s[:2]:>2s}}" + "".join(f"{{x:12.6f}}" for x in r) + "   "
              + "".join(f"{{v:11.3E}}" for v in f))
    with open("GEOMETRY", "w") as g:
        for (_, r), f in zip(atoms, fion):
            g.write("".join(f"{{x:20.12f}}" for x in r) + "".join(f"{{v:20.12f}}" for v in f) + "\\n")
    """
)

CELL = [[10.0, 0.0, 0.0], [0.0, 11.0, 0.0], [0.5, 0.0, 12.0]]
# N, Si, N, Si, N: species interleaved, so CPMD's order differs from eOn's.
ATOMS = [
    (7, 4.1, 5.2, 6.3),
    (14, 6.0, 5.0, 6.9),
    (7, 5.5, 7.1, 5.6),
    (14, 3.9, 4.0, 7.7),
    (7, 6.6, 6.2, 5.2),
]


def write_exchange(path, cell, atoms):
    rows = ["\t".join(f"{v:.19f}" for v in row) for row in cell]
    rows += [f"{z}\t{x:.19f}\t{y:.19f}\t{w:.19f}" for z, x, y, w in atoms]
    (path / "from_eon_to_extpot").write_text("\n".join(rows) + "\n")


def run_wrapper(workdir, **extra):
    env = dict(os.environ)
    env.update(
        CPMD_BIN=str(workdir.parent / "cpmd.x"),
        CPMD_LAUNCH="",
        PP_LIBRARY_PATH=str(workdir.parent),
        CPMD_PP="14:Si_MT_BLYP.psp:LMAX=D,7:N_MT_BLYP.psp:LMAX=P",
        CPMD_CUTOFF="70",
    )
    env.pop("CPMD_RESTART_POOL", None)
    env.update(extra)
    return subprocess.run(
        [sys.executable, str(WRAPPER)], cwd=workdir, env=env,
        capture_output=True, text=True,
    )


@pytest.fixture
def workdir(tmp_path):
    fake = tmp_path / "cpmd.x"
    fake.write_text(FAKE_CPMD)
    fake.chmod(0o755)
    work = tmp_path / "extpot_1_0"
    work.mkdir()
    write_exchange(work, CELL, ATOMS)
    return work


def expected():
    centre = [sum(CELL[r][c] for r in range(3)) / 2 for c in range(3)]
    energy = sum(
        K_HA * sum(((x - c) / A0_CPMD) ** 2 for x, c in zip((x, y, w), centre))
        for _, x, y, w in ATOMS
    )
    forces = [
        [-2 * K_HA * (v - c) / A0_CPMD * HA_EV / A0_2018 for v, c in zip((x, y, w), centre)]
        for _, x, y, w in ATOMS
    ]
    return energy * HA_EV, forces


def read_result(workdir):
    rows = [l.split() for l in (workdir / "from_extpot_to_eon").read_text().splitlines()]
    return float(rows[0][0]), [[float(v) for v in r] for r in rows[1:]]


def test_energy_forces_units_sign_and_order(workdir):
    out = run_wrapper(workdir, CPMD_CONVERGENCE="1e-7", CPMD_MAXITER="123")
    assert out.returncode == 0, out.stderr
    energy, forces = read_result(workdir)
    e_ref, f_ref = expected()
    assert energy == pytest.approx(e_ref, rel=0, abs=E_ROUND_EV)
    assert len(forces) == len(ATOMS)
    for got, want in zip(forces, f_ref):
        assert got == pytest.approx(want, rel=1e-9, abs=1e-9)
    deck = (workdir / "cpmd.inp").read_text()
    assert " CONVERGENCE ORBITALS\n  1.000E-07\n" in deck
    assert " MAXITER\n  123\n" in deck
    assert "SYMMETRY" not in deck
    # The cell goes to CPMD as given, rows are lattice vectors.
    assert " CELL VECTORS\n  10.0000000000 0.0000000000 0.0000000000\n" in deck
    assert "  0.5000000000 0.0000000000 12.0000000000\n" in deck
    # CPMD reads species by species, in the order they first appear.
    assert deck.index("*N_MT_BLYP.psp") < deck.index("*Si_MT_BLYP.psp")
    assert (workdir / "argv").read_text().split() == ["cpmd.inp", str(workdir.parent)]


def test_second_call_restarts_from_the_first(workdir):
    assert run_wrapper(workdir).returncode == 0
    assert "RESTART" not in (workdir / "cpmd.inp").read_text()
    (workdir / "RESTART.1").write_text("stand-in")
    (workdir / "LATEST").write_text("./RESTART.1\n 1\n")
    assert run_wrapper(workdir).returncode == 0
    assert " RESTART WAVEFUNCTION LATEST" in (workdir / "cpmd.inp").read_text()


def test_no_convergence_retries_cold_with_pcg(workdir):
    (workdir / "no_convergence_once").write_text("")
    (workdir / "RESTART.1").write_text("stand-in")
    (workdir / "LATEST").write_text("./RESTART.1\n 1\n")
    out = run_wrapper(workdir)
    assert out.returncode == 0, out.stderr
    deck = (workdir / "cpmd.inp").read_text()
    assert " PCG MINIMIZE" in deck and "RESTART" not in deck
    assert (workdir / "retried_cold").exists()
    energy, _ = read_result(workdir)
    assert energy == pytest.approx(expected()[0], rel=0, abs=E_ROUND_EV)


def test_zero_cell_is_refused(workdir):
    write_exchange(workdir, [[0.0] * 3] * 3, ATOMS)
    out = run_wrapper(workdir)
    assert out.returncode != 0
    assert "zero cell" in out.stderr
    assert not (workdir / "from_extpot_to_eon").exists()
