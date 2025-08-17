#!/usr/bin/env python3

import ase.io
from ase.calculators.nwchem import NWChem


def run_nwchem_with_ase():
    atoms = ase.io.read("pos.con")

    # Set up NWChem calculator with desired parameters
    calc = NWChem(
        label="nwchem_run",
        basis="3-21G",
        scf={"thresh": 1e-8, "maxiter": 2000},
        task="gradient",
    )

    atoms.calc = calc

    # Run a calculation to get energy and forces
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()

    print(f"Calculated Energy (eV): {energy}")
    print(f"Calculated Forces (eV/Å):\n{forces}")

    # Optionally, save results to a trajectory or log file
    atoms.write("nwchem_run.xyz")


if __name__ == "__main__":
    run_nwchem_with_ase()
