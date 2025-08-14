# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "ase",
# ]
# ///

import numpy as np
from ase import Atoms
from ase.io import write


def main():
    """
    Creates a simple Si-Au dimer and writes it to EON .con and
    LAMMPS data formats.
    """
    print("Creating Si-Au dimer structure...")

    # Define the atoms and their positions
    symbols = ["Si", "Au"]
    positions = [[0.0, 0.0, 0.0], [1.2, 1.3, 1.4]]
    distance = np.linalg.norm(np.array(positions[0]) - np.array(positions[1]))
    print(f"Distance between atoms: {distance:.4f} Å")

    # Create the ASE Atoms object
    atoms = Atoms(symbols=symbols, positions=positions)

    atoms.set_cell([20, 20, 20])
    atoms.pbc = True
    atoms.center()

    # Write the EON .con file
    print("Writing to system.con...")
    write("pos.con", atoms, format="eon")

    # Write the LAMMPS data file
    # 'atomic' atom_style, is the simplest.
    print("Writing to data.lammps...")
    write("data.lammps", atoms, format="lammps-data", atom_style="atomic")

    print("Structure generation complete.")


if __name__ == "__main__":
    main()
