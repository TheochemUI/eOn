from ase.calculators.nwchem import NWChem
from ase.calculators.socketio import *
import os
import copy
import psutil
import sys
from pathlib import Path

run_path = Path.cwd() / "runs"  # Where things are run
scratch_path = Path.cwd() / "nwchem_aux"  # scratch from NWChem
nwchem_path = (
    "/home/rgoswami/Git/Github/TheochemUI/nwchem_and_more/nwchem/bin/LINUX64/nwchem "
)
memory = "2 gb"  # Probably needs no change
# Multiplicity of the molecule. 1 for singlets, 2 for doublets.
mult = 1
label = "bah"
# run_dir = run_path / label
scratch_dir = scratch_path / label
nwchem_kwargs = dict(
    label=label,
    # directory=run_dir,
    perm=scratch_dir,
    scratch=scratch_dir,
    command=f"mpirun -n {psutil.cpu_count(logical=False)} {nwchem_path} PREFIX.nwi > PREFIX.nwo",
    # command=f'{nwchem_path} PREFIX.nwi > PREFIX.nwo',
    memory=memory,
    scf=dict(
        nopen=mult - 1,
        thresh=1e-8,
        maxiter=200,
    ),
    basis="3-21G",
    task="gradient",
    driver=dict(
        socket=f"unix {label}",
    ),
)
nwchem = NWChem(**nwchem_kwargs)
# Example
from ase.build import molecule

h2 = molecule("H2")
# h2.calc = nwchem
# print(h2.get_potential_energy())
# nwcalc = SocketIOCalculator(nwchem, unixsocket=label, log=sys.stdout)
# # Perform calculations based on your needs
# h2.calc=nwcalc
# energy = h2.get_potential_energy()
# forces = h2.get_forces()
# print(f"Energy: {energy:.6f} eV")
# print(f"Forces:\n{forces}")
# nwcalc.close()
with SocketIOCalculator(launch_client=nwchem, unixsocket=label, log=sys.stdout) as nwcalc:
    h2.calc = nwcalc
    # Perform calculations based on your needs
    nwcalc.calculate(h2)
    energy = h2.get_potential_energy()
    forces = h2.get_forces()
    print(f"Energy: {energy:.6f} eV")
    print(f"Forces:\n{forces}")
