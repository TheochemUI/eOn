from ase.calculators.nwchem import NWChem
from ase.calculators.socketio import SocketClient
import psutil
import sys
from pathlib import Path

run_path = Path.cwd() / "runs"  # Where things are run
scratch_path = Path.cwd() / "nwchem_aux"  # scratch from NWChem
nwchem_path = "/home/rgoswami/micromamba/envs/eongit/bin/nwchem"
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
    command=f'mpirun -n {psutil.cpu_count(logical=False)} {nwchem_path} PREFIX.nwi > PREFIX.nwo',
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
h2.calc = nwchem
client = SocketClient(unixsocket=label)
# client.calculate(h2, use_stress=False)
for i, _ in enumerate(client.irun(h2, use_stress=False)):
    print(".")
    # print(f'step {i}: {h2.get_potential_energy()}')
    # print(h2.get_forces())
