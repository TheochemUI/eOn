from ase.calculators.nwchem import NWChem
from ase.calculators.socketio import SocketIOCalculator
import psutil
import sys
from pathlib import Path

# Example
from ase.build import molecule
h2 = molecule("H2")
with SocketIOCalculator(unixsocket='bah', log=sys.stdout) as nwcalc:
    h2.calc = nwcalc
    print(h2.get_potential_energy())
    print(h2.get_forces())
