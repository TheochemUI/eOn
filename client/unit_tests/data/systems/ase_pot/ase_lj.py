"""ASE LJ calculator for eOn test suite."""
import numpy as np
from ase import Atoms
from ase.calculators.lj import LennardJones


def ase_calc():
    return LennardJones(epsilon=1.0, sigma=1.0, rc=10.0, ro=0.0, smooth=False)


def _calculate(R, atomicNrs, box, calc):
    system = Atoms(symbols=atomicNrs, positions=R, pbc=True, cell=box)
    system.calc = calc
    forces = system.get_forces()
    energy = system.get_potential_energy()
    return energy, forces


def batch_calculate(Rs, atomicNrs, boxes, calc):
    n = Rs.shape[0]
    energies = np.empty(n)
    forces = np.empty((n, Rs.shape[1], 3))
    for i in range(n):
        energies[i], forces[i] = _calculate(Rs[i], atomicNrs[i], boxes[i], calc)
    return energies, forces
