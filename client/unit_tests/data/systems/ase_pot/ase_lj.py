"""ASE LJ calculator for eOn test suite."""

from ase import Atoms
from ase.calculators.lj import LennardJones
import numpy as np


def ase_calc():
    return LennardJones(epsilon=1.0, sigma=1.0, rc=10.0, ro=0.0, smooth=False)


def _calculate(R, atomicNrs, box, calc):
    system = Atoms(symbols=atomicNrs, positions=R, pbc=True, cell=box)
    system.calc = calc
    forces = system.get_forces()
    energy = system.get_potential_energy()
    return energy, forces


def batch_calculate(Rs, atomicNrs, boxes, calc):
    num_structures = np.shape(Rs)[0]
    energies = np.empty(num_structures)
    forces = np.empty((num_structures, atomicNrs, 3))
    for structure_idx in range(num_structures):
        energies[structure_idx], forces[structure_idx, :] = _calculate(
            Rs[structure_idx, :, :],
            atomicNrs[structure_idx, :],
            boxes[structure_idx, :, :],
            calc,
        )
    return energies, forces
