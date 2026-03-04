"""ASE LJ calculator for eOn test suite."""
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
