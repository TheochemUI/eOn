#!/home/chill/local/bin/gpaw-python
import ase
from gpaw import GPAW
import gpaw
from mpi4py import MPI
import time
import numpy
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
parent_comm = comm.Get_parent()

while True:
    if rank == 0:
        die = numpy.array(0,'i')
        parent_comm.Recv([die, MPI.INT])
    else:
        die = None
    die = comm.bcast(die)
    if die:
        sys.exit(0)

    if rank == 0:
        natoms = numpy.array(0, 'l')
        parent_comm.Recv([natoms, MPI.LONG])
        #print "got natoms:", natoms
        atomic_numbers = numpy.zeros(natoms, 'i')
        positions = numpy.zeros(3*natoms, 'd')
        cell = numpy.zeros(9, 'd')
        pbc = numpy.array(0, 'i')

        parent_comm.Recv([atomic_numbers, MPI.INT])
        parent_comm.Recv([positions, MPI.DOUBLE])
        parent_comm.Recv([cell, MPI.DOUBLE])
        parent_comm.Recv([pbc, MPI.INT])
        #print "atomic numbers:", atomic_numbers
        #print "positions:", positions
        #print "cell:",cell 
        #print "pbc:", pbc
        if pbc == 1:
            pbc = True
        else:
            pbc = False

        cell.shape = (3,3)
        positions.shape = (natoms,3)

        atomic_symbols = ''.join([ ase.chemical_symbols[i] for i in atomic_numbers])
        atoms = ase.Atoms(atomic_symbols, positions=positions, cell=cell, pbc=pbc)
    else:
        atoms = None
    atoms = comm.bcast(atoms)

    # gpaw calculator:
    convergence = {'energy':0.001}
    calc = GPAW(xc='LDA', txt=None, convergence=convergence)
    atoms.set_calculator(calc)

    calculation_failed  = numpy.array(0,'i')
    try:
        f1 = atoms.get_forces()
        e1 = atoms.get_potential_energy()
    except gpaw.KohnShamConvergenceError:
        calculation_failed  = numpy.array(1,'i')

    if rank == 0:
        parent_comm.Send([calculation_failed, MPI.INT])
        if not calculation_failed:
            parent_comm.Send([e1, MPI.DOUBLE])
            parent_comm.Send([f1, MPI.DOUBLE])
