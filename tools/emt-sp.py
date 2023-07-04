#!/usr/bin/env gpaw-python
#XXX: gpaw-python has a MPI_Barrier call. This means
#     that all other MPI processes have to call the barrier
#     as well.

import ase
import gpaw
import time
import numpy
import sys
from gpaw.mpi import world
from ase.utils import devnull

def create_gpaw(comm):
    from gpaw import GPAW, FermiDirac, PoissonSolver, setup_paths
    from gpaw import Mixer
    setup_paths.insert(0,'.')
    # gpaw calculator:
    convergence = {
                    'energy':0.001,
                    'density':1e-2,
                    'eigenstates':1e-5
                  }
    #calc = GPAW(xc='PBE',
    #            h=.30,
    #            nbands=-8,
    #            txt='gpaw_%i.txt'%world.rank,
    #            convergence=convergence,
    #            occupations=FermiDirac(width=0.05),
    #            mixer = Mixer(beta=0.10, nmaxold=5, weight=100.0),
    #            communicator=comm)
    from ase.calculators.emt import EMT
    calc = EMT()
    return calc

rank = world.rank

# process_type can be one of three values:
# 0: server
# 1: client
# 2: potential
# Each independent program identifies itself with one of these three.

process_type = numpy.array((2,), dtype='i')
process_types = numpy.empty(world.size, dtype='i')

world.all_gather(process_type, process_types)

servers = 0
clients = 0
potentials = 0
for i, t in enumerate(process_types):
    if t == 0:
        server_rank = i
        servers += 1
    elif t == 1:
        clients += 1
    elif t == 2:
        potentials += 1

potential_group_size = potentials/clients
potential_ranks = numpy.empty(potentials, dtype='i')
client_ranks = numpy.empty(clients, dtype='i')
j = 0
k = 0
first_potential_rank = None
for i in range(world.size):
    if process_types[i] == 1:
        client_ranks[k] = i
        k += 1
    elif process_types[i] == 2:
        if first_potential_rank == None:
            first_potential_rank = i
        potential_ranks[j] = i
        j += 1

my_potential_rank = world.rank-first_potential_rank

my_client_rank = client_ranks[my_potential_rank/potential_group_size]
print("pot: rank: %i my_client_rank: %i" % (world.rank, my_client_rank))

for i in range(clients):
    s = potential_group_size
    new_comm = world.new_communicator(potential_ranks[i*s:i*s+s])
    if new_comm != None:
        my_comm = new_comm

first_time = True
while True:
    natoms = numpy.array((0,), 'i')
    if my_comm.rank == 0:
        world.receive(natoms, my_client_rank, tag=0)
    my_comm.broadcast(natoms, 0)

    atomic_numbers = numpy.zeros(natoms, 'i')
    positions = numpy.zeros(3*natoms, 'd')
    cell = numpy.zeros(9, 'd')
    pbc = numpy.array((0,), 'i')
    if my_comm.rank == 0:
        world.receive(atomic_numbers, my_client_rank, tag=0)
        world.receive(positions,      my_client_rank, tag=0)
        world.receive(cell,           my_client_rank, tag=0)
        world.receive(pbc,            my_client_rank, tag=0)
    my_comm.broadcast(atomic_numbers, 0)
    my_comm.broadcast(positions,      0)
    my_comm.broadcast(cell,           0)
    my_comm.broadcast(pbc,            0)

    if pbc == 1:
        pbc = True
    else:
        pbc = False

    if first_time:
        cell.shape = (3,3)
        positions.shape = (natoms,3)

        atomic_symbols = ''.join([ ase.chemical_symbols[int(i)] for i in atomic_numbers])
        atoms = ase.Atoms(atomic_symbols, positions=positions, cell=cell, pbc=pbc)

        calc = create_gpaw(my_comm)
        atoms.set_calculator(calc)
    else:
        atoms.set_positions(positions)

    calculation_failed  = numpy.array((0,),'i')
    try:
        f1 = atoms.get_forces()
        e1 = atoms.get_potential_energy()
        e1 = numpy.array([e1,])
    except gpaw.KohnShamConvergenceError:
        calculation_failed  = numpy.array(1,'i')

    if my_comm.rank == 0:
        world.send(calculation_failed, my_client_rank, tag=0)
        if not calculation_failed:
            world.send(e1, my_client_rank, tag=0)
            world.send(f1, my_client_rank, tag=0)
