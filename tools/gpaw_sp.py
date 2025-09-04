#!/usr/bin/env gpaw-python
#XXX: gpaw-python has a MPI_Barrier call. This means
#     that all other MPI processes have to call the barrier
#     as well.

import ase
import gpaw
import time
import numpy
import sys
import os
from gpaw.mpi import world
from ase.utils import devnull

def create_gpaw(comm):
    from gpaw import GPAW, FermiDirac, PoissonSolver, setup_paths
    from gpaw import Mixer
    setup_paths.insert(0,'.')
    calc = GPAW(xc='PBE',
                occupations=FermiDirac(width=0.02),
                mixer = Mixer(beta=0.10, nmaxold=5, weight=100.0),
                communicator=comm)
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
client_ranks = []
potential_ranks = []
for i,t in enumerate(process_types):
    if t == 0:
        servers += 1
    elif t == 1:
        clients += 1
        client_ranks.append(i)
    elif t == 2:
        potentials += 1
        potential_ranks.append(i)
if "EON_CLIENT_STANDALONE" in os.environ:
    os.environ["EON_NUMBER_OF_CLIENTS"] = "1"
if "EON_NUMBER_OF_CLIENTS" not in os.environ:
    sys.stderr.write("you must set the env var EON_NUMBER_OF_CLIENTS\n");
    sys.exit(1)
clients = int(os.environ["EON_NUMBER_OF_CLIENTS"])

potential_group_size = potentials/clients

my_client_rank = client_ranks[potential_ranks.index(rank)/potential_group_size]
#print "pot: rank: %i my_client_rank: %i" % (world.rank, my_client_rank)

for i in range(clients):
    s = potential_group_size
    new_comm = world.new_communicator(numpy.array(potential_ranks[i*s:i*s+s], dtype='i'))
    if new_comm != None:
        my_comm = new_comm

first_time = True
nforce_calls = 0
performance_log = ""
while True:
    nforce_calls += 1
    natoms = numpy.array((0,), 'i')
    if my_comm.rank == 0:
        world.receive(natoms, my_client_rank, tag=0)
    my_comm.broadcast(natoms, 0)

    atomic_numbers = numpy.zeros(natoms, 'i')
    positions = numpy.zeros(3*natoms, 'd')
    cell = numpy.zeros(9, 'd')
    pbc = numpy.array((0,), 'i')
    logdir = numpy.zeros(1024, 'l')
    if my_comm.rank == 0:
        world.receive(atomic_numbers, my_client_rank, tag=0)
        world.receive(positions,      my_client_rank, tag=0)
        world.receive(cell,           my_client_rank, tag=0)
        world.receive(pbc,            my_client_rank, tag=0)
        world.receive(logdir,        my_client_rank, tag=0)
    my_comm.broadcast(atomic_numbers, 0)
    my_comm.broadcast(positions,      0)
    my_comm.broadcast(cell,           0)
    my_comm.broadcast(pbc,            0)
    my_comm.broadcast(logdir,         0)

    t0 = time.time()
    tmp = []
    for x in logdir:
        if x == 0: break
        tmp.append(chr(x))
    logdir = ''.join(tmp)

    #XXX: is this really needed? Why did it happen...?
    while not os.path.isdir(logdir):
        sys.stderr("error: logdir: %s didn't exist somehow\n" % logdir)
        time.sleep(0.1)

    if pbc == 1:
        pbc = True
    else:
        pbc = False

    cell.shape = (3,3)
    positions.shape = (natoms,3)
    if first_time:
        atomic_symbols = ''.join([ ase.chemical_symbols[int(i)] for i in atomic_numbers])
        atoms = ase.Atoms(atomic_symbols, positions=positions, cell=cell, pbc=pbc)

        calc = create_gpaw(my_comm)
        atoms.set_calculator(calc)
        first_time = False
    else:
        atoms.set_positions(positions)
    logfile = os.path.join(logdir, "gpaw_%i.txt"%nforce_calls)
    calc.set(txt=logfile)

    calculation_failed  = numpy.array((0,),'i')
    try:
        f1 = atoms.get_forces()
        e1 = atoms.get_potential_energy()
        e1 = numpy.array([e1,])
    except gpaw.KohnShamConvergenceError:
        calculation_failed  = numpy.array(1,'i')
    t1 = time.time()

    if my_comm.rank == 0:
        performance_log = os.path.join(logdir, "performance.txt")
        fperformance = open(performance_log, "a+")
        fperformance.write("%i %.3f\n" % (nforce_calls, (t1-t0)))
        fperformance.close()
        world.send(calculation_failed, my_client_rank, tag=0)
        if not calculation_failed:
            world.send(e1, my_client_rank, tag=0)
            world.send(f1, my_client_rank, tag=0)
