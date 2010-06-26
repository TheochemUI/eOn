#!/usr/bin/env python
import os
import sys
from optparse import OptionParser

import pathfix
import io
import statelist

def get_trajectory(trajectory_path):
    f = open(trajectory_path)
    states = []
    for line in f:
        fields = line.split()
        statenr = int(fields[0])
        states.append(statenr)
    return states

def dynamics(path_results, path_states, movie_path, unique=False):
    trajectory_path = os.path.join(path_results, "dynamics.txt")

    trajectory = get_trajectory(trajectory_path)

    num_steps = len(trajectory)

    if unique:
        trajectory = sorted(trajectory)
        trajectory = set(trajectory)

    if num_steps == 0:
        print "error: there have been no dynamics steps"
        sys.exit(1)

    atoms = io.loadcon(os.path.join(path_results, "reactant.con"))
    io.saveposcar(movie_path, atoms, 'a')

    for state in trajectory:
        atoms = io.loadcon(os.path.join(path_states, str(state), "reactant.con"))
        io.saveposcar(movie_path, atoms, 'a')

if __name__ == "__main__":
    optpar = OptionParser(usage = "usage: %prog [options]")
    optpar.add_option("-t", "--type", action="store", dest="type", default = "dynamics", help="specify the type of movie to make")
    (options, args) = optpar.parse_args()

    # XXX: evil
    import config

    if options.type == "dynamics":
        dynamics(config.path_results, config.path_states, "movie.poscar")
    elif options.type == "states":
        dynamics(config.path_results, config.path_states, "movie.poscar", True)
