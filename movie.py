#!/usr/bin/env python
import os
import sys

import io

class Graph:
    def __init__(self,name=""):
        self.name = name
        self.list_neighbor = {}
        self.list_node = {}

    def add_node(self,node):
        self.list_node[node] = True

    def add_edge(self,node,nodebis):
        try :
            self.list_neighbor[node].append(nodebis)
        except :
            self.list_neighbor[node] = []
            self.list_neighbor[node].append(nodebis)
        try :
            self.list_neighbor[nodebis].append(node)
        except :
            self.list_neighbor[nodebis] = []
            self.list_neighbor[nodebis].append(node)

    def neighbors(self,node):
        try :
            return self.list_neighbor[node]
        except :
            return []

    def nodes(self):
        return self.list_node.keys()

    def shortest_path(self, node1, node2, path=[]):
        path = path + [node1]

        if node1 == node2:
            return path

        if node1 not in self.list_node or node2 not in self.list_node:
            return None

        shortest = None

        for node in self.neighbors(node1):
            if node not in path:
                newpath = self.shortest_path(node, node2, path)
                if newpath:
                    if not shortest or len(newpath) < len(shortest):
                        shortest = newpath
        return shortest

    def delete_edge(self,node,nodebis):
        self.list_neighbor[node].remove(nodebis)
        self.list_neighbor[nodebis].remove(node)

    def delete_node(self,node):
        del self.list_node[node]
        try :
            for nodebis in self.list_neighbor[node] :
                self.list_neighbor[nodebis].remove(node)
            del self.list_neighbor[node]
        except :
            return "error"

def get_trajectory(trajectory_path, unique=False):
    f = open(trajectory_path)
    trajectory = []
    for line in f:
        fields = line.split()
        statenr = int(fields[0])
        trajectory.append(statenr)

    if unique:
        trajectory = sorted(trajectory)
        trajectory = list(set(trajectory))

    return trajectory

def dynamics(path_root, states, unique=False):
    trajectory_path = os.path.join(path_root, "dynamics.txt")
    trajectory = get_trajectory(trajectory_path, unique)

    atoms_list = []

    if len(trajectory) == 0:
        print "error: there have been no dynamics steps"
        sys.exit(1)

    for n in trajectory:
        state = states.get_state(n)
        reactant = state.get_reactant()
        atoms_list.append(reactant)

    return atoms_list

def fastestpath(path_root, states):
    trajectory_path = os.path.join(path_root, "dynamics.txt")
    trajectory = get_trajectory(trajectory_path, True)

    G = Graph()

    for statenr in trajectory:
        state = states.get_state(statenr)
        G.add_node(state)
        ptable = state.get_process_table()
        for i,p in ptable.iteritems():
            if p['product'] != -1:
                neighbor_state = states.get_state(p['product'])
                G.add_node(neighbor_state)
                G.add_edge(state, neighbor_state) 

    nodes = map(states.get_state, trajectory)

    atoms_list = [states.get_state(0).get_reactant()]

    state_pairs = [ nodes[i:i+2] for i in range(0, len(nodes)-1) ]
    for s1, s2 in state_pairs:
        path = G.shortest_path(s1,s2)[1:]
        for s in path:
            atoms_list.append(s.get_reactant())
    return atoms_list

def save_movie(atoms_list, movie_path):
    for atoms in atoms_list:
        io.saveposcar(movie_path, atoms, 'a')

def make_movie(movie_type, path_root, states):
    if movie_type == 'dynamics':
        atoms_list = dynamics(path_root, states)
    elif movie_type == 'states':
        atoms_list = dynamics(path_root, states, True)
    elif movie_type == 'fastestpath':
        atoms_list = fastestpath(path_root, states)
    else:
        print "unknown MOVIE_TYPE"
        sys.exit(1)
    movie_path = "movie.poscar"
    if os.path.isfile(movie_path):
        print "file %s already exists" % movie_path
        sys.exit(1)
    save_movie(atoms_list, movie_path)
    print "saved %i frames to %s" % (len(atoms_list), movie_path)

#if __name__ == "__main__":
#    optpar = OptionParser(usage = "usage: %prog [options]")
#    optpar.add_option("-t", "--type", action="store", dest="type", default = "dynamics", help="specify the type of movie to make valid options are dynamics, states, and fastestpath")
#    (options, args) = optpar.parse_args()
#
#    # XXX: evil
#    import config
#
#    if options.type == "dynamics":
#        dynamics(config.path_results, config.path_states, "movie.poscar")
#    elif options.type == "fastestpath":
#        fastestpath(config.path_states, "movie.poscar")
#    elif options.type == "states":
#        dynamics(config.path_results, config.path_states, "movie.poscar", True)
