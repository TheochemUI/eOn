
import os
import sys
import glob

from eon import fileio as io

def make_movie(movie_type, path_root, states, separate_files=False):
    movie_path = "movie.poscar"
    if movie_type == 'dynamics':
        atoms_list = dynamics(path_root, states)
        movie_path = "dynamics.poscar"
    elif movie_type == 'states':
        atoms_list = dynamics(path_root, states, True)
        movie_path = "states.poscar"
    elif movie_type == 'fastestpath':
        atoms_list = fastest_path(path_root, states)
        movie_path = "fastestpath.poscar"
    elif movie_type == 'fastestfullpath':
        atoms_list = fastest_path(path_root, states, full=True)
        movie_path = "fastestfullpath.poscar"
    elif movie_type.split(',')[0] == 'processes':
        try:
            statenr = int(movie_type.split(',')[1])
        except ValueError:
            print("state number must be an integer")
            sys.exit(1)
        except IndexError:
            print("must give a state number")
            sys.exit(1)
        print("making process movie for state %i" % statenr)
        if len(movie_type.split(',')) > 2:
            limit = int(movie_type.split(',')[2])
        else:
            limit = 0
        atoms_list = processes(states, statenr, limit)
        movie_path = "processes_%i.poscar" % statenr
    elif movie_type == 'graph':
        s = dot(path_root, states)
        if os.path.isfile("graph.dot"):
            print("File %s already exists" % "graph.dot")
            sys.exit(1)
        f = open("graph.dot",'w')
        f.write(s)
        f.close()
        print("If you have graphviz installed:")
        print("dot graph.dot -Tpng -o graph.png")
        sys.exit(0)
    else:
        print("Unknown MOVIE_TYPE")
        sys.exit(1)
    if separate_files:
        movie_dir = "movies"
        if not os.path.exists(movie_dir):
            os.mkdir(movie_dir)
        movie_path = os.path.join(movie_dir, movie_path)
        # Delete existing files with the same root.
        for old_poscar in glob.glob(movie_path + ".*"):
            os.unlink(old_poscar)
        # Write new movie.
        for i, atoms in enumerate(atoms_list):
            path_i = movie_path + (".%010d" % i)
            io.saveposcar(path_i, atoms, 'w')
        print("Saved %i frames to %s.*" % (len(atoms_list), movie_path))
    else:
        # Delete existing file.
        if os.path.isfile(movie_path):
            os.unlink(movie_path)
        # Write movie.
        for atoms in atoms_list:
            io.saveposcar(movie_path, atoms, 'a')
        print("Saved %i frames to %s" % (len(atoms_list), movie_path))


def get_trajectory(trajectory_path):
    f = open(trajectory_path)
    trajectory = [0]
    f.readline()
    f.readline()
    for line in f:
        fields = line.split()
        statenr = int(fields[3])
        trajectory.append(statenr)

    return trajectory

def processes(states, statenr, limit):
    try:
        state = states.get_state(statenr)
    except IOError:
        print("error: Cannot make movie for non-existant state")
        sys.exit(1)

    process_table = state.get_process_table()
    for k,v in list(process_table.items()):
        process_table[k]['id'] = k
        process_table[k]['reactant'] = state.get_process_reactant(k)
        process_table[k]['saddle'] = state.get_process_saddle(k)
        process_table[k]['product'] = state.get_process_product(k)
    processes = list(process_table.values())
    sorted_processes = sorted(processes, key=lambda a: a['rate'])
    sorted_processes.reverse()

    atoms_list = []
    print("%4s %16s %16s %16s" % ("ID", "Rate", "Barrier", "Prefactor"))
    print("-------------------------------------------------------")
    for p in sorted_processes:
        atoms_list.append(p['reactant'])
        atoms_list.append(p['saddle'])
        atoms_list.append(p['product'])
        print("%4i %16.5e %16.5f %16.5e" % (p['id'], p['rate'] ,p['barrier'], p['prefactor']))
        limit -= 1
        if limit == 0:
            break
    return atoms_list

def dot(path_root, states):
    G = make_graph(states)
    return G.dot()

def dynamics(path_root, states, unique=False):
    if not unique:
        trajectory_path = os.path.join(path_root, "dynamics.txt")
        trajectory = get_trajectory(trajectory_path)
    else:
        trajectory = list(range(states.get_num_states()))

    atoms_list = []

    if len(trajectory) == 0:
        print("error: There have been no dynamics steps")
        sys.exit(1)

    for n in trajectory:
        state = states.get_state(n)
        reactant = state.get_reactant()
        atoms_list.append(reactant)

    return atoms_list

def make_graph(states):
    # Build the graph.
    G = Graph()
    for statenr in range(states.get_num_states()):
        state = states.get_state(statenr)
        G.add_node(state)
        ptable = state.get_process_table()
        for i,p in list(ptable.items()):
            if p['product'] != -1:
                neighbor_state = states.get_state(p['product'])
                G.add_node(neighbor_state)
                G.add_edge(state, neighbor_state, weight=1.0/p['rate'])
    return G

def fastest_path(path_root, states, full=False):
    G = make_graph(states)

    state_list = [states.get_state(0)]
    if full:
        nodes = sorted(G.nodes(), lambda a,b: a.number-b.number)
        state_pairs = [ nodes[i:i+2] for i in range(0, len(nodes)-1) ]
        for s1, s2 in state_pairs:
            path = G.shortest_path(s1,s2)[1:]
            for s in path:
                state_list.append(s)
    else:
        end = states.get_state(states.get_num_states()-1)
        state_list = G.shortest_path(state_list[0], end)

    atoms_list = []
    for i in range(len(state_list)):
        atoms_list.append(state_list[i].get_reactant())
        if len(state_list)-1 == i:
            continue
        atoms_list.append(state_list[i].get_process_saddle(
                                        get_fastest_process_id(state_list[i],
                                                               state_list[i+1])))
    time = 0.0
    for i in range(len(state_list)-1):
        ratesum = 0.0
        for j in list(state_list[i].get_process_table().values()):
            ratesum += j['rate']

        print(time, state_list[i].number)
        time += 1/ratesum

    return atoms_list

def get_fastest_process_id(state1, state2):
    ptable = state1.get_process_table()
    fastest = None
    for i,p in list(ptable.items()):
        if p['product'] == state2.number:
            if not fastest or fastest[1] < p['rate']:
                fastest = (i, p['rate'])
    return fastest[0]

def get_fastest_process_rate(state1, state2):
    ptable = state1.get_process_table()
    fastest = None
    for i,p in list(ptable.items()):
        if p['product'] == state2.number:
            if not fastest or fastest[1] < p['rate']:
                fastest = (i, p['rate'])
    return fastest[1]

class Graph:
    def __init__(self, name=""):
        self.name = name
        self.graph = {}

    def __repr__(self):
        return str(self.graph)

    def dot(self):
        unique_edges = set()
        for node,edgedict in list(self.graph.items()):
            edgelist = list(edgedict.keys())
            for edge in edgelist:
                if node.number < edge.number:
                    unique_edges.add( (node.number, edge.number) )
                else:
                    unique_edges.add( (edge.number, node.number) )

        s = "graph akmc {\n"
        for edge in unique_edges:
            s += "%s -- %s;\n" % edge
        s += "}\n"
        return s

    def add_node(self, node):
        if node not in self.graph:
            self.graph[node] = {}

    def add_edge(self, node1, node2, weight=1.0):
        if node1 not in self.graph:
            self.add_node(node1)
        if node2 not in self.graph:
            self.add_node(node2)

        self.graph[node1][node2] = weight

    def neighbors(self,node):
        if node in self.graph:
            return list(self.graph[node].keys())
        else:
            return None

    def nodes(self):
        return list(self.graph.keys())

    def dijkstra(self, start, end=None):
        G = self.graph
        D = {}	# dictionary of final distances
        P = {}	# dictionary of predecessors
        Q = priorityDictionary()   # est.dist. of non-final vert.
        Q[start] = 0

        for v in Q:
            D[v] = Q[v]
            if v == end: break

            for w in G[v]:
                vwLength = D[v] + G[v][w]
                if w in D:
                    if vwLength < D[w]:
                        raise ValueError("Dijkstra: found better path to already-final vertex")
                elif w not in Q or vwLength < Q[w]:
                    Q[w] = vwLength
                    P[w] = v

        return (D,P)

    def shortest_path(self, start, end):
        D,P = self.dijkstra(start,end)
        path = []
        while 1:
            path.append(end)
            if end == start: break
            end = P[end]
        path.reverse()
        return path

class priorityDictionary(dict):
    def __init__(self):
        '''Initialize priorityDictionary by creating binary heap
           of pairs (value,key).  Note that changing or removing a dict entry will
           not remove the old pair from the heap until it is found by smallest() or
           until the heap is rebuilt.
        '''
        self.__heap = []
        dict.__init__(self)

    def smallest(self):
        '''Find smallest item after removing deleted items from heap.'''
        if len(self) == 0:
            raise IndexError("smallest of empty priorityDictionary")
        heap = self.__heap
        while heap[0][1] not in self or self[heap[0][1]] != heap[0][0]:
            lastItem = heap.pop()
            insertionPoint = 0
            while 1:
                smallChild = 2*insertionPoint+1
                if smallChild+1 < len(heap) and \
                        heap[smallChild] > heap[smallChild+1]:
                    smallChild += 1
                if smallChild >= len(heap) or lastItem <= heap[smallChild]:
                    heap[insertionPoint] = lastItem
                    break
                heap[insertionPoint] = heap[smallChild]
                insertionPoint = smallChild
        return heap[0][1]

    def __iter__(self):
        '''Create destructive sorted iterator of priorityDictionary.'''
        def iterfn():
            while len(self) > 0:
                x = self.smallest()
                yield x
                del self[x]
        return iterfn()

    def __setitem__(self,key,val):
        '''Change value stored in dictionary and add corresponding
           pair to heap.  Rebuilds the heap if the number of deleted items grows
           too large, to avoid memory leakage.
        '''
        dict.__setitem__(self,key,val)
        heap = self.__heap
        if len(heap) > 2 * len(self):
            self.__heap = [(v,k) for k,v in list(self.items())]
            self.__heap.sort()  # builtin sort likely faster than O(n) heapify
        else:
            newPair = (val,key)
            insertionPoint = len(heap)
            heap.append(None)
            while insertionPoint > 0 and \
                    newPair < heap[(insertionPoint-1)//2]:
                heap[insertionPoint] = heap[(insertionPoint-1)//2]
                insertionPoint = (insertionPoint-1)//2
            heap[insertionPoint] = newPair

    def setdefault(self,key,val):
        '''Reimplement setdefault to call our customized __setitem__.'''
        if key not in self:
            self[key] = val
        return self[key]
