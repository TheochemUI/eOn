import os
import numpy
import atoms
import io

class NotImplementedError(Exception):
    pass

class DisplaceError(Exception):
    pass

class Displace:
    def __init__(self, reactant, std_dev, radius):
        '''Reactant is an Atoms object. std_dev is the standard deviation
        of the normal distribution used to create the random displacements.
        radius is the distance to neighbors that will also be displaced.'''
        self.reactant = reactant
        self.std_dev = std_dev
        self.radius = radius

        #temporary numpy array of same size as self.reactant.r
        self.temp_array = numpy.zeros(self.reactant.r.shape)

        self.neighbors_list = atoms.sweep_and_prune(self.reactant, self.radius)

    def make_displacement(self, path):
        '''Writes the displacement_passed.con and mode_passed.dat to
        path.'''
        raise NotImplementedError

    def get_displacement(self, atom_index):
        '''Returns a displacement to be added to self.reactant.r'''
        #get neighboring atoms to the atom_index atom
        displaced_atoms = self.neighbors_list[atom_index]
        #add the selected atom to the list
        displaced_atoms.append(atom_index)

        displacement = numpy.zeros(self.reactant.r.shape)
        for i in xrange(len(displaced_atoms)):
            #don't displace frozen atoms
            if not self.reactant.free[displaced_atoms[i]]:
                continue
            #Displace one of the free atoms by a gaussian distributed
            #random number with a standard deviation of self.std_dev.
            displacement[displaced_atoms[i]] = numpy.random.normal(scale = self.std_dev, size=3)
        
        return displacement


    def save_files(self, path, displacement):
        self.temp_array = self.reactant.r.copy()
        self.reactant.r += displacement
        io.savecon(os.path.join(path, "displacement_passed.con"), self.reactant)
        self.reactant.r = self.temp_array
        mode = displacement/numpy.linalg.norm(displacement)
        mode = displacement
        io.save_mode(os.path.join(path, "mode_passed.dat"), displacement, self.reactant)

#class KDB()
#class Recycle()

class Undercoordinated(Displace):
    def __init__(self, reactant, max_coordination, std_dev=0.05, radius=5.0):
        Displace.__init__(self, reactant, std_dev, radius)

        self.max_coordination = max_coordination
        self.undercoordinated_atoms = []

        # XXX: 3.3 should be replaced by a fancy method that uses bond distances.
        self.coordination_distance = 3.3

        #cns is an array of the coordination numbers for each of the atoms.
        cns = atoms.coordination_numbers(self.reactant, self.coordination_distance)

        #Only allow displacements for atoms <= the maximum coordination and
        #that are free.
        self.undercoordinated_atoms = [ i for i in range(len(cns)) 
                if cns[i] <= self.max_coordination and 
                    self.reactant.free[i] == 1]
        if len(self.undercoordinated_atoms) == 0:
            errmsg = "No free atoms have a coordination of %i or less." 
            errmsg = errmsg % self.max_coordination
            raise DisplaceError(errmsg)

    def make_displacement(self, path):
        """Select an undercoordinated atom and displace all atoms in a radius 
        about it."""
        # TODO: We should make sure that the amount of I/O to disk
        #       is what we think it should be: about 100 kB or so per
        #       make_displacement().
        epicenter = self.undercoordinated_atoms[numpy.random.randint(len(self.undercoordinated_atoms))] 
        displacement = self.get_displacement(epicenter)
        self.save_files(path, displacement)

class Leastcoordinated(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0):
        Displace.__init__(self, reactant, std_dev, radius)

        self.leastcoordinated_atoms = []

        # XXX: 3.3 should be replaced by a fancy method that uses bond distances.
        self.coordination_distance = 3.3

        self.leastcoordinated_atoms = atoms.least_coordinated(self.reactant, 
                self.coordination_distance)
        self.leastcoordinated_atoms = [ i for i in self.leastcoordinated_atoms
                                        if self.reactant.free[i] == 1]

        if len(self.leastcoordinated_atoms) == 0:
            errmsg = "The least coordinated atoms are all frozen."
            raise DisplaceError(errmsg)

    def make_displacement(self, path):
        """Select an undercoordinated atom and displace all atoms in a radius about it."""
        # TODO: We should make sure that the amount of I/O to disk
        #       is what we think it should be: about 100 kB or so per
        #       make_displacement().
        epicenter = self.leastcoordinated_atoms[numpy.random.randint(len(self.leastcoordinated_atoms))] 
        displacement = self.get_displacement(epicenter)
        self.save_files(path, displacement)

class Random(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0):
        Displace.__init__(self, reactant, std_dev, radius)

        #each item in this list is the index of a free atom
        self.free_atoms = [ i for i in range(len(self.reactant.free)) 
                if self.reactant.free[i] ]

        if len(self.free_atoms) == 0: 
            raise DisplaceError("There are no free atoms in the reactant.")

    def make_displacement(self, path):
        """Select a random atom and displace all atoms in a radius about it."""
        #chose a random atom
        epicenter = self.free_atoms[numpy.random.randint(len(self.free_atoms))] 
        displacement = self.get_displacement(epicenter)
        self.save_files(path, displacement)

if __name__ == '__main__':
    import sys
    import time

    if len(sys.argv) < 3:
        print "%s: reactant.con outpath" % sys.argv[0]
        sys.exit(1)

    reactant = io.loadcon(sys.argv[1])
    d = Random(reactant, 0.05, 5.0)
    #d = Undercoordinated(reactant, 11, 0.05, 5.0)
    t0 = time.time()
    ntimes = 1000
    for i in xrange(ntimes):
        d.make_displacement(sys.argv[2])
    dt = time.time()-t0
    print "%.2f displacements per second" % (float(ntimes)/dt)
