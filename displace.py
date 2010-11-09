import os, re
import numpy
import atoms
import io
from math import cos, sin

class NotImplementedError(Exception):
    pass

class DisplaceError(Exception):
    pass

class Displace:
    def __init__(self, reactant, std_dev, radius, hole_epicenters):
        '''Reactant is an Atoms object. std_dev is the standard deviation
        of the normal distribution used to create the random displacements.
        radius is the distance to neighbors that will also be displaced.'''
        self.reactant = reactant
        self.std_dev = std_dev
        self.radius = radius
        self.hole_epicenters = hole_epicenters

        #temporary numpy array of same size as self.reactant.r
        self.temp_array = numpy.zeros(self.reactant.r.shape)

        self.neighbors_list = atoms.sweep_and_prune(self.reactant, self.radius)

    def make_displacement(self):
        '''Writes the displacement_passed.con and mode_passed.dat to
        path.'''
        raise NotImplementedError

    def get_displacement(self, atom_index):
        '''Returns a displacement to be added to self.reactant.r'''
        #get neighboring atoms to the atom_index atom
        #and add the selected atom to the list
        displaced_atoms = self.neighbors_list[atom_index] + [atom_index]

        displacement = numpy.zeros(self.reactant.r.shape)
        for i in range(len(displaced_atoms)):
            #don't displace frozen atoms
            if not self.reactant.free[displaced_atoms[i]]:
                continue
            #Displace one of the free atoms by a gaussian distributed
            #random number with a standard deviation of self.std_dev.
            displacement[displaced_atoms[i]] = numpy.random.normal(scale = self.std_dev, size=3)
        
        displacement_atoms = self.reactant.copy()
        displacement_atoms.r += displacement
        return displacement_atoms, displacement

    def filter_epicenters(self, epicenters):
        '''Returns the epicenters that lie in the hole defined by Displace.hole_epicenters.
           If Displace.hole_epicenters is None, all of the epicenters are accepted.'''
        cutoff = 9.9 #The size of the hole. TODO: this should be parameterized.
        if self.hole_epicenters == None:
            return epicenters
        new_epicenters = []
        for e in epicenters:
            for h in self.hole_epicenters:
                if numpy.linalg.norm(atoms.pbc(self.reactant.r[h] - self.reactant.r[e], self.reactant.box)) < cutoff:
                    new_epicenters.append(e)
                    break
        return new_epicenters

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
    def __init__(self, reactant, max_coordination, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.max_coordination = max_coordination
        self.undercoordinated_atoms = []

        self.coordination_distance = cutoff

        #cns is an array of the coordination numbers for each of the atoms.
        cns = atoms.coordination_numbers(self.reactant, self.coordination_distance)

        #Only allow displacements for atoms <= the maximum coordination and
        #that are free.
        self.undercoordinated_atoms = [ i for i in range(len(cns)) 
                if cns[i] <= self.max_coordination and 
                    self.reactant.free[i] == 1]
                    
        self.undercoordinated_atoms = self.filter_epicenters(self.undercoordinated_atoms)
                    
        if len(self.undercoordinated_atoms) == 0:
            errmsg = "No free atoms have a coordination of %i or less." 
            errmsg = errmsg % self.max_coordination
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an undercoordinated atom and displace all atoms in a radius 
        about it."""
        # TODO: We should make sure that the amount of I/O to disk
        #       is what we think it should be: about 100 kB or so per
        #       make_displacement().
        epicenter = self.undercoordinated_atoms[numpy.random.randint(len(self.undercoordinated_atoms))] 
        return self.get_displacement(epicenter)

class Leastcoordinated(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.leastcoordinated_atoms = []

        self.coordination_distance = cutoff

        self.leastcoordinated_atoms = atoms.least_coordinated(self.reactant, 
                self.coordination_distance)
        self.leastcoordinated_atoms = [ i for i in self.leastcoordinated_atoms
                                        if self.reactant.free[i] == 1]

        self.leastcoordinated_atoms = self.filter_epicenters(self.leastcoordinated_atoms)

        if len(self.leastcoordinated_atoms) == 0:
            errmsg = "The least coordinated atoms are all frozen."
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an undercoordinated atom and displace all atoms in a radius about it."""
        # TODO: We should make sure that the amount of I/O to disk
        #       is what we think it should be: about 100 kB or so per
        #       make_displacement().
        epicenter = self.leastcoordinated_atoms[numpy.random.randint(len(self.leastcoordinated_atoms))] 
        return self.get_displacement(epicenter)

class Random(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        #each item in this list is the index of a free atom
        self.free_atoms = [ i for i in range(len(self.reactant.free)) 
                if self.reactant.free[i] ]
                
        self.free_atoms = self.filter_epicenters(self.free_atoms)

        if len(self.free_atoms) == 0: 
            raise DisplaceError("There are no free atoms in the reactant.")

    def make_displacement(self):
        """Select a random atom and displace all atoms in a radius about it."""
        #chose a random atom
        epicenter = self.free_atoms[numpy.random.randint(len(self.free_atoms))] 
        return self.get_displacement(epicenter)

class Water(Displace):
    '''Displace molecules of water without streatching them'''
    def __init__(self, reactant, stdev_translation, stdev_rotation):
        '''reactant: structure to be displaced\n'''\
        '''stdev_translation: translational standard deviation (Angstrom)\n'''\
        '''stdev_rotation: rotational standard deviation (radian)'''
        self.reactant=reactant
        self.stdev_translation=stdev_translation
        self.stdev_rotation=stdev_rotation
        for i, name in enumerate(reactant.names):
            if not re.search('^H', name): break
        # For water assume that all the hydrogen are listed first, then all the oxygen
        self.n_water=i/2

    def make_displacement(self):
        '''Returns Atom object containing displaced structure and an array containing the displacement'''
        return self.get_displacement()

    def get_displacement(self):
        '''Returns Atom object containing displaced structure and an array containing the displacement'''
        free=self.reactant.free
        displaced_atoms = self.reactant.copy()
        i=numpy.random.uniform(0, self.n_water)
        h1=i*2
        h2=i*2+1
        o=i+self.n_water*2
        #don't displace if any of the three atoms is fixed
        #if not (free[h1] and free[h2] and free[o]):
        #    continue
        #Displace one of the free atoms by a gaussian distributed
        #random number with a standard deviation of self.std_dev.
        disp=numpy.random.normal(scale = self.stdev_translation, size=3)
        displaced_atoms.r[h1]+=disp
        displaced_atoms.r[h2]+=disp
        displaced_atoms.r[o]+=disp
        rh1=displaced_atoms.r[h1]
        rh2=displaced_atoms.r[h2]
        ro=displaced_atoms.r[o]
        disp=numpy.random.normal(scale = self.stdev_rotation, size=3)
        rh1, rh2, ro=self.rotate_water(rh1, rh2, ro, disp[0], disp[1], disp[2])
        displaced_atoms.r[h1]=rh1
        displaced_atoms.r[h2]=rh2
        displaced_atoms.r[o]=ro
        displacement=displaced_atoms.r - self.reactant.r
        return displaced_atoms, displacement

    ## Rotate a molecule of water.
    # Rotate around the centre of gravity of the molecule.
    # @param[in] "hydrogen1, hydrogen2, oxygen" numpy.array: coordinates of atoms.
    # @param[in] "psi, theta, phi" float: Angle in @em radians of rotation around <em> x, y, z </em> axes.
    # @param[in] "hydrogenMass, oxygenMass" float: masses of the hydrogen and oxygen atoms.
    # @return "hydrogen1, hydrogen2, oxygen" coordinates of atoms after rotation
    # The rotations uses the @em x, y, z convention (pitch-roll-yaw). The fisrt rotation to take place is around z, then y, then x.
    # Equations and notations are from: http://mathworld.wolfram.com/EulerAngles.html .
    @staticmethod
    def rotate_water(hydrogen1, hydrogen2, oxygen, psi, theta, phi, hydrogen_mass=1.0, oxygen_mass=16.0):
        G=(hydrogen_mass*(hydrogen1+hydrogen2)+oxygen_mass*oxygen)/(hydrogen_mass*2+oxygen_mass)
        rot=numpy.array([
            [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)],
            [sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), cos(theta)*sin(psi)],
            [cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(theta)*cos(psi)]  ])
        rh1=numpy.tensordot(rot, (hydrogen1-G), 1)+G
        rh2=numpy.tensordot(rot, (hydrogen2-G), 1)+G
        ro=numpy.tensordot(rot, (oxygen-G), 1)+G
        return rh1, rh2, ro

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
