##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------
import logging
logger = logging.getLogger('displace')

import os, re
from math import cos, sin
import numpy

import atoms
import fileio as io
import config

class DisplacementManager:
    def __init__(self, reactant, moved_atoms):
        self.reactant = reactant
        # TODO: Remove all the pointless config.* crap
        if config.displace_random_weight > 0:
            self.random = Random(self.reactant,
                                 config.disp_magnitude, config.disp_radius,
                                 hole_epicenters=moved_atoms)
        if config.displace_under_coordinated_weight > 0:
            self.under = Undercoordinated(self.reactant, 
                                          config.disp_max_coord,
                                          config.disp_magnitude, config.disp_radius,
                                          hole_epicenters=moved_atoms,
                                          cutoff=config.comp_neighbor_cutoff,
                                          use_covalent=config.comp_use_covalent,
                                          covalent_scale=config.comp_covalent_scale)
        if config.displace_least_coordinated_weight > 0:
            self.least = Leastcoordinated(self.reactant,
                                          config.disp_magnitude, config.disp_radius,
                                          hole_epicenters=moved_atoms,
                                          cutoff=config.comp_neighbor_cutoff,
                                          use_covalent=config.comp_use_covalent,
                                          covalent_scale=config.comp_covalent_scale)
        if config.displace_not_FCC_HCP_weight > 0:
            self.not_FCC_HCP = NotFCCorHCP(self.reactant, 
                                           config.disp_magnitude,
                                           config.disp_radius,
                                           hole_epicenters=moved_atoms,
                                           cutoff=config.comp_neighbor_cutoff,
                                           use_covalent=config.comp_use_covalent,
                                           covalent_scale=config.comp_covalent_scale)
        if config.displace_not_TCP_BCC_weight > 0:
            self.not_TCP_BCC = NotTCPorBCC(self.reactant, 
                                           config.disp_magnitude,
                                           config.disp_radius,
                                           hole_epicenters=moved_atoms,
                                           cutoff=config.comp_neighbor_cutoff,
                                           use_covalent=config.comp_use_covalent,
                                           covalent_scale=config.comp_covalent_scale)
        if config.displace_listed_weight > 0:
            self.listed = ListedAtoms(self.reactant, 
                                      config.disp_magnitude, config.disp_radius,
                                      hole_epicenters=moved_atoms,
                                      cutoff=config.comp_neighbor_cutoff,
                                      use_covalent=config.comp_use_covalent,
                                      covalent_scale=config.comp_covalent_scale,
                                      displace_all=config.displace_all_listed)
        if config.displace_water_weight > 0:
            self.water = Water(self.reactant,
                               config.stdev_translation, config.stdev_rotation,
                               config.molecule_list, config.disp_at_random)
        # ### TShacked start
        if config.displace_not_TCP_weight > 0:
            self.not_TCP = NotTCP(self.reactant, 
                                           config.disp_magnitude,
                                           config.disp_radius,
                                           hole_epicenters=moved_atoms,
                                           cutoff=config.comp_neighbor_cutoff,
                                           use_covalent=config.comp_use_covalent,
                                           covalent_scale=config.comp_covalent_scale)
        # ### TShacked end
        total = 0.0
        total += config.displace_random_weight
        total += config.displace_listed_weight
        total += config.displace_not_FCC_HCP_weight
        total += config.displace_not_TCP_BCC_weight
        total += config.displace_under_coordinated_weight
        total += config.displace_least_coordinated_weight
        total += config.displace_water_weight
        # ### TShacked start
        total += config.displace_not_TCP_weight
        # ### TShacked end
        # If no fractions are defined, do 100% random displacements.
        if total == 0.0:
            total = 1.0
            self.plist = [1.0/total]
            self.random = Random(self.reactant, 
                                 config.disp_magnitude, config.disp_radius,
                                 hole_epicenters=moved_atoms)
        else:
            self.plist = [config.displace_random_weight/total]
        self.plist.append(self.plist[-1] + config.displace_listed_weight/total)
        self.plist.append(self.plist[-1] + config.displace_not_FCC_HCP_weight/total)
        self.plist.append(self.plist[-1] + config.displace_not_TCP_BCC_weight/total)
        self.plist.append(self.plist[-1] + config.displace_under_coordinated_weight/total)
        self.plist.append(self.plist[-1] + config.displace_least_coordinated_weight/total)
        self.plist.append(self.plist[-1] + config.displace_water_weight/total)
        # ### TShacked start
        self.plist.append(self.plist[-1] + config.displace_not_TCP_weight/total)
        # ### TShacked end

    def make_displacement(self):
        # ### TShacked start
        disp_types = ["random", "listed", "not_FCC_HCP", "not_TCP_BCC", "under", "least", "water","not_TCP"]
        # ### TShacked end
        r = numpy.random.random_sample()
        i = 0
        while self.plist[i] < r:
            i += 1
        disp_type = disp_types[i]
        if disp_type == "random":
            logger.debug("Made random displacement")
            return self.random.make_displacement()
        elif disp_type == "listed":
            logger.debug("Made listed atom displacement")
            return self.listed.make_displacement()
        elif disp_type == "under":
            logger.debug("Made under-coordinated displacement")
            return self.under.make_displacement()
        elif disp_type == "least":
            logger.debug("Made least-coordinated displacement")
            return self.least.make_displacement()
        elif disp_type == "not_FCC_HCP":
            logger.debug("Made not-FCC-or-HCP displacement")
            return self.not_FCC_HCP.make_displacement()
        elif disp_type == "not_TCP_BCC":
            logger.debug("Made not-TCP-or-BCC displacement")
            return self.not_TCP_BCC.make_displacement()
        elif disp_type == "water":
            logger.debug("Made water displacement")
            return self.water.make_displacement()
        # ### TShacked start
        elif disp_type == "not_TCP":
            logger.debug("Made not-TCP displacement")
            return self.not_TCP.make_displacement()
        # ### TShacked end
        raise DisplaceError()


class NotImplementedError(Exception):
    pass

class DisplaceError(Exception):
    pass

class Displace:
    def __init__(self, reactant, std_dev, radius, hole_epicenters):
        '''Reactant is an Atoms object. std_dev is the standard deviation
           of the normal distribution used to create the random displacements.
           radius is the distance to neighbors that will also be displaced.
        '''
        self.reactant = reactant
        self.std_dev = std_dev
        self.radius = radius
        self.hole_epicenters = hole_epicenters

        # temporary numpy array of same size as self.reactant.r
        self.temp_array = numpy.zeros(self.reactant.r.shape)

        self.neighbors_list = None

    def make_displacement(self):
        '''Writes the displacement_passed.con and mode_passed.dat to path.'''
        raise NotImplementedError

    def get_displacement(self, atom_index):
        '''Returns a displacement to be added to self.reactant.r'''
        if self.neighbors_list is None:
            self.neighbors_list = atoms.neighbor_list(self.reactant, self.radius, config.comp_brute_neighbors)
        displacement_norm = 0.0
        displacement = numpy.zeros(self.reactant.r.shape)
        if hasattr(atom_index, '__getitem__'):
            logger.debug("Displacement epicenters: ", atom_index)
            neighbors = [ self.neighbors_list[i] for i in xrange(len(self.neighbors_list)) if i in atom_index ]
            #flatten
            neighbors = sum(neighbors,[])
            neighbors = numpy.array(list(set(neighbors)), dtype=int)

            displaced_atoms = numpy.append(atom_index, neighbors)
        else:
            logger.debug("Displacement epicenter: %d" % atom_index)
            displaced_atoms = [atom_index] + self.neighbors_list[atom_index]

        # ensures that the total displacement vector exceeds a given length
        while (displacement_norm <= config.disp_min_norm):
            displacement = numpy.zeros(self.reactant.r.shape)
            for i in range(len(displaced_atoms)):
                # don't displace frozen atoms
                if not self.reactant.free[displaced_atoms[i]]:
                    continue
                # Displace one of the free atoms by a gaussian distributed
                # random number with a standard deviation of self.std_dev.
                displacement[displaced_atoms[i]] = numpy.random.normal(scale = self.std_dev, size=3)
                if config.displace_1d:
                    displacement[displaced_atoms[i]] = displacement[displaced_atoms[i]] * [1,0,0]
            displacement_norm = numpy.sqrt(numpy.sum(displacement.flatten()**2))

        displacement_atoms = self.reactant.copy()
        displacement_atoms.r += displacement
        
        if config.random_mode:
            displacement *= 0.0
            for i in range(len(displaced_atoms)):
                if not self.reactant.free[displaced_atoms[i]]:
                    continue
                displacement[displaced_atoms[i]] = numpy.random.normal(scale = self.std_dev, size=3)
                if config.displace_1d:
                    displacement[displaced_atoms[i]] = displacement[displaced_atoms[i]] * [1,0,0]
        
        displacement /= numpy.linalg.norm(displacement)
        return displacement_atoms, displacement/numpy.linalg.norm(displacement)

    def filter_epicenters(self, epicenters):
        '''Returns the epicenters that lie in the hole defined by Displace.hole_epicenters.
           If Displace.hole_epicenters is None, all of the epicenters are accepted.'''
        if self.hole_epicenters is None:
            return epicenters
        new_epicenters = []
        for e in epicenters:
            if e in self.hole_epicenters:
                new_epicenters.append(e)
        return new_epicenters

class Undercoordinated(Displace):
    def __init__(self, reactant, max_coordination, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.max_coordination = max_coordination
        self.undercoordinated_atoms = []

        self.coordination_distance = cutoff
        self.initialized = False

    def init(self):
        self.initialized = True
        # cns is an array of the coordination numbers for each of the atoms.
        cns = atoms.coordination_numbers(self.reactant, self.coordination_distance)

        # Only allow displacements for atoms <= the maximum coordination and that are free.
        self.undercoordinated_atoms = [ i for i in range(len(cns)) 
                if cns[i] <= self.max_coordination and 
                    self.reactant.free[i] == 1]

        self.undercoordinated_atoms = self.filter_epicenters(self.undercoordinated_atoms)

        if len(self.undercoordinated_atoms) == 0:
            errmsg = "No free atoms have a coordination of %i or less"
            errmsg = errmsg % self.max_coordination
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an undercoordinated atom and displace all atoms in a radius about it."""
        # TODO: We should make sure that the amount of I/O to disk is what we think it should be:
        #       about 100 kB or so per make_displacement().
        if not self.initialized:
            self.init()
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
            errmsg = "The least coordinated atoms are all frozen"
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an undercoordinated atom and displace all atoms in a radius about it."""
        epicenter = self.leastcoordinated_atoms[numpy.random.randint(len(self.leastcoordinated_atoms))] 
        return self.get_displacement(epicenter)

class NotFCCorHCP(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.not_HCP_or_FCC_atoms = []

        self.coordination_distance = cutoff

        self.not_HCP_or_FCC_atoms = atoms.not_HCP_or_FCC(self.reactant, 
                self.coordination_distance)

        self.not_HCP_or_FCC_atoms = [ i for i in self.not_HCP_or_FCC_atoms
                                        if self.reactant.free[i] == 1]

        self.not_HCP_or_FCC_atoms = self.filter_epicenters(self.not_HCP_or_FCC_atoms)

        if len(self.not_HCP_or_FCC_atoms) == 0:
            errmsg = "The atoms without FCC or HCP coordination are all frozen"
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an atom without HCP or FCC coordination and displace all atoms in a radius about it."""
        epicenter = self.not_HCP_or_FCC_atoms[numpy.random.randint(len(self.not_HCP_or_FCC_atoms))] 
        return self.get_displacement(epicenter)

class NotTCPorBCC(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.not_TCP_or_BCC_atoms = []

        self.coordination_distance = cutoff

        self.not_TCP_or_BCC_atoms = atoms.not_TCP_or_BCC(self.reactant, 
                self.coordination_distance)

        self.not_TCP_or_BCC_atoms = [ i for i in self.not_TCP_or_BCC_atoms
                                        if self.reactant.free[i] == 1]

        self.not_TCP_or_BCC_atoms = self.filter_epicenters(self.not_TCP_or_BCC_atoms)

        if len(self.not_TCP_or_BCC_atoms) == 0:
            errmsg = "The atoms without BCC or TCP coordination are all frozen."
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an atom without TCP or BCC coordination and displace all atoms in a radius about it."""
        epicenter = self.not_TCP_or_BCC_atoms[numpy.random.randint(len(self.not_TCP_or_BCC_atoms))] 
        return self.get_displacement(epicenter)

# ### TShacked start
class NotTCP(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.not_TCP_atoms = []

        self.coordination_distance = cutoff

        self.not_TCP_atoms = atoms.not_TCP(self.reactant, 
                self.coordination_distance)

        self.not_TCP_atoms = [ i for i in self.not_TCP_atoms
                                        if self.reactant.free[i] == 1]

        self.not_TCP_atoms = self.filter_epicenters(self.not_TCP_atoms)

        if len(self.not_TCP_atoms) == 0:
            errmsg = "The atoms without TCP coordination are all frozen"
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an atom without HCP or FCC coordination and displace all atoms in a radius about it."""
        epicenter = self.not_TCP_atoms[numpy.random.randint(len(self.not_TCP_atoms))] 
        return self.get_displacement(epicenter)
# ### TShacked end

class ListedAtoms(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None, cutoff=3.3, use_covalent=False, covalent_scale=1.3, displace_all=False):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        self.displace_all = displace_all
        # each item in this list is the index of a free atom
        self.listed_atoms = [ i for i in config.disp_listed_atoms 
                if self.reactant.free[i] ]

        self.listed_atoms = self.filter_epicenters(self.listed_atoms)

        if len(self.listed_atoms) == 0:
            raise DisplaceError("Listed atoms are all frozen")

    def make_displacement(self):
        """Select a listed atom and displace all atoms in a radius about it."""
        # chose a random atom from the supplied list
        if self.displace_all:
            epicenter = self.listed_atoms
        else:
            epicenter = self.listed_atoms[numpy.random.randint(len(self.listed_atoms))]
        return self.get_displacement(epicenter)

class Random(Displace):
    def __init__(self, reactant, std_dev=0.05, radius=5.0, hole_epicenters=None):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters)

        # each item in this list is the index of a free atom
        self.free_atoms = [ i for i in range(len(self.reactant.free))
                if self.reactant.free[i] ]

        self.free_atoms = self.filter_epicenters(self.free_atoms)

        if len(self.free_atoms) == 0: 
            raise DisplaceError("There are no free atoms in the reactant")

    def make_displacement(self):
        """Select a random atom and displace all atoms in a radius about it."""
        # chose a random atom
        epicenter = self.free_atoms[numpy.random.randint(len(self.free_atoms))]
        return self.get_displacement(epicenter)

class Water(Displace):
    '''Displace molecules of water without streatching them.'''
    def __init__(self, reactant, stdev_translation, stdev_rotation, molecule_list=[], random=0):
        """reactant: structure to be displaced\n"""\
        """stdev_translation: translational standard deviation (Angstrom)\n"""\
        """stdev_rotation: rotational standard deviation (radian)"""\
        """molecules: list of indices of the molecules to displace or None to displace all the molecules"""\
        """random: if 0 displace all molecules in molecule_list, if 'random > 0' picked up at random in 'molecule_list' a number of molecules equal to the number soted in 'random' and displace only these"""
        assert(isinstance(molecule_list, list))
        self.reactant = reactant
        self.stdev_translation = stdev_translation
        self.stdev_rotation = stdev_rotation
        for i, name in enumerate(reactant.names):
            if not re.search('^H', name): break
        # For water assume that all the hydrogen are listed first, then all the oxygen
        self.n_water = i/2
        if len(molecule_list) == 0: molecule_list=range(self.n_water)
        self.molecule_list = molecule_list
        self.random = random

    def make_displacement(self):
        '''Returns Atom object containing displaced structure and an array containing the displacement.'''
        return self.get_displacement()

    def get_displacement(self):
        '''Returns Atom object containing displaced structure and an array containing the displacement.'''
        free = self.reactant.free
        displaced_atoms = self.reactant.copy()
        if self.random > 0:
            n=len(self.molecule_list)
            molecule_list = list()
            for i in range(self.random):
                i = int(numpy.random.uniform(0, n))
                molecule_list.append(self.molecule_list[i])
        else: molecule_list = self.molecule_list
        for i in molecule_list:
            #print 'displacing', i
            h1 = i*2
            h2 = i*2+1
            o = i+self.n_water*2
            #don't displace if any of the three atoms is fixed
            #if not (free[h1] and free[h2] and free[o]):
            #    continue
            #Displace one of the free atoms by a gaussian distributed
            #random number with a standard deviation of self.std_dev.
            disp = numpy.random.normal(scale = self.stdev_translation, size=3)
            displaced_atoms.r[h1] += disp
            displaced_atoms.r[h2] += disp
            displaced_atoms.r[o] += disp
            rh1 = displaced_atoms.r[h1]
            rh2 = displaced_atoms.r[h2]
            ro = displaced_atoms.r[o]
            disp = numpy.random.normal(scale = self.stdev_rotation, size = 3)
            rh1, rh2, ro = self.rotate_water(rh1, rh2, ro, disp[0], disp[1], disp[2])
            displaced_atoms.r[h1] = rh1
            displaced_atoms.r[h2] = rh2
            displaced_atoms.r[o] = ro
        displacement = displaced_atoms.r - self.reactant.r
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
    def rotate_water(hydrogen1, hydrogen2, oxygen, psi, theta, phi, hydrogen_mass = 1.0, oxygen_mass = 16.0):
        G = (hydrogen_mass*(hydrogen1 + hydrogen2) + oxygen_mass*oxygen)/(hydrogen_mass*2.0 + oxygen_mass)
        rot = numpy.array([
            [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)],
            [sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), cos(theta)*sin(psi)],
            [cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(theta)*cos(psi)]  ])
        rh1 = numpy.tensordot(rot, (hydrogen1-G), 1) + G
        rh2 = numpy.tensordot(rot, (hydrogen2-G), 1) + G
        ro = numpy.tensordot(rot, (oxygen-G), 1) + G
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
