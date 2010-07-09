""" The atoms module. """

import numpy
import logging
logger = logging.getLogger('atoms')
#TODO: Add list of covalent radii from tsse/jmol/kdb. (Rye can do)



class Atoms:
    """ The Atoms class. """

    def __init__(self, n_atoms):
        self.r = numpy.zeros((n_atoms,3))
        self.free = numpy.ones(n_atoms)
        self.box = numpy.zeros((3,3))
        self.names = ['']*n_atoms

    def __len__(self):
        '''
        Returns the number of atoms in the object'''
        return len(self.r)
        
    def copy(self):
        p = Atoms(len(self))
        p.r = self.r.copy()
        p.free = self.free.copy()
        p.box = self.box.copy()
        p.names = self.names[:]
        return p

    def __deepcopy__(self):
        return self.copy()

def pbc(r, box, ibox = None):
    """
    Applies periodic boundary conditions.
    Parameters:
        r:      the vector the boundary conditions are applied to
        box:    the box that defines the boundary conditions
        ibox:   the inverse of the box. This will be calcluated if not provided.
    """
    if ibox == None:    
        ibox = numpy.linalg.inv(box)
    vdir = numpy.dot(r, ibox)
    vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    return numpy.dot(vdir, box)
    
def per_atom_norm(v, box, ibox = None):
    '''
    Returns a length N numpy array containing per atom distance
        v:      an Nx3 numpy array
        box:    box matrix that defines the boundary conditions
        ibox:   the inverse of the box. will be calculated if not provided
    '''
    diff = pbc(v, box, ibox)
    return numpy.array([numpy.linalg.norm(d) for d in diff])

def per_atom_norm_gen(v, box, ibox = None):
    '''
    Returns a generator which yields the distance between pairs of atoms
        v:      an Nx3 numpy array
        box:    box matrix that defines the boundary conditions
        ibox:   the inverse of the box. will be calculated if not provided
    '''
    diff = pbc(v, box, ibox)
    for d in diff:
        yield numpy.linalg.norm(d)

def identical(atoms1, atoms2, epsilon_r):
    '''
    Determines whether two structures are identical if atoms of the same element are
    considered indistinguishable.
        atoms1:     first atoms object for comparison
        atoms2:     second atoms object for comparison
        epsilon_r:  distance (in angstroms) that two atoms must be seperated by in order to be considered different
    '''
    #XXX: n^2
    if len(atoms1) != len(atoms2):
        return False
    
    for i in range(3):
        for j in range(3):
            #XXX: Hardcoded comparison criteria
            if abs(atoms1.box[i][j] - atoms2.box[i][j]) > .0001:
                logger.warning("Identical returned false because boxes were not the same")
                return False
    box = atoms1.box 
    ibox = numpy.linalg.inv(box)

    mismatch = []
    pan = per_atom_norm(atoms1.r - atoms2.r, box, ibox)
    for i in range(len(pan)):
        if pan[i] > epsilon_r:
            mismatch.append(i)
    
    for i in mismatch:
        pan = per_atom_norm(atoms1.r - atoms2.r[i], box, ibox)
        minpan = 1e300
        minj = 0
        for j in range(len(pan)):
            if i == j:
                continue
            minpan = min(minpan, pan[j])
        if not minpan < epsilon_r:
            return False

    return True
            
            
def sweep_and_prune(p_in, cutoff, strict = True, bc = True):
    """ Returns a list of nearest neighbors within cutoff for each atom. 
        Parameters:
            p_in:   Atoms object
            cutoff: the radius within which two atoms are considered to intersect.
            strict: perform an actual distance check if True
            bc:     include neighbors across pbc's """
    #TODO: Get rid of 'cutoff' and use the covalent bond radii. (Rye can do)
        # Do we want to use covalent radii? I think the displace class wants to allow for user-defined cutoffs.
            # We should have both options available. -Rye
    #TODO: Make work for nonorthogonal boxes.
    p = p_in.copy()
    p.r = pbc(p.r, p.box)
    p.r -= numpy.array([min(p.r[:,0]), min(p.r[:,1]), min(p.r[:,2])]) 
    numatoms = len(p)
    coord_list = []
    for i in range(numatoms):
        coord_list.append([i, p.r[i]])
    for axis in range(3):
        sorted_axis = sorted(coord_list, key = lambda foo: foo[1][axis])
        intersect_axis = []
        for i in range(numatoms):
            intersect_axis.append([])
        for i in range(numatoms):
            done = False
            j = i + 1
            if not bc and j >= numatoms:
                done = True
            while not done:
                j = j % numatoms
                if j == i:
                    done = True
                dist = abs(sorted_axis[j][1][axis] - sorted_axis[i][1][axis])
                if p.box[axis][axis] - sorted_axis[i][1][axis] < cutoff:
                    dist = min(dist, (p.box[axis][axis] - sorted_axis[i][1][axis]) + sorted_axis[j][1][axis])
                if dist < cutoff:
                    intersect_axis[sorted_axis[i][0]].append(sorted_axis[j][0]) 
                    intersect_axis[sorted_axis[j][0]].append(sorted_axis[i][0]) 
                    j += 1
                    if not bc and j >= numatoms:
                        done = True
                else:
                    done = True
        if axis == 0:
            intersect = []
            for i in range(numatoms):
                intersect.append([])
                intersect[i] = intersect_axis[i]
        else:
            for i in range(numatoms):
                intersect[i] = list(set(intersect[i]).intersection(intersect_axis[i]))
    if strict:
        ibox = numpy.linalg.inv(p.box)
        for i in range(numatoms):
            l = intersect[i][:]
            for j in l:
                dist = numpy.linalg.norm(pbc(p.r[i] - p.r[j], p.box, ibox))
                if dist > cutoff:
                    intersect[i].remove(j)
                    intersect[j].remove(i)
    return intersect
    

    
def coordination_numbers(p, cutoff):
    """ Returns a list of coordination numbers for each atom in p """
    nl = sweep_and_prune(p, cutoff)
    return [len(l) for l in nl]
    


def least_coordinated(p, cutoff):
    """ Returns a list of atom indices in p with the lowest coordination numbers for unfrozen atoms"""
    cn = coordination_numbers(p, cutoff)
    maxcoord = max(cn)
    mincoord = min(cn)
    while mincoord <= maxcoord:
        least = []
        for i in range(len(cn)):
            if cn[i] <= mincoord and p.free[i]:
                least.append(i)
        if len(least) > 0:
            return least
        mincoord += 1

    
    
    
