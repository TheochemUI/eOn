
""" The atoms module. """
from eon.config import config

from math import cos, sin, acos
import numpy
import logging
logger = logging.getLogger('atoms')

class Atoms:
    """ The Atoms class. """

    def __init__(self, n_atoms):
        self.r = numpy.zeros((n_atoms,3))
        self.free = numpy.ones(n_atoms)
        self.box = numpy.zeros((3,3))
        self.names = ['']*n_atoms
        self.mass = numpy.zeros(n_atoms)

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
        p.mass = self.mass.copy()
        return p

    def free_r(self):
        nfree = sum(self.free)
        temp = numpy.zeros((nfree, 3))
        index = 0
        for i in range(len(self.r)):
            if self.free[i]:
                temp[index] = self.r[i]
                index += 1
        return temp

    def append(self, r, free, name, mass):
        self.r = numpy.append(self.r, [r], 0)
        self.free = numpy.append(self.free, free)
        self.names.append(name)
        self.mass = numpy.append(self.mass, mass)

def pbc(r, box, ibox = None):
    """
    Applies periodic boundary conditions.
    Parameters:
        r:      the vector the boundary conditions are applied to
        box:    the box that defines the boundary conditions
        ibox:   the inverse of the box. This will be calcluated if not provided.
    """
    if ibox is None:
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
    return numpy.sqrt(numpy.sum(diff**2.0, axis=1))

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

def get_process_atoms(r, p, epsilon_r=0.2, nshells=1):
    '''
    Given the reactant and product configurations of a process, return
    the atoms that move significantly and their neighbors along the trajectory.
    '''
    mobileAtoms = []
    r2p = per_atom_norm(p.r - r.r, r.box)
    for i in range(len(r)):
        if r2p[i] > epsilon_r:
            mobileAtoms.append(i)
    if len(mobileAtoms) == 0:
        mobileAtoms.append(list(r2p).index(max(r2p)))
    neighborAtoms = []
    for atom in mobileAtoms:
        r1 = elements[r.names[atom]]["radius"]
        for i in range(len(r)):
            if i in mobileAtoms or i in neighborAtoms:
                continue
            r2 = elements[r.names[i]]["radius"]
            maxDist = (r1 + r2) * (1.0 + 0.2)*nshells
            if numpy.linalg.norm(pbc(p.r[atom] - p.r[i], p.box)) < maxDist:
                neighborAtoms.append(i)
    return mobileAtoms + neighborAtoms

def identical(atoms1, atoms2):
    '''
    Determines whether two structures are identical if atoms of the same
    element are considered indistinguishable.
           atoms1:  first atoms object for comparison
           atoms2:  second atoms object for comparison
        epsilon_r:  distance (in angstroms) that two atoms must be seperated by
                    in order to be considered different
    '''
    #XXX: n^2
    epsilon_r = config.comp_eps_r
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
        elif atoms1.names[i] != atoms2.names[i]:
            return False

    for i in mismatch:
        pan = per_atom_norm(atoms1.r - atoms2.r[i], box, ibox)
        minpan = 1e300
        minj = 0
        for j in range(len(pan)):
            if i == j:
                continue
            if pan[j] < minpan:
                minpan = pan[j]
                minj = j
        if not (minpan < epsilon_r and atoms1.names[minj] == atoms2.names[i]):
            return False
    return True


def brute_neighbor_list(p, cutoff):
    nl = []
    ibox = numpy.linalg.inv(p.box)
    for a in range(len(p)):
        nl.append([])
        for b in range(len(p)):
            if b != a:
                dist = numpy.linalg.norm(pbc(p.r[a] - p.r[b], p.box, ibox))
                if dist < cutoff:
                    nl[a].append(b)
    return nl


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

def neighbor_list(p, cutoff, brute=False):
    if brute:
        nl = brute_neighbor_list(p, cutoff)
    else:
        nl = sweep_and_prune(p, cutoff)
    return nl

def neighbor_list_vectors(p, cutoff, brute=False):
    '''Points from center to neighbor'''
    nl = neighbor_list(p, cutoff, brute=brute)
    nl_vec = []
    
    ibox = numpy.linalg.inv(p.box)
    for center_index in range(len(p)):
        nl_vec.append([])
        for neighbor_index in nl[center_index]:
            vec = pbc(p.r[neighbor_index] - p.r[center_index], p.box, ibox)
            nl_vec[center_index].append(vec)
    return nl_vec

def coordination_numbers(p, cutoff, brute=False):
    """ Returns a list of coordination numbers for each atom in p """
    nl = neighbor_list(p, cutoff, brute)
    return [len(l) for l in nl]

def least_coordinated(p, cutoff, brute=False):
    """ Returns a list of atom indices in p with the lowest coordination numbers
        for unfrozen atoms"""
    cn = coordination_numbers(p, cutoff, brute)
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


def match(a,b,eps_r,neighbor_cutoff,indistinguishable):
    if len(a)!=len(b):
        return False

    if config.comp_check_rotation:
        if indistinguishable and config.comp_use_identical:
            return get_mappings(a,b,eps_r,neighbor_cutoff)
        else:
            return rot_match(a,b)
    else:
        if indistinguishable and config.comp_use_identical:
            return identical(a,b)
        else:
            diff = pbc(a.r-b.r, a.box)
            return numpy.max(numpy.sum(diff**2.0, axis=1)) < eps_r**2.0
            #return max(per_atom_norm(a.r-b.r, a.box))<eps_r


def point_energy_match(file_a, energy_a, file_b, energy_b):
    import eon.fileio as io
    if abs(energy_a - energy_b) > config.comp_eps_e:
        return False
    a = io.loadcon(file_a)
    b = io.loadcon(file_b)
    if match(a, b, config.comp_eps_r, config.comp_neighbor_cutoff, False):
        return True

def points_energies_match(file_a, energy_a, files_b, energies_b):
    for i in range(len(files_b)):
        if point_energy_match(file_a, energy_a, files_b[i], energies_b[i]):
            return i
    return None


def rot_match(a, b):
    if not (a.free.all() and b.free.all()):
        logger.warning("Comparing structures with frozen atoms with rotational matching; check_rotation may be set incorrectly")
    acm = sum(a.r)/len(a)
    bcm = sum(b.r)/len(b)

    # SVD algorithm ("Least-Squares Fitting of Two 3-D Point Sets" by Arun, Huang, and Blostein)
    ta_r = a.r - acm
    tb_r = b.r - bcm
    H = numpy.dot(ta_r.transpose() , tb_r)
    U, S, V = numpy.linalg.svd(H)
    R = numpy.dot(U, V)
    if numpy.linalg.det(R) < 0:
        V[2] *= -1
        R = numpy.dot(U, V)
    ta_r = numpy.dot(ta_r, R)
    dist = max(numpy.linalg.norm(ta_r - tb_r, axis=1))
    return dist < config.comp_eps_r

    ''' Quaternion algorithm
    ta = a.copy()
    tb = b.copy()
    ta.r -= acm
    tb.r -= bcm

    #Horn, J. Opt. Soc. Am. A, 1987
    m = numpy.dot(tb.r.transpose(), ta.r)
    sxx = m[0][0]
    sxy = m[0][1]
    sxz = m[0][2]
    syx = m[1][0]
    syy = m[1][1]
    syz = m[1][2]
    szx = m[2][0]
    szy = m[2][1]
    szz = m[2][2]

    n = numpy.zeros((4,4))
    n[0][1] = syz-szy
    n[0][2] = szx-sxz
    n[0][3] = sxy-syx

    n[1][2] = sxy+syx
    n[1][3] = szx+sxz

    n[2][3] = syz + szy

    n += n.transpose()

    n[0][0] = sxx + syy + szz
    n[1][1] = sxx-syy-szz
    n[2][2] = -sxx + syy -szz
    n[3][3] = -sxx -syy + szz

    w,v = numpy.linalg.eig(n)
    maxw = 0
    maxv = 0
    for i in range(len(w)):
        if w[i] > maxw:
            maxw = w[i]
            maxv = v[:,i]

    R = numpy.zeros((3,3))

    aa = maxv[0]**2
    bb = maxv[1]**2
    cc = maxv[2]**2
    dd = maxv[3]**2
    ab = maxv[0]*maxv[1]
    ac = maxv[0]*maxv[2]
    ad = maxv[0]*maxv[3]
    bc = maxv[1]*maxv[2]
    bd = maxv[1]*maxv[3]
    cd = maxv[2]*maxv[3]

    R[0][0] = aa + bb - cc - dd
    R[0][1] = 2*(bc-ad)
    R[0][2] = 2*(bd+ac)
    R[1][0] = 2*(bc+ad)
    R[1][1] = aa - bb + cc - dd
    R[1][2] = 2*(cd-ab)
    R[2][0] = 2*(bd-ac)
    R[2][1] = 2*(cd+ab)
    R[2][2] = aa - bb - cc + dd
    tb.r = numpy.dot(tb.r, R.transpose())

    dist = max(per_atom_norm(ta.r - tb.r, ta.box))
    return dist < config.comp_eps_r
    '''

    ### This gives the RMSD faster, but does not give the optimial rotation
    ### this could be amended by solving for the eigenvector corresponding to the largest eigenvalue
    ### Theobald, Acta Crystallographica A, 2005
    #
    ##could be faster if done explicitly
    #c0 = numpy.linalg.det(k)
    #c1 = -8*numpy.linalg.det(m)
    #c2 = -2*numpy.trace(numpy.dot(m.transpose(), m))

    #ga = numpy.trace(numpy.dot(ta.r.transpose(), ta.r))
    #gb = numpy.trace(numpy.dot(tb.r.transpose(), tb.r))
    #
    #lold = 0.0
    #l = (ga + gb)/2.0
    #while abs(lold - l) > 0.00001:
    #    lold = l
    #    l -= (l**4 + c2*l**2 + c1*l + c0)/(4*l**3 + 2*c2*l + c1)
    #rmsd = sqrt((ga + gb - 2*l)/len(a))
    #return rmsd < config.comp_rot_rmsd


def rotm(axis, theta):
    '''
    Gives the matrix representing a rotation of theta radians about axis
    '''
    u = axis[0]
    v = axis[1]
    w = axis[2]
    u2 = u*u
    v2 = v*v
    w2 = w*w
    ct = cos(theta)
    st = sin(theta)
    mag = numpy.linalg.norm(axis)
    if (mag*mag == 0 or theta == 0.0):
        return numpy.identity(3)
    return numpy.array([
        [u2 +(v2 +w2)*ct, u*v*(1-ct)-w*mag*st, u*w*(1-ct)+v*mag*st],
        [u*v*(1-ct)+w*mag*st, v2 +(u2 +w2)*ct, v*w*(1-ct)-u*mag*st],
        [u*w*(1-ct)-v*mag*st, v*w*(1-ct)+u*mag*st, w2 +(v2 +u2)*ct]
        ])/(mag*mag)


def cna(p, cutoff, brute=False):
    """ Returns a list of cna numbers for all atoms in p
        Inspired by the CNA code provided by Asap (wiki.fysik.dtu.dk/asap)"""
    can_values = numpy.zeros(len(p))
    nr_FCC = numpy.zeros(len(p))
    nr_HCP = numpy.zeros(len(p))
    nl = neighbor_list(p, cutoff, brute)

    # loops over all the atoms
    for a2 in range(len(p)):
        nl_a2 = nl[a2]
        # loops over the atoms neighboring a2
        for n2 in range(len(nl_a2)):
            a1 = nl_a2[n2];
            if a1 < a2:
                common = []
                nl_a1 = nl[a1]
                # loops over the atoms neighboring a1
                for n1 in range(len(nl_a1)):
                    a3 = nl_a1[n1]
                    # checks if atom a_3 is a common neighbor to a1 and a2
                    for m2 in range(len(nl_a2)):
                        if a3 == nl_a2[m2]:
                            common.append(a3)
                # determines the connectivity of common neighbors
                if len(common) == 4:
                    bonds_nr = 0
                    bonds_sum = 0
                    for j2 in range(1,4):
                        nl_j2 = nl[common[j2]]
                        for j1 in range(j2):
                            for n in range(len(nl_j2)):
                                if common[j1] == nl_j2[n]:
                                    bonds_nr += 1
                                    bonds_sum += j1 + j2

                    if bonds_nr == 2:
                        if bonds_sum == 6:
                            nr_FCC[a1] += 1
                            nr_FCC[a2] += 1
                        else:
                            nr_HCP[a1] += 1
                            nr_HCP[a2] += 1

    # 2: fcc (421), 1: hcp (422), 0: other
    for i in range(len(p)):
        if len(nl[i]) == 12:
            if nr_FCC[i] == 12:
                can_values[i] = 2
            elif (nr_FCC[i] == 6) and (nr_HCP[i] == 6):
                can_values[i] = 1
    return can_values

def not_HCP_or_FCC(p, cutoff, brute=False):
    """ Returns a list of indices for the atoms with cna = 0 """
    not_cna = []
    cna_numbers = cna(p, cutoff, brute)
    for i in range(len(cna_numbers)):
        if cna_numbers[i] == 0:
            not_cna.append(i)
    return not_cna


# ### TShacked start
def cnat(p, cutoff, brute=False):
    """ Returns a list of cna numbers for all atoms in p
        Inspired by the CNA code provided by Asap (wiki.fysik.dtu.dk/asap)"""
    can_values = numpy.zeros(len(p))
    nr_5 = numpy.zeros(len(p))
    nr_6 = numpy.zeros(len(p))
    nl = neighbor_list(p, cutoff, brute)

    # loops over all the atoms
    for a2 in range(len(p)):
        nl_a2 = nl[a2]
        # loops over the atoms neighboring a2
        for n2 in range(len(nl_a2)):
            a1 = nl_a2[n2];
            if a1 < a2:
                common = []
                nl_a1 = nl[a1]
                # loops over the atoms neighboring a1
                for n1 in range(len(nl_a1)):
                    a3 = nl_a1[n1]
                    # checks if atom a_3 is a common neighbor to a1 and a2
                    for m2 in range(len(nl_a2)):
                        if a3 == nl_a2[m2]:
                            common.append(a3)
                # determines the connectivity of common neighbors
                #print a2, n2, len(common)
                if len(common) in [5,6]:
                    bonds_nr = 0
                    bonds_sum = 0
                    for j2 in range(1,len(common)):
                        nl_j2 = nl[common[j2]]
                        for j1 in range(j2):
                            for n in range(len(nl_j2)):
                                if common[j1] == nl_j2[n]:
                                    bonds_nr += 1
                                    bonds_sum += j1 + j2

                    if bonds_nr == 5 and len(common) == 5:
                        nr_5[a1] += 1
                        nr_5[a2] += 1
                    elif bonds_nr == 6 and len(common) == 6:
                        nr_6[a1] += 1
                        nr_6[a2] += 1

    # 1: CN12, 2: CN14, 3: CN15, 4: CN16, 5: BCC, 0: other
    for i in range(len(p)):
        if nr_5[i] == 12:
            if len(nl[i]) == 12:
                can_values[i] = 1
            if (len(nl[i]) == 14) and (nr_6[i] == 2):
                can_values[i] = 2
            if (len(nl[i]) == 15) and (nr_6[i] == 3):
                can_values[i] = 3
            if (len(nl[i]) == 16) and (nr_6[i] == 4):
                can_values[i] = 4
        if nr_5[i] == 0:
            if nr_6[i] == 8 and len(nl[i]) == 14:
                can_values[i] = 5
    return can_values

def cnar(p, cutoff, brute=False):
    """ Returns a list of cna numbers for all atoms in p
        Inspired by the CNA code provided by Asap (wiki.fysik.dtu.dk/asap)"""
    # not compatible with older python versions
    #cna = {i:{} for i in range(len(p))}

    # make a dict of dicts for each atom.
    cna = {}
    for i in range(len(p)):
        cna[i] = {}

    nl = neighbor_list(p, cutoff, brute)

    def codeString(j,k,l):
        return "%d,%d,%d" % (j,k,l)

    # loops over all the atoms
    for a2 in range(len(p)):
        nl_a2 = nl[a2]
        # loops over the atoms neighboring a2
        for n2 in range(len(nl_a2)):
            a1 = nl_a2[n2];
            if a1 < a2: # prevent double counting?
                common = []
                nl_a1 = nl[a1]
                # loops over the atoms neighboring a1
                for n1 in range(len(nl_a1)):
                    a3 = nl_a1[n1]
                    # checks if atom a_3 is a common neighbor to a1 and a2
                    for m2 in range(len(nl_a2)):
                        if a3 == nl_a2[m2]:
                            common.append(a3)
                # determines the connectivity of common neighbors
                bonds_nr = 0
                bonds_sum = 0
                for j2 in range(1,len(common)):
                    nl_j2 = nl[common[j2]]
                    for j1 in range(j2):
                        for n in range(len(nl_j2)):
                            if common[j1] == nl_j2[n]:
                                bonds_nr += 1
                                bonds_sum += j1 + j2
                code = codeString(len(common), bonds_nr, bonds_sum)
                if not code in cna[a1]:
                    cna[a1][code] = 0
                cna[a1][code] += 1
                if not code in cna[a2]:
                    cna[a2][code] = 0
                cna[a2][code] += 1

    return cna


def not_TCP(p, cutoff, brute=False):
    """ Returns a list of indices for the atoms with cna = 0 """
    not_cna = []
    cna_numbers = cnat(p, cutoff, brute)
    for i in range(len(cna_numbers)):
        if cna_numbers[i] == 0 or cna_numbers[i] == 5:
            not_cna.append(i)
    return not_cna

def not_TCP_or_BCC(p, cutoff, brute=False):
    """ Returns a list of indices for the atoms with cna = 0 """
    not_cna = []
    cna_numbers = cnat(p, cutoff, brute)
    for i in range(len(cna_numbers)):
        if cna_numbers[i] == 0:
            not_cna.append(i)
    return not_cna


import sys
sys.setrecursionlimit(10000)
def get_mappings(a, b, eps_r, neighbor_cutoff, mappings = None):
    """ A recursive depth-first search for a complete set of mappings from atoms
        in configuration a to atoms in configuration b. Do not use the mappings
        argument, this is only used internally for recursion.

        Returns None if no mapping was found, or a dictionary mapping atom
        indices a to atom indices b.

        Note: If a and b are mirror images, this function will still return a
        mapping from a to b, even though it may not be possible to align them
        through translation and rotation. """
    # If this is the top-level user call, create and loop through top-level
    # mappings.
    if mappings is None:
        # Find the least common coordination number in b.
        bCoordinations = coordination_numbers(b, neighbor_cutoff)
        bCoordinationsCounts = {}
        for coordination in bCoordinations:
            if coordination in bCoordinationsCounts:
                bCoordinationsCounts[coordination] += 1
            else:
                bCoordinationsCounts[coordination] = 1
        bLeastCommonCoordination = list(bCoordinationsCounts.keys())[0]
        for coordination in list(bCoordinationsCounts.keys()):
            if bCoordinationsCounts[coordination] < bCoordinationsCounts[bLeastCommonCoordination]:
                bLeastCommonCoordination = coordination
        # Find one atom in a with the least common coordination number in b.
        # If it does not exist, return None.
        aCoordinations = coordination_numbers(a, neighbor_cutoff)
        try:
            aAtom = aCoordinations.index(bLeastCommonCoordination)
        except ValueError:
            return None
        # Create a mapping from the atom chosen from a to each of the atoms with
        # the least common coordination number in b, and recurse.
        for i in range(len(bCoordinations)):
            if bCoordinations[i] == bLeastCommonCoordination:
                # Make sure the element types are the same.
                if a.names[aAtom] != b.names[i]:
                    continue
                mappings = get_mappings(a, b, eps_r, neighbor_cutoff, {aAtom:i})
                # If the result is not none, then we found a successful mapping.
                if mappings is not None:
                    return mappings
        # There were no mappings.
        return None

    # This is a recursed invocation of this function.
    else:
        # Find an atom from a that has not yet been mapped.
        unmappedA = 0
        while unmappedA < len(a):
            if unmappedA not in list(mappings.keys()):
                break
            unmappedA += 1
        # Calculate the distances from unmappedA to all mapped a atoms.
        distances = {}
        for i in list(mappings.keys()):
            distances[i] = numpy.linalg.norm(pbc(a.r[unmappedA] - a.r[i], a.box))
        # Loop over each unmapped b atom. Compare the distances between it and
        # the mapped b atoms to the corresponding distances between unmappedA
        # and the mapped atoms. If everything is similar, create a new mapping
        # and recurse.
        for bAtom in range(len(b)):
            if bAtom not in list(mappings.values()):
                for aAtom in distances:
                    # Break if type check fails.
                    if b.names[bAtom] != a.names[unmappedA]:
                        break
                    # Break if distance check fails
                    bDist = numpy.linalg.norm(pbc(b.r[bAtom] - b.r[mappings[aAtom]], b.box))
                    if abs(distances[aAtom] - bDist) > eps_r:
                        break
                else:
                    # All distances were good, so create a new mapping.
                    newMappings = mappings.copy()
                    newMappings[unmappedA] = bAtom
                    # If this is now a complete mapping from a to b, return it.
                    if len(newMappings) == len(a):
                        return newMappings
                    # Otherwise, recurse.
                    newMappings = get_mappings(a, b, eps_r, neighbor_cutoff, newMappings)
                    # Pass any successful mapping up the recursion chain.
                    if newMappings is not None:
                        return newMappings
        # There were no mappings.
        return None

def get_rotation_matrix(axis, theta):
    axis = axis / numpy.linalg.norm(axis)
    t = theta
    T = 1.0 - cos(t)
    rx, ry, rz = axis
    rotmat = numpy.zeros((3, 3))
    rotmat[0][0] = T*rx*rx + cos(t)
    rotmat[0][1] = T*ry*rx + rz*sin(t)
    rotmat[0][2] = T*rz*rx - ry*sin(t)
    rotmat[1][0] = T*rx*ry - rz*sin(t)
    rotmat[1][1] = T*ry*ry + cos(t)
    rotmat[1][2] = T*rz*ry + rx*sin(t)
    rotmat[2][0] = T*rx*rz + ry*sin(t)
    rotmat[2][1] = T*ry*rz - rx*sin(t)
    rotmat[2][2] = T*rz*rz + cos(t)
    return rotmat

def rotate(r, axis, center, angle):
    new_r = r.copy()
    if abs(angle) == 0.0:
        return new_r
    rotmat = get_rotation_matrix(axis, angle)
    center = center.copy()
    new_r -= center
    new_r = numpy.dot(new_r, rotmat)
    new_r += center
    return new_r


def internal_motion(a, b):
    """ Takes two atoms objects and returns the motion from a to b that is
    entirely internal - no rotation or translation, in the form of a new atoms
    object. """
    b = b.copy()
    b.r -= a.r[0] - b.r[0]
    a0a1 = (a.r[1] - a.r[0]) / numpy.linalg.norm(a.r[1] - a.r[0])
    b0b1 = (b.r[1] - b.r[0]) / numpy.linalg.norm(b.r[1] - b.r[0])
    axis1 = numpy.cross(b0b1, a0a1) / numpy.linalg.norm(numpy.cross(b0b1, a0a1))
    theta1 = acos((a0a1*b0b1).sum())
    b.r = rotate(b.r, axis1, a.r[0], theta1)
    axis2 = (a.r[2] - a.r[0]) / numpy.linalg.norm(a.r[2] - a.r[0])
    va = a.r[2] - ((a.r[2] - a.r[0]) * axis2).sum() * axis2
    va = va / numpy.linalg.norm(va)
    vb = b.r[2] - ((b.r[2] - a.r[0]) * axis2).sum() * axis2
    vb = vb / numpy.linalg.norm(vb)
    theta2 = acos((va * vb).sum())
    b.r = rotate(b.r, axis2, a.r[0], theta2)
    return b


elements = {}
numElements = 119
elements[  0] = elements[ 'Xx'] = {'symbol':  'Xx', 'name':       'unknown', 'mass':   1.00000000, 'radius':  1.0000, 'color': [1.000, 0.078, 0.576], 'number': 0}
elements[  1] = elements[  'H'] = {'symbol':   'H', 'name':      'hydrogen', 'mass':   1.00794000, 'radius':  0.3100, 'color': [1.000, 1.000, 1.000], 'number': 1}
elements[  2] = elements[ 'He'] = {'symbol':  'He', 'name':        'helium', 'mass':   4.00260200, 'radius':  0.2800, 'color': [0.851, 1.000, 1.000], 'number': 2}
elements[  3] = elements[ 'Li'] = {'symbol':  'Li', 'name':       'lithium', 'mass':   6.94100000, 'radius':  1.2800, 'color': [0.800, 0.502, 1.000], 'number': 3}
elements[  4] = elements[ 'Be'] = {'symbol':  'Be', 'name':     'beryllium', 'mass':   9.01218200, 'radius':  0.9600, 'color': [0.761, 1.000, 0.000], 'number': 4}
elements[  5] = elements[  'B'] = {'symbol':   'B', 'name':         'boron', 'mass':  10.81100000, 'radius':  0.8400, 'color': [1.000, 0.710, 0.710], 'number': 5}
elements[  6] = elements[  'C'] = {'symbol':   'C', 'name':        'carbon', 'mass':  12.01070000, 'radius':  0.7300, 'color': [0.565, 0.565, 0.565], 'number': 6}
elements[  7] = elements[  'N'] = {'symbol':   'N', 'name':      'nitrogen', 'mass':  14.00670000, 'radius':  0.7100, 'color': [0.188, 0.314, 0.973], 'number': 7}
elements[  8] = elements[  'O'] = {'symbol':   'O', 'name':        'oxygen', 'mass':  15.99940000, 'radius':  0.6600, 'color': [1.000, 0.051, 0.051], 'number': 8}
elements[  9] = elements[  'F'] = {'symbol':   'F', 'name':      'fluorine', 'mass':  18.99840320, 'radius':  0.5700, 'color': [0.565, 0.878, 0.314], 'number': 9}
elements[ 10] = elements[ 'Ne'] = {'symbol':  'Ne', 'name':          'neon', 'mass':  20.17970000, 'radius':  0.5800, 'color': [0.702, 0.890, 0.961], 'number': 10}
elements[ 11] = elements[ 'Na'] = {'symbol':  'Na', 'name':        'sodium', 'mass':  22.98976928, 'radius':  1.6600, 'color': [0.671, 0.361, 0.949], 'number': 11}
elements[ 12] = elements[ 'Mg'] = {'symbol':  'Mg', 'name':     'magnesium', 'mass':  24.30500000, 'radius':  1.4100, 'color': [0.541, 1.000, 0.000], 'number': 12}
elements[ 13] = elements[ 'Al'] = {'symbol':  'Al', 'name':      'aluminum', 'mass':  26.98153860, 'radius':  1.2100, 'color': [0.749, 0.651, 0.651], 'number': 13}
elements[ 14] = elements[ 'Si'] = {'symbol':  'Si', 'name':       'silicon', 'mass':  28.08550000, 'radius':  1.1100, 'color': [0.941, 0.784, 0.627], 'number': 14}
elements[ 15] = elements[  'P'] = {'symbol':   'P', 'name':    'phosphorus', 'mass':  30.97376200, 'radius':  1.0700, 'color': [1.000, 0.502, 0.000], 'number': 15}
elements[ 16] = elements[  'S'] = {'symbol':   'S', 'name':        'sulfur', 'mass':  32.06500000, 'radius':  1.0500, 'color': [1.000, 1.000, 0.188], 'number': 16}
elements[ 17] = elements[ 'Cl'] = {'symbol':  'Cl', 'name':      'chlorine', 'mass':  35.45300000, 'radius':  1.0200, 'color': [0.122, 0.941, 0.122], 'number': 17}
elements[ 18] = elements[ 'Ar'] = {'symbol':  'Ar', 'name':         'argon', 'mass':  39.94800000, 'radius':  1.0600, 'color': [0.502, 0.820, 0.890], 'number': 18}
elements[ 19] = elements[  'K'] = {'symbol':   'K', 'name':     'potassium', 'mass':  39.09830000, 'radius':  2.0300, 'color': [0.561, 0.251, 0.831], 'number': 19}
elements[ 20] = elements[ 'Ca'] = {'symbol':  'Ca', 'name':       'calcium', 'mass':  40.07800000, 'radius':  1.7600, 'color': [0.239, 1.000, 0.000], 'number': 20}
elements[ 21] = elements[ 'Sc'] = {'symbol':  'Sc', 'name':      'scandium', 'mass':  44.95591200, 'radius':  1.7000, 'color': [0.902, 0.902, 0.902], 'number': 21}
elements[ 22] = elements[ 'Ti'] = {'symbol':  'Ti', 'name':      'titanium', 'mass':  47.86700000, 'radius':  1.6000, 'color': [0.749, 0.761, 0.780], 'number': 22}
elements[ 23] = elements[  'V'] = {'symbol':   'V', 'name':      'vanadium', 'mass':  50.94150000, 'radius':  1.5300, 'color': [0.651, 0.651, 0.671], 'number': 23}
elements[ 24] = elements[ 'Cr'] = {'symbol':  'Cr', 'name':      'chromium', 'mass':  51.99610000, 'radius':  1.3900, 'color': [0.541, 0.600, 0.780], 'number': 24}
elements[ 25] = elements[ 'Mn'] = {'symbol':  'Mn', 'name':     'manganese', 'mass':  54.93804500, 'radius':  1.3900, 'color': [0.611, 0.478, 0.780], 'number': 25}
elements[ 26] = elements[ 'Fe'] = {'symbol':  'Fe', 'name':          'iron', 'mass':  55.84500000, 'radius':  1.3200, 'color': [0.878, 0.400, 0.200], 'number': 26}
elements[ 27] = elements[ 'Co'] = {'symbol':  'Co', 'name':        'cobalt', 'mass':  58.69340000, 'radius':  1.2600, 'color': [0.941, 0.565, 0.627], 'number': 27}
elements[ 28] = elements[ 'Ni'] = {'symbol':  'Ni', 'name':        'nickel', 'mass':  58.93319500, 'radius':  1.2400, 'color': [0.314, 0.816, 0.314], 'number': 28}
elements[ 29] = elements[ 'Cu'] = {'symbol':  'Cu', 'name':        'copper', 'mass':  63.54600000, 'radius':  1.3200, 'color': [0.784, 0.502, 0.200], 'number': 29}
elements[ 30] = elements[ 'Zn'] = {'symbol':  'Zn', 'name':          'zinc', 'mass':  65.38000000, 'radius':  1.2200, 'color': [0.490, 0.502, 0.690], 'number': 30}
elements[ 31] = elements[ 'Ga'] = {'symbol':  'Ga', 'name':       'gallium', 'mass':  69.72300000, 'radius':  1.2200, 'color': [0.761, 0.561, 0.561], 'number': 31}
elements[ 32] = elements[ 'Ge'] = {'symbol':  'Ge', 'name':     'germanium', 'mass':  72.64000000, 'radius':  1.2000, 'color': [0.400, 0.561, 0.561], 'number': 32}
elements[ 33] = elements[ 'As'] = {'symbol':  'As', 'name':       'arsenic', 'mass':  74.92160000, 'radius':  1.1900, 'color': [0.741, 0.502, 0.890], 'number': 33}
elements[ 34] = elements[ 'Se'] = {'symbol':  'Se', 'name':      'selenium', 'mass':  78.96000000, 'radius':  1.2000, 'color': [1.000, 0.631, 0.000], 'number': 34}
elements[ 35] = elements[ 'Br'] = {'symbol':  'Br', 'name':       'bromine', 'mass':  79.90400000, 'radius':  1.2000, 'color': [0.651, 0.161, 0.161], 'number': 35}
elements[ 36] = elements[ 'Kr'] = {'symbol':  'Kr', 'name':       'krypton', 'mass':  83.79800000, 'radius':  1.1600, 'color': [0.361, 0.722, 0.820], 'number': 36}
elements[ 37] = elements[ 'Rb'] = {'symbol':  'Rb', 'name':      'rubidium', 'mass':  85.46780000, 'radius':  2.2000, 'color': [0.439, 0.180, 0.690], 'number': 37}
elements[ 38] = elements[ 'Sr'] = {'symbol':  'Sr', 'name':     'strontium', 'mass':  87.62000000, 'radius':  1.9500, 'color': [0.000, 1.000, 0.000], 'number': 38}
elements[ 39] = elements[  'Y'] = {'symbol':   'Y', 'name':       'yttrium', 'mass':  88.90585000, 'radius':  1.9000, 'color': [0.580, 1.000, 1.000], 'number': 39}
elements[ 40] = elements[ 'Zr'] = {'symbol':  'Zr', 'name':     'zirconium', 'mass':  91.22400000, 'radius':  1.7500, 'color': [0.580, 0.878, 0.878], 'number': 40}
elements[ 41] = elements[ 'Nb'] = {'symbol':  'Nb', 'name':       'niobium', 'mass':  92.90638000, 'radius':  1.6400, 'color': [0.451, 0.761, 0.788], 'number': 41}
elements[ 42] = elements[ 'Mo'] = {'symbol':  'Mo', 'name':    'molybdenum', 'mass':  95.96000000, 'radius':  1.5400, 'color': [0.329, 0.710, 0.710], 'number': 42}
elements[ 43] = elements[ 'Tc'] = {'symbol':  'Tc', 'name':    'technetium', 'mass':  98.00000000, 'radius':  1.4700, 'color': [0.231, 0.620, 0.620], 'number': 43}
elements[ 44] = elements[ 'Ru'] = {'symbol':  'Ru', 'name':     'ruthenium', 'mass': 101.07000000, 'radius':  1.4600, 'color': [0.141, 0.561, 0.561], 'number': 44}
elements[ 45] = elements[ 'Rh'] = {'symbol':  'Rh', 'name':       'rhodium', 'mass': 102.90550000, 'radius':  1.4200, 'color': [0.039, 0.490, 0.549], 'number': 45}
elements[ 46] = elements[ 'Pd'] = {'symbol':  'Pd', 'name':     'palladium', 'mass': 106.42000000, 'radius':  1.3900, 'color': [0.000, 0.412, 0.522], 'number': 46}
elements[ 47] = elements[ 'Ag'] = {'symbol':  'Ag', 'name':        'silver', 'mass': 107.86820000, 'radius':  1.4500, 'color': [0.753, 0.753, 0.753], 'number': 47}
elements[ 48] = elements[ 'Cd'] = {'symbol':  'Cd', 'name':       'cadmium', 'mass': 112.41100000, 'radius':  1.4400, 'color': [1.000, 0.851, 0.561], 'number': 48}
elements[ 49] = elements[ 'In'] = {'symbol':  'In', 'name':        'indium', 'mass': 114.81800000, 'radius':  1.4200, 'color': [0.651, 0.459, 0.451], 'number': 49}
elements[ 50] = elements[ 'Sn'] = {'symbol':  'Sn', 'name':           'tin', 'mass': 118.71000000, 'radius':  1.3900, 'color': [0.400, 0.502, 0.502], 'number': 50}
elements[ 51] = elements[ 'Sb'] = {'symbol':  'Sb', 'name':      'antimony', 'mass': 121.76000000, 'radius':  1.3900, 'color': [0.620, 0.388, 0.710], 'number': 51}
elements[ 52] = elements[ 'Te'] = {'symbol':  'Te', 'name':     'tellurium', 'mass': 127.60000000, 'radius':  1.3800, 'color': [0.831, 0.478, 0.000], 'number': 52}
elements[ 53] = elements[  'I'] = {'symbol':   'I', 'name':        'iodine', 'mass': 126.90470000, 'radius':  1.3900, 'color': [0.580, 0.000, 0.580], 'number': 53}
elements[ 54] = elements[ 'Xe'] = {'symbol':  'Xe', 'name':         'xenon', 'mass': 131.29300000, 'radius':  1.4000, 'color': [0.259, 0.620, 0.690], 'number': 54}
elements[ 55] = elements[ 'Cs'] = {'symbol':  'Cs', 'name':        'cesium', 'mass': 132.90545190, 'radius':  2.4400, 'color': [0.341, 0.090, 0.561], 'number': 55}
elements[ 56] = elements[ 'Ba'] = {'symbol':  'Ba', 'name':        'barium', 'mass': 137.32700000, 'radius':  2.1500, 'color': [0.000, 0.788, 0.000], 'number': 56}
elements[ 57] = elements[ 'La'] = {'symbol':  'La', 'name':     'lanthanum', 'mass': 138.90547000, 'radius':  2.0700, 'color': [0.439, 0.831, 1.000], 'number': 57}
elements[ 58] = elements[ 'Ce'] = {'symbol':  'Ce', 'name':        'cerium', 'mass': 140.11600000, 'radius':  2.0400, 'color': [1.000, 1.000, 0.780], 'number': 58}
elements[ 59] = elements[ 'Pr'] = {'symbol':  'Pr', 'name':  'praseodymium', 'mass': 140.90765000, 'radius':  2.0300, 'color': [0.851, 1.000, 0.780], 'number': 59}
elements[ 60] = elements[ 'Nd'] = {'symbol':  'Nd', 'name':     'neodymium', 'mass': 144.24200000, 'radius':  2.0100, 'color': [0.780, 1.000, 0.780], 'number': 60}
elements[ 61] = elements[ 'Pm'] = {'symbol':  'Pm', 'name':    'promethium', 'mass': 145.00000000, 'radius':  1.9900, 'color': [0.639, 1.000, 0.780], 'number': 61}
elements[ 62] = elements[ 'Sm'] = {'symbol':  'Sm', 'name':      'samarium', 'mass': 150.36000000, 'radius':  1.9800, 'color': [0.561, 1.000, 0.780], 'number': 62}
elements[ 63] = elements[ 'Eu'] = {'symbol':  'Eu', 'name':      'europium', 'mass': 151.96400000, 'radius':  1.9800, 'color': [0.380, 1.000, 0.780], 'number': 63}
elements[ 64] = elements[ 'Gd'] = {'symbol':  'Gd', 'name':    'gadolinium', 'mass': 157.25000000, 'radius':  1.9600, 'color': [0.271, 1.000, 0.780], 'number': 64}
elements[ 65] = elements[ 'Tb'] = {'symbol':  'Tb', 'name':       'terbium', 'mass': 158.92535000, 'radius':  1.9400, 'color': [0.189, 1.000, 0.780], 'number': 65}
elements[ 66] = elements[ 'Dy'] = {'symbol':  'Dy', 'name':    'dysprosium', 'mass': 162.50000000, 'radius':  1.9200, 'color': [0.122, 1.000, 0.780], 'number': 66}
elements[ 67] = elements[ 'Ho'] = {'symbol':  'Ho', 'name':       'holmium', 'mass': 164.93032000, 'radius':  1.9200, 'color': [0.000, 1.000, 0.612], 'number': 67}
elements[ 68] = elements[ 'Er'] = {'symbol':  'Er', 'name':        'erbium', 'mass': 167.25900000, 'radius':  1.8900, 'color': [0.000, 0.902, 0.459], 'number': 68}
elements[ 69] = elements[ 'Tm'] = {'symbol':  'Tm', 'name':       'thulium', 'mass': 168.93421000, 'radius':  1.9000, 'color': [0.000, 0.831, 0.322], 'number': 69}
elements[ 70] = elements[ 'Yb'] = {'symbol':  'Yb', 'name':     'ytterbium', 'mass': 173.05400000, 'radius':  1.8700, 'color': [0.000, 0.749, 0.220], 'number': 70}
elements[ 71] = elements[ 'Lu'] = {'symbol':  'Lu', 'name':      'lutetium', 'mass': 174.96680000, 'radius':  1.8700, 'color': [0.000, 0.671, 0.141], 'number': 71}
elements[ 72] = elements[ 'Hf'] = {'symbol':  'Hf', 'name':       'hafnium', 'mass': 178.49000000, 'radius':  1.7500, 'color': [0.302, 0.761, 1.000], 'number': 72}
elements[ 73] = elements[ 'Ta'] = {'symbol':  'Ta', 'name':      'tantalum', 'mass': 180.94788000, 'radius':  1.7000, 'color': [0.302, 0.651, 1.000], 'number': 73}
elements[ 74] = elements[  'W'] = {'symbol':   'W', 'name':      'tungsten', 'mass': 183.84000000, 'radius':  1.6200, 'color': [0.129, 0.580, 0.839], 'number': 74}
elements[ 75] = elements[ 'Re'] = {'symbol':  'Re', 'name':       'rhenium', 'mass': 186.20700000, 'radius':  1.5100, 'color': [0.149, 0.490, 0.671], 'number': 75}
elements[ 76] = elements[ 'Os'] = {'symbol':  'Os', 'name':        'osmium', 'mass': 190.23000000, 'radius':  1.4400, 'color': [0.149, 0.400, 0.588], 'number': 76}
elements[ 77] = elements[ 'Ir'] = {'symbol':  'Ir', 'name':       'iridium', 'mass': 192.21700000, 'radius':  1.4100, 'color': [0.090, 0.329, 0.529], 'number': 77}
elements[ 78] = elements[ 'Pt'] = {'symbol':  'Pt', 'name':      'platinum', 'mass': 195.08400000, 'radius':  1.3600, 'color': [0.816, 0.816, 0.878], 'number': 78}
elements[ 79] = elements[ 'Au'] = {'symbol':  'Au', 'name':          'gold', 'mass': 196.96656900, 'radius':  1.3600, 'color': [1.000, 0.820, 0.137], 'number': 79}
elements[ 80] = elements[ 'Hg'] = {'symbol':  'Hg', 'name':       'mercury', 'mass': 200.59000000, 'radius':  1.3200, 'color': [0.722, 0.722, 0.816], 'number': 80}
elements[ 81] = elements[ 'Tl'] = {'symbol':  'Tl', 'name':      'thallium', 'mass': 204.38330000, 'radius':  1.4500, 'color': [0.651, 0.329, 0.302], 'number': 81}
elements[ 82] = elements[ 'Pb'] = {'symbol':  'Pb', 'name':          'lead', 'mass': 207.20000000, 'radius':  1.4600, 'color': [0.341, 0.349, 0.380], 'number': 82}
elements[ 83] = elements[ 'Bi'] = {'symbol':  'Bi', 'name':       'bismuth', 'mass': 208.98040000, 'radius':  1.4800, 'color': [0.620, 0.310, 0.710], 'number': 83}
elements[ 84] = elements[ 'Po'] = {'symbol':  'Po', 'name':      'polonium', 'mass': 210.00000000, 'radius':  1.4000, 'color': [0.671, 0.361, 0.000], 'number': 84}
elements[ 85] = elements[ 'At'] = {'symbol':  'At', 'name':      'astatine', 'mass': 210.00000000, 'radius':  1.5000, 'color': [0.459, 0.310, 0.271], 'number': 85}
elements[ 86] = elements[ 'Rn'] = {'symbol':  'Rn', 'name':         'radon', 'mass': 220.00000000, 'radius':  1.5000, 'color': [0.259, 0.510, 0.588], 'number': 86}
elements[ 87] = elements[ 'Fr'] = {'symbol':  'Fr', 'name':      'francium', 'mass': 223.00000000, 'radius':  2.6000, 'color': [0.259, 0.000, 0.400], 'number': 87}
elements[ 88] = elements[ 'Ra'] = {'symbol':  'Ra', 'name':        'radium', 'mass': 226.00000000, 'radius':  2.2100, 'color': [0.000, 0.490, 0.000], 'number': 88}
elements[ 89] = elements[ 'Ac'] = {'symbol':  'Ac', 'name':      'actinium', 'mass': 227.00000000, 'radius':  2.1500, 'color': [0.439, 0.671, 0.980], 'number': 89}
elements[ 90] = elements[ 'Th'] = {'symbol':  'Th', 'name':       'thorium', 'mass': 231.03588000, 'radius':  2.0600, 'color': [0.000, 0.729, 1.000], 'number': 90}
elements[ 91] = elements[ 'Pa'] = {'symbol':  'Pa', 'name':  'protactinium', 'mass': 232.03806000, 'radius':  2.0000, 'color': [0.000, 0.631, 1.000], 'number': 91}
elements[ 92] = elements[  'U'] = {'symbol':   'U', 'name':       'uranium', 'mass': 237.00000000, 'radius':  1.9600, 'color': [0.000, 0.561, 1.000], 'number': 92}
elements[ 93] = elements[ 'Np'] = {'symbol':  'Np', 'name':     'neptunium', 'mass': 238.02891000, 'radius':  1.9000, 'color': [0.000, 0.502, 1.000], 'number': 93}
elements[ 94] = elements[ 'Pu'] = {'symbol':  'Pu', 'name':     'plutonium', 'mass': 243.00000000, 'radius':  1.8700, 'color': [0.000, 0.420, 1.000], 'number': 94}
elements[ 95] = elements[ 'Am'] = {'symbol':  'Am', 'name':     'americium', 'mass': 244.00000000, 'radius':  1.8000, 'color': [0.329, 0.361, 0.949], 'number': 95}
elements[ 96] = elements[ 'Cm'] = {'symbol':  'Cm', 'name':        'curium', 'mass': 247.00000000, 'radius':  1.6900, 'color': [0.471, 0.361, 0.890], 'number': 96}
elements[ 97] = elements[ 'Bk'] = {'symbol':  'Bk', 'name':     'berkelium', 'mass': 247.00000000, 'radius':  1.6600, 'color': [0.541, 0.310, 0.890], 'number': 97}
elements[ 98] = elements[ 'Cf'] = {'symbol':  'Cf', 'name':   'californium', 'mass': 251.00000000, 'radius':  1.6800, 'color': [0.631, 0.212, 0.831], 'number': 98}
elements[ 99] = elements[ 'Es'] = {'symbol':  'Es', 'name':   'einsteinium', 'mass': 252.00000000, 'radius':  1.6500, 'color': [0.702, 0.122, 0.831], 'number': 99}
elements[100] = elements[ 'Fm'] = {'symbol':  'Fm', 'name':       'fermium', 'mass': 257.00000000, 'radius':  1.6700, 'color': [0.702, 0.122, 0.729], 'number': 100}
elements[101] = elements[ 'Md'] = {'symbol':  'Md', 'name':   'mendelevium', 'mass': 258.00000000, 'radius':  1.7300, 'color': [0.702, 0.051, 0.651], 'number': 101}
elements[102] = elements[ 'No'] = {'symbol':  'No', 'name':      'nobelium', 'mass': 259.00000000, 'radius':  1.7600, 'color': [0.741, 0.051, 0.529], 'number': 102}
elements[103] = elements[ 'Lr'] = {'symbol':  'Lr', 'name':    'lawrencium', 'mass': 266.00000000, 'radius':  1.6100, 'color': [0.780, 0.000, 0.400], 'number': 103}
elements[104] = elements[ 'Rf'] = {'symbol':  'Rf', 'name': 'rutherfordium', 'mass': 267.00000000, 'radius':  1.5700, 'color': [0.800, 0.000, 0.349], 'number': 104}
elements[105] = elements[ 'Db'] = {'symbol':  'Db', 'name':       'dubnium', 'mass': 268.00000000, 'radius':  1.4900, 'color': [0.820, 0.000, 0.310], 'number': 105}
elements[106] = elements[ 'Sg'] = {'symbol':  'Sg', 'name':    'seaborgium', 'mass': 269.00000000, 'radius':  1.4300, 'color': [0.851, 0.000, 0.271], 'number': 106}
elements[107] = elements[ 'Bh'] = {'symbol':  'Bh', 'name':       'bohrium', 'mass': 270.00000000, 'radius':  1.4100, 'color': [0.878, 0.000, 0.220], 'number': 107}
elements[108] = elements[ 'Hs'] = {'symbol':  'Hs', 'name':       'hassium', 'mass': 270.00000000, 'radius':  1.3400, 'color': [0.902, 0.000, 0.180], 'number': 108}
elements[109] = elements[ 'Mt'] = {'symbol':  'Mt', 'name':    'meitnerium', 'mass': 278.00000000, 'radius':  1.2900, 'color': [0.922, 0.000, 0.149], 'number': 109}
elements[110] = elements[ 'Ds'] = {'symbol':  'Ds', 'name':  'darmstadtium', 'mass': 281.00000000, 'radius':  1.2800, 'color': [0.922, 0.000, 0.149], 'number': 110}
elements[111] = elements[ 'Rg'] = {'symbol':  'Rg', 'name':   'roentgenium', 'mass': 282.00000000, 'radius':  1.2100, 'color': [0.922, 0.000, 0.149], 'number': 111}
elements[112] = elements[ 'Cn'] = {'symbol':  'Cn', 'name':   'copernicium', 'mass': 285.00000000, 'radius':  1.2200, 'color': [0.922, 0.000, 0.149], 'number': 112}
elements[113] = elements[ 'Nh'] = {'symbol':  'Nh', 'name':      'nihonium', 'mass': 286.00000000, 'radius':  1.3600, 'color': [0.922, 0.000, 0.149], 'number': 113}
elements[114] = elements[ 'Fl'] = {'symbol':  'Fl', 'name':     'flerovium', 'mass': 289.00000000, 'radius':  1.4300, 'color': [0.922, 0.000, 0.149], 'number': 114}
elements[115] = elements[ 'Mc'] = {'symbol':  'Mc', 'name':     'moscovium', 'mass': 290.00000000, 'radius':  1.5800, 'color': [0.922, 0.000, 0.149], 'number': 115}
elements[116] = elements[ 'Lv'] = {'symbol':  'Lv', 'name':   'livermorium', 'mass': 293.00000000, 'radius':  1.6600, 'color': [0.922, 0.000, 0.149], 'number': 116}
elements[117] = elements[ 'Ts'] = {'symbol':  'Ts', 'name':    'tennessine', 'mass': 294.00000000, 'radius':  1.5600, 'color': [0.922, 0.000, 0.149], 'number': 117}
elements[118] = elements[ 'Og'] = {'symbol':  'Og', 'name':     'oganesson', 'mass': 294.00000000, 'radius':  1.5700, 'color': [0.922, 0.000, 0.149], 'number': 118}
