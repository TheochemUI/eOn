"""Atoms / structure helpers for the eOn Python server.

Storage and I/O
    :class:`eon.structure.Structure` (alias :class:`Atoms`) is the in-memory
    working type. On disk, structures are ``readcon.ConFrame`` values
    (see :func:`eon.structure.Structure.from_conframe`).

Geometry kernels
    PBC, neighbor lists (vesin), and process-atom selection live under
    :mod:`eon.geometry`. This module re-exports them for compatibility and
    keeps structure-matching / CNA helpers plus a radius/color overlay.
"""
import numpy
import logging
logger = logging.getLogger("atoms")

from eon.structure import Atoms, Structure  # noqa: F401


def _readcon_z_helpers():
    """readcon>=0.14.9 exports symbol_to_atomic_number / atomic_number_to_symbol."""
    import readcon

    to_z = getattr(readcon, "symbol_to_atomic_number", None) or getattr(
        readcon, "symbol_to_z", None
    )
    to_sym = getattr(readcon, "atomic_number_to_symbol", None) or getattr(
        readcon, "z_to_symbol", None
    )
    if to_z is None or to_sym is None:
        ver = getattr(readcon, "__version__", "unknown")
        raise ImportError(
            f"readcon {ver} has no symbol/Z helpers; need readcon>=0.14.9"
        )
    return to_z, to_sym


def atomic_number(symbol_or_z):
    """Z for a symbol. Integers pass through. Symbol lookup is readcon."""
    if isinstance(symbol_or_z, (int, numpy.integer)):
        return int(symbol_or_z)
    if symbol_or_z in {"Xx", "X"}:
        return 0
    to_z, _ = _readcon_z_helpers()
    z = int(to_z(symbol_or_z))
    if z > 0:
        return z
    raise KeyError(f"unknown element {symbol_or_z!r}")


def symbol_for_z(z):
    """Chemical symbol for Z. Lookup is readcon; Z=0 is the unknown placeholder."""
    z = int(z)
    if z == 0:
        return "Xx"
    _, to_sym = _readcon_z_helpers()
    symbol = to_sym(z)
    if symbol and symbol not in {"X", "Xx"}:
        return str(symbol)
    raise KeyError(f"unknown Z {z!r}")


from eon.geometry import (  # noqa: F401
    box_to_length_angle,
    brute_neighbor_list,
    coordination_numbers,
    get_process_atoms,
    least_coordinated,
    length_angle_to_box,
    neighbor_list,
    neighbor_list_pairs,
    neighbor_list_vectors,
    pbc,
    per_atom_norm,
    per_atom_norm_gen,
)

# --- structure comparison / CNA (unchanged algorithms) ---

def identical(atoms1, atoms2, epsilon_r):
    """True if two structures match when same-element atoms are interchangeable.

    Parameters
    ----------
    atoms1, atoms2 : Structure
        Configurations to compare (same box within 1e-4).
    epsilon_r : float
        Max allowed MIC displacement (Å) for a pair to count as the same site.
    """
    if len(atoms1) != len(atoms2):
        return False

    for i in range(3):
        for j in range(3):
            if abs(atoms1.box[i][j] - atoms2.box[i][j]) > 0.0001:
                logger.warning(
                    "Identical returned false because boxes were not the same"
                )
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


def match(a, b, eps_r, neighbor_cutoff, indistinguishable,
          check_rotation=False, use_identical=False):
    if len(a) != len(b):
        return False

    if check_rotation:
        if indistinguishable and use_identical:
            return get_mappings(a, b, eps_r, neighbor_cutoff)
        else:
            return rot_match(a, b, eps_r)
    else:
        if indistinguishable and use_identical:
            return identical(a, b, eps_r)
        else:
            diff = pbc(a.r - b.r, a.box)
            return numpy.max(numpy.sum(diff**2.0, axis=1)) < eps_r**2.0
            #return max(per_atom_norm(a.r-b.r, a.box))<eps_r


def point_energy_match(file_a, energy_a, file_b, energy_b, eps_e, eps_r,
                       neighbor_cutoff, check_rotation=False, use_identical=False):
    import eon.fileio as io
    if abs(energy_a - energy_b) > eps_e:
        return False
    a = io.loadcon(file_a)
    b = io.loadcon(file_b)
    if match(a, b, eps_r, neighbor_cutoff, False,
             check_rotation=check_rotation, use_identical=use_identical):
        return True
    return False


def points_energies_match(file_a, energy_a, files_b, energies_b, eps_e, eps_r,
                          neighbor_cutoff, check_rotation=False, use_identical=False):
    for i in range(len(files_b)):
        if point_energy_match(file_a, energy_a, files_b[i], energies_b[i],
                              eps_e, eps_r, neighbor_cutoff,
                              check_rotation=check_rotation,
                              use_identical=use_identical):
            return i
    return None


def rot_match(a, b, eps_r):
    if not (a.free.all() and b.free.all()):
        logger.warning("Comparing structures with frozen atoms with rotational matching; check_rotation may be set incorrectly")
    if len(a) == 0:
        return len(b) == 0
    try:
        from pyeonclient import _core

        ira = getattr(_core, "ira_match", None)
        if ira is not None:
            z1 = numpy.asarray([atomic_number(s) for s in a.names], dtype=numpy.int64)
            z2 = numpy.asarray([atomic_number(s) for s in b.names], dtype=numpy.int64)
            hd, err = ira(
                numpy.ascontiguousarray(a.r, dtype=float),
                z1,
                numpy.ascontiguousarray(b.r, dtype=float),
                z2,
                float(eps_r),
            )
            if err == 0:
                return hd < eps_r
    except Exception:
        pass
    return _rot_match_kabsch(a, b, eps_r)


def _rot_match_kabsch(a, b, eps_r):
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
    return dist < eps_r


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
    ct = numpy.cos(theta)
    st = numpy.sin(theta)
    mag = numpy.linalg.norm(axis)
    if (mag*mag == 0 or theta == 0.0):
        return numpy.identity(3)
    return numpy.array([
        [u2 +(v2 +w2)*ct, u*v*(1-ct)-w*mag*st, u*w*(1-ct)+v*mag*st],
        [u*v*(1-ct)+w*mag*st, v2 +(u2 +w2)*ct, v*w*(1-ct)-u*mag*st],
        [u*w*(1-ct)-v*mag*st, v*w*(1-ct)+u*mag*st, w2 +(v2 +u2)*ct]
        ])/(mag*mag)


# Same labels as eonc::EpiCenters::cna (client/EpiCenters.cpp).
CNA_FCC = 0
CNA_HCP = 1
CNA_OTHER = 2


def _adjacency_from_pairs(n, i, j):
    """Unique-index adjacency from vesin pair arrays.

    Each row is increasing atom index, matching ``EpiCenters.cpp`` insertion
    order so the CNA bond-sum still distinguishes 421 (fcc) from 422 (hcp).
    """
    i = numpy.asarray(i, dtype=numpy.int64)
    j = numpy.asarray(j, dtype=numpy.int64)
    mask = i != j
    i, j = i[mask], j[mask]
    nl = [[] for _ in range(n)]
    if i.size == 0:
        return nl
    packed = numpy.unique(i * numpy.int64(n) + j)
    for a, b in zip((packed // n).tolist(), (packed % n).tolist()):
        nl[int(a)].append(int(b))
    return nl


def _neighbor_lists(p, cutoff):
    """CNA adjacency from :func:`eon.geometry.neighbor_list_pairs`."""
    i, j, _s = neighbor_list_pairs(p, cutoff)
    return _adjacency_from_pairs(len(p), i, j)


def _neighbor_sets(nl):
    """Per-atom neighbor sets for O(1) membership; list order stays on *nl*."""
    return [set(nbs) for nbs in nl]


def _common_neighbors_ordered(nl_a1, nbs_a2):
    """Common neighbors of a1 and a2, in a1 neighbor-list order.

    Bond-sum (FCC 421 vs HCP 422) depends on that order, so this must not
    be an unordered set intersection.
    """
    return [a3 for a3 in nl_a1 if a3 in nbs_a2]


def _common_bond_stats(nl_sets, common):
    bonds_nr = 0
    bonds_sum = 0
    for j2 in range(1, len(common)):
        nbs = nl_sets[common[j2]]
        for j1 in range(j2):
            if common[j1] in nbs:
                bonds_nr += 1
                bonds_sum += j1 + j2
    return bonds_nr, bonds_sum


def cna(p, cutoff, brute=False):
    """Common-neighbor labels for every atom in *p*.

    Labels match the C++ client (``EpiCenters::cna``): 0 fcc (421), 1 hcp
    (422), 2 other. Inspired by the CNA code provided by Asap
    (wiki.fysik.dtu.dk/asap).
    """
    n = len(p)
    can_values = numpy.full(n, CNA_OTHER, dtype=int)
    nr_FCC = numpy.zeros(n, dtype=int)
    nr_HCP = numpy.zeros(n, dtype=int)
    nl = _neighbor_lists(p, cutoff)
    nl_sets = _neighbor_sets(nl)

    for a2 in range(n):
        nbs_a2 = nl_sets[a2]
        for a1 in nl[a2]:
            if a1 >= a2:
                continue
            common = _common_neighbors_ordered(nl[a1], nbs_a2)
            if len(common) != 4:
                continue
            bonds_nr, bonds_sum = _common_bond_stats(nl_sets, common)
            if bonds_nr == 2:
                if bonds_sum == 6:
                    nr_FCC[a1] += 1
                    nr_FCC[a2] += 1
                else:
                    nr_HCP[a1] += 1
                    nr_HCP[a2] += 1

    for i in range(n):
        if len(nl[i]) == 12:
            if nr_FCC[i] == 12:
                can_values[i] = CNA_FCC
            elif nr_FCC[i] == 6 and nr_HCP[i] == 6:
                can_values[i] = CNA_HCP
    return can_values

def not_HCP_or_FCC(p, cutoff, brute=False):
    """Indices of atoms that are neither fcc nor hcp (CNA other)."""
    cna_numbers = cna(p, cutoff, brute)
    return [i for i, label in enumerate(cna_numbers) if label == CNA_OTHER]


def cnat(p, cutoff, brute=False):
    """ Returns a list of cna numbers for all atoms in p
        Inspired by the CNA code provided by Asap (wiki.fysik.dtu.dk/asap)"""
    can_values = numpy.zeros(len(p))
    nr_5 = numpy.zeros(len(p))
    nr_6 = numpy.zeros(len(p))
    nl = _neighbor_lists(p, cutoff)
    nl_sets = _neighbor_sets(nl)

    # loops over all the atoms
    for a2 in range(len(p)):
        nbs_a2 = nl_sets[a2]
        for a1 in nl[a2]:
            if a1 < a2:
                common = _common_neighbors_ordered(nl[a1], nbs_a2)
                # determines the connectivity of common neighbors
                if len(common) in [5,6]:
                    bonds_nr, bonds_sum = _common_bond_stats(nl_sets, common)

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

    nl = _neighbor_lists(p, cutoff)
    nl_sets = _neighbor_sets(nl)

    def codeString(j,k,l):
        return "%d,%d,%d" % (j,k,l)

    # loops over all the atoms
    for a2 in range(len(p)):
        nbs_a2 = nl_sets[a2]
        for a1 in nl[a2]:
            if a1 < a2: # prevent double counting?
                common = _common_neighbors_ordered(nl[a1], nbs_a2)
                bonds_nr, bonds_sum = _common_bond_stats(nl_sets, common)
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
def get_mappings(a, b, eps_r, neighbor_cutoff, mappings=None):
    """Depth-first search for a complete atom mapping from a onto b.

    Returns None if no mapping was found, or a dict mapping a indices to b
    indices. Mirror images still map; this does not test proper rotation.
    """
    if mappings is None:
        b_coord = coordination_numbers(b, neighbor_cutoff)
        counts = {}
        for c in b_coord:
            counts[c] = counts.get(c, 0) + 1
        least = min(counts, key=counts.get)
        a_coord = coordination_numbers(a, neighbor_cutoff)
        try:
            a_atom = list(a_coord).index(least)
        except ValueError:
            return None
        for i, c in enumerate(b_coord):
            if c != least or a.names[a_atom] != b.names[i]:
                continue
            found = get_mappings(a, b, eps_r, neighbor_cutoff, {a_atom: i})
            if found is not None:
                return found
        return None

    n = len(a)
    mapped_a = mappings
    mapped_b = set(mapped_a.values())
    unmapped_a = next((i for i in range(n) if i not in mapped_a), n)
    mapped_idx = list(mapped_a.keys())
    diffs = pbc(a.r[unmapped_a] - a.r[mapped_idx], a.box)
    dists = numpy.linalg.norm(numpy.atleast_2d(diffs), axis=1)
    b_of_mapped = [mapped_a[i] for i in mapped_idx]
    want_name = a.names[unmapped_a]
    for b_atom in range(len(b)):
        if b_atom in mapped_b or b.names[b_atom] != want_name:
            continue
        b_diffs = pbc(b.r[b_atom] - b.r[b_of_mapped], b.box)
        b_dists = numpy.linalg.norm(numpy.atleast_2d(b_diffs), axis=1)
        if numpy.max(numpy.abs(dists - b_dists)) > eps_r:
            continue
        new_map = mapped_a.copy()
        new_map[unmapped_a] = b_atom
        if len(new_map) == n:
            return new_map
        found = get_mappings(a, b, eps_r, neighbor_cutoff, new_map)
        if found is not None:
            return found
    return None

def get_rotation_matrix(axis, theta):
    axis = axis / numpy.linalg.norm(axis)
    t = theta
    ct = numpy.cos(t)
    st = numpy.sin(t)
    T = 1.0 - ct
    rx, ry, rz = axis
    rotmat = numpy.zeros((3, 3))
    rotmat[0][0] = T*rx*rx + ct
    rotmat[0][1] = T*ry*rx + rz*st
    rotmat[0][2] = T*rz*rx - ry*st
    rotmat[1][0] = T*rx*ry - rz*st
    rotmat[1][1] = T*ry*ry + ct
    rotmat[1][2] = T*rz*ry + rx*st
    rotmat[2][0] = T*rx*rz + ry*st
    rotmat[2][1] = T*ry*rz - rx*st
    rotmat[2][2] = T*rz*rz + ct
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
    theta1 = numpy.arccos((a0a1*b0b1).sum())
    b.r = rotate(b.r, axis1, a.r[0], theta1)
    axis2 = (a.r[2] - a.r[0]) / numpy.linalg.norm(a.r[2] - a.r[0])
    va = a.r[2] - ((a.r[2] - a.r[0]) * axis2).sum() * axis2
    va = va / numpy.linalg.norm(va)
    vb = b.r[2] - ((b.r[2] - a.r[0]) * axis2).sum() * axis2
    vb = vb / numpy.linalg.norm(vb)
    theta2 = numpy.arccos((va * vb).sum())
    b.r = rotate(b.r, axis2, a.r[0], theta2)
    return b


# Presentation overlay keyed by Z. Symbol / name / mass live in readcon.
_RADIUS = (
    1.0000,  # 0
    0.3100,  # 1
    0.2800,  # 2
    1.2800,  # 3
    0.9600,  # 4
    0.8400,  # 5
    0.7300,  # 6
    0.7100,  # 7
    0.6600,  # 8
    0.5700,  # 9
    0.5800,  # 10
    1.6600,  # 11
    1.4100,  # 12
    1.2100,  # 13
    1.1100,  # 14
    1.0700,  # 15
    1.0500,  # 16
    1.0200,  # 17
    1.0600,  # 18
    2.0300,  # 19
    1.7600,  # 20
    1.7000,  # 21
    1.6000,  # 22
    1.5300,  # 23
    1.3900,  # 24
    1.3900,  # 25
    1.3200,  # 26
    1.2600,  # 27
    1.2400,  # 28
    1.3200,  # 29
    1.2200,  # 30
    1.2200,  # 31
    1.2000,  # 32
    1.1900,  # 33
    1.2000,  # 34
    1.2000,  # 35
    1.1600,  # 36
    2.2000,  # 37
    1.9500,  # 38
    1.9000,  # 39
    1.7500,  # 40
    1.6400,  # 41
    1.5400,  # 42
    1.4700,  # 43
    1.4600,  # 44
    1.4200,  # 45
    1.3900,  # 46
    1.4500,  # 47
    1.4400,  # 48
    1.4200,  # 49
    1.3900,  # 50
    1.3900,  # 51
    1.3800,  # 52
    1.3900,  # 53
    1.4000,  # 54
    2.4400,  # 55
    2.1500,  # 56
    2.0700,  # 57
    2.0400,  # 58
    2.0300,  # 59
    2.0100,  # 60
    1.9900,  # 61
    1.9800,  # 62
    1.9800,  # 63
    1.9600,  # 64
    1.9400,  # 65
    1.9200,  # 66
    1.9200,  # 67
    1.8900,  # 68
    1.9000,  # 69
    1.8700,  # 70
    1.8700,  # 71
    1.7500,  # 72
    1.7000,  # 73
    1.6200,  # 74
    1.5100,  # 75
    1.4400,  # 76
    1.4100,  # 77
    1.3600,  # 78
    1.3600,  # 79
    1.3200,  # 80
    1.4500,  # 81
    1.4600,  # 82
    1.4800,  # 83
    1.4000,  # 84
    1.5000,  # 85
    1.5000,  # 86
    2.6000,  # 87
    2.2100,  # 88
    2.1500,  # 89
    2.0600,  # 90
    2.0000,  # 91
    1.9600,  # 92
    1.9000,  # 93
    1.8700,  # 94
    1.8000,  # 95
    1.6900,  # 96
    1.6600,  # 97
    1.6800,  # 98
    1.6500,  # 99
    1.6700,  # 100
    1.7300,  # 101
    1.7600,  # 102
    1.6100,  # 103
    1.5700,  # 104
    1.4900,  # 105
    1.4300,  # 106
    1.4100,  # 107
    1.3400,  # 108
    1.2900,  # 109
    1.2800,  # 110
    1.2100,  # 111
    1.2200,  # 112
    1.3600,  # 113
    1.4300,  # 114
    1.5800,  # 115
    1.6600,  # 116
    1.5600,  # 117
    1.5700,  # 118
)
_COLOR = (
    (1.000, 0.078, 0.576),  # 0
    (1.000, 1.000, 1.000),  # 1
    (0.851, 1.000, 1.000),  # 2
    (0.800, 0.502, 1.000),  # 3
    (0.761, 1.000, 0.000),  # 4
    (1.000, 0.710, 0.710),  # 5
    (0.565, 0.565, 0.565),  # 6
    (0.188, 0.314, 0.973),  # 7
    (1.000, 0.051, 0.051),  # 8
    (0.565, 0.878, 0.314),  # 9
    (0.702, 0.890, 0.961),  # 10
    (0.671, 0.361, 0.949),  # 11
    (0.541, 1.000, 0.000),  # 12
    (0.749, 0.651, 0.651),  # 13
    (0.941, 0.784, 0.627),  # 14
    (1.000, 0.502, 0.000),  # 15
    (1.000, 1.000, 0.188),  # 16
    (0.122, 0.941, 0.122),  # 17
    (0.502, 0.820, 0.890),  # 18
    (0.561, 0.251, 0.831),  # 19
    (0.239, 1.000, 0.000),  # 20
    (0.902, 0.902, 0.902),  # 21
    (0.749, 0.761, 0.780),  # 22
    (0.651, 0.651, 0.671),  # 23
    (0.541, 0.600, 0.780),  # 24
    (0.611, 0.478, 0.780),  # 25
    (0.878, 0.400, 0.200),  # 26
    (0.941, 0.565, 0.627),  # 27
    (0.314, 0.816, 0.314),  # 28
    (0.784, 0.502, 0.200),  # 29
    (0.490, 0.502, 0.690),  # 30
    (0.761, 0.561, 0.561),  # 31
    (0.400, 0.561, 0.561),  # 32
    (0.741, 0.502, 0.890),  # 33
    (1.000, 0.631, 0.000),  # 34
    (0.651, 0.161, 0.161),  # 35
    (0.361, 0.722, 0.820),  # 36
    (0.439, 0.180, 0.690),  # 37
    (0.000, 1.000, 0.000),  # 38
    (0.580, 1.000, 1.000),  # 39
    (0.580, 0.878, 0.878),  # 40
    (0.451, 0.761, 0.788),  # 41
    (0.329, 0.710, 0.710),  # 42
    (0.231, 0.620, 0.620),  # 43
    (0.141, 0.561, 0.561),  # 44
    (0.039, 0.490, 0.549),  # 45
    (0.000, 0.412, 0.522),  # 46
    (0.753, 0.753, 0.753),  # 47
    (1.000, 0.851, 0.561),  # 48
    (0.651, 0.459, 0.451),  # 49
    (0.400, 0.502, 0.502),  # 50
    (0.620, 0.388, 0.710),  # 51
    (0.831, 0.478, 0.000),  # 52
    (0.580, 0.000, 0.580),  # 53
    (0.259, 0.620, 0.690),  # 54
    (0.341, 0.090, 0.561),  # 55
    (0.000, 0.788, 0.000),  # 56
    (0.439, 0.831, 1.000),  # 57
    (1.000, 1.000, 0.780),  # 58
    (0.851, 1.000, 0.780),  # 59
    (0.780, 1.000, 0.780),  # 60
    (0.639, 1.000, 0.780),  # 61
    (0.561, 1.000, 0.780),  # 62
    (0.380, 1.000, 0.780),  # 63
    (0.271, 1.000, 0.780),  # 64
    (0.189, 1.000, 0.780),  # 65
    (0.122, 1.000, 0.780),  # 66
    (0.000, 1.000, 0.612),  # 67
    (0.000, 0.902, 0.459),  # 68
    (0.000, 0.831, 0.322),  # 69
    (0.000, 0.749, 0.220),  # 70
    (0.000, 0.671, 0.141),  # 71
    (0.302, 0.761, 1.000),  # 72
    (0.302, 0.651, 1.000),  # 73
    (0.129, 0.580, 0.839),  # 74
    (0.149, 0.490, 0.671),  # 75
    (0.149, 0.400, 0.588),  # 76
    (0.090, 0.329, 0.529),  # 77
    (0.816, 0.816, 0.878),  # 78
    (1.000, 0.820, 0.137),  # 79
    (0.722, 0.722, 0.816),  # 80
    (0.651, 0.329, 0.302),  # 81
    (0.341, 0.349, 0.380),  # 82
    (0.620, 0.310, 0.710),  # 83
    (0.671, 0.361, 0.000),  # 84
    (0.459, 0.310, 0.271),  # 85
    (0.259, 0.510, 0.588),  # 86
    (0.259, 0.000, 0.400),  # 87
    (0.000, 0.490, 0.000),  # 88
    (0.439, 0.671, 0.980),  # 89
    (0.000, 0.729, 1.000),  # 90
    (0.000, 0.631, 1.000),  # 91
    (0.000, 0.561, 1.000),  # 92
    (0.000, 0.502, 1.000),  # 93
    (0.000, 0.420, 1.000),  # 94
    (0.329, 0.361, 0.949),  # 95
    (0.471, 0.361, 0.890),  # 96
    (0.541, 0.310, 0.890),  # 97
    (0.631, 0.212, 0.831),  # 98
    (0.702, 0.122, 0.831),  # 99
    (0.702, 0.122, 0.729),  # 100
    (0.702, 0.051, 0.651),  # 101
    (0.741, 0.051, 0.529),  # 102
    (0.780, 0.000, 0.400),  # 103
    (0.800, 0.000, 0.349),  # 104
    (0.820, 0.000, 0.310),  # 105
    (0.851, 0.000, 0.271),  # 106
    (0.878, 0.000, 0.220),  # 107
    (0.902, 0.000, 0.180),  # 108
    (0.922, 0.000, 0.149),  # 109
    (0.922, 0.000, 0.149),  # 110
    (0.922, 0.000, 0.149),  # 111
    (0.922, 0.000, 0.149),  # 112
    (0.922, 0.000, 0.149),  # 113
    (0.922, 0.000, 0.149),  # 114
    (0.922, 0.000, 0.149),  # 115
    (0.922, 0.000, 0.149),  # 116
    (0.922, 0.000, 0.149),  # 117
    (0.922, 0.000, 0.149),  # 118
)
numElements = len(_RADIUS)


class _ElementStyle:
    """Radius/color by Z; ``elements[symbol]`` resolves Z through readcon."""

    def __getitem__(self, key):
        z = key if isinstance(key, int) else atomic_number(key)
        if z < 0 or z >= len(_RADIUS):
            z = 0
        return {
            "radius": _RADIUS[z],
            "color": list(_COLOR[z]),
            "number": z,
            "symbol": symbol_for_z(z),
        }

    def __contains__(self, key):
        try:
            self[key]
            return True
        except Exception:
            return False


elements = _ElementStyle()
