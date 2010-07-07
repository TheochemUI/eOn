'''
Con(figuration) i/o library
'''

import numpy

import atoms

def loadcon(filein):
    '''
    Load a con file
        filein - may be either a filename or a file-like object
    '''
    if hasattr(filein, 'readline'):
        con = filein
    else:
        con = open(filein, 'r')
    con.readline() #Line 1: comment
    con.readline() #Line 2: comment
    # determine how many dimensions
    tmp = numpy.array(con.readline().split()) #Line 3: Box lengths
    for i in range(len(tmp)):
        dim=i+1
        try: float(tmp[i])
        except:
            dim=i
            break
    #handle the box   
    BoxLength=numpy.zeros(dim)
    for i in range(dim):
        BoxLength[i]=float(tmp[i])
    con.readline()
    #BoxAngle=numpy.array([ float(f) for f in con.readline().split()[0:dim] ]) #Line 4: Box angles
    boxtemp=numpy.zeros((dim,dim),'d')
    for i in range(dim):
        # how do box angles work in con files??
        boxtemp[i][i]=BoxLength[i]
    con.readline() #Line 5: comment
    con.readline() #Line 6: comment
    num_types = int(con.readline().split()[0]) #Line 7: number of atom types
    num_each_type = con.readline().split() #line 8: number of each type of atom
    mass_of_type = con.readline().split() #line 9: mass of each type of atom
    num_atoms = 0
    for i in range(num_types):
        num_each_type[i] = int(num_each_type[i])
        mass_of_type[i] = float(mass_of_type[i])
        num_atoms += num_each_type[i]
    a = atoms.Atoms(num_atoms)
    a.box = boxtemp
    index = 0
    for i in range(num_types):
        name = con.readline().strip()
        con.readline() #skip meaningless line
        for j in range(num_each_type[i]):
            vals = con.readline().split()
            for k in range(dim):
                a.r[index][k] = float(vals[k])
            a.names[index] = name
            if not int(vals[dim])==0:
                a.free[index]=0       
            index += 1
    return a



def savecon(fileout, p, w = 'w'):
    '''
    Save a con file
        fileout - can be either a file name or a file-like object
        p    - information (in the form of an atoms object) to save
        w        - write/append flag
    '''
    if hasattr(fileout, 'write'):
        con = fileout
    else:
        con = open(fileout, w)
    print >> con, "Generated by eOn"
    print >> con
    dim = len(p.r[0])
    # Note: This will only work for orthogonal boxes
    BoxDiag=numpy.zeros((dim), 'd')
    for i in range(dim): 
        BoxDiag[i] = p.box[i][i]
    print >> con, " ".join(['%.5f' % s for s in BoxDiag])
    Angle = numpy.zeros((dim), 'd') + 90.0
    print >> con, " ".join(['%.5f' % s for s in Angle])
    print >> con
    print >> con
    atom_count = {}
    name_order = []
    for i in range(len(p)):
        name = p.names[i]
        if name not in name_order:
            name_order.append(name)
        if name in atom_count:
            atom_count[name] += 1
        else:
            atom_count[name] = 1
    print >> con, len(name_order)
    print >> con, " ".join([str(atom_count[i]) for i in name_order])
    print >> con, " ".join(["1.0" for i in name_order])
    index = 0
    for i in range(len(name_order)):
        print >> con, name_order[i]
        print >> con, "Coordinates of Component", i+1
        for j in range(atom_count[name_order[i]]):
            print >> con, p.r[index][0], p.r[index][1], p.r[index][2], int(not p.free[index]), index+1
            index += 1


def load_mode(modefilein):
    ''' 
    Reads a mode.dat file into an N by 3 numpy array
        modefilein - may be either a file-like object of a filename
    '''
    if hasattr(modefilein, 'readline'):
        f = modefilein
    else:
        f = open(modefilein, 'r')
    natoms = int(f.readline().split()[0])/3
    ret = numpy.zeros((natoms, 3))
    for i in range(natoms):
        line = f.readline().strip().split()
        for j in range(3):
            ret[i][j] = float(line[j])
    return ret

def save_mode(modefileout, displace_vector, reactant):
    '''
    Saves an Nx3 numpy array into a mode.dat file. 
        modefileout - may be either a filename or file-like object
        displace_vector - the mode (Nx3 numpy array)
        reactant - an atoms object that this mode is from. This is used to determine which atoms are free
    '''
    if hasattr(modefileout, 'write'):
        f = modefileout
    else:
        f = open(modefileout, 'w')
    f.write("%d %d\n" % (len(reactant) * 3, 3 * int(reactant.free.sum())))
    for i in range(len(displace_vector)):
        f.write("%.5f %.5f %.5f\n" % (displace_vector[i][0], 
            displace_vector[i][1], displace_vector[i][2]))
    if type(modefileout) is not file:
        f.close()


def save_results_dat(fileout, results):
    '''
    Saves a results.dat file from a dictionary
    '''
    if hasattr(fileout, 'write'):
        f = fileout
    else:
        f = open(fileout, 'w')
    
    for key in results:
        print >> f, results[key], key

def parse_results_dat(filein):
    '''
    Reads a results.dat file and gives a dictionary of the values contained
    therein
    '''
    if hasattr(filein, 'readline'):
        f = filein
    else:
        f = open(filein)
    results = {}
    for line in f:
        line = line.split()
        if len(line) < 2:
            continue
        if '.' in line[0]:
            results[line[1]] = float(line[0])
        else:
            results[line[1]] = int(line[0])
    f.close()
    #input validation
    if not type(results['termination_reason']) == int:
        raise TypeError('termination_reason')
    if not type(results['random_seed']) == int:
        raise TypeError('random_seed')
    if not type(results['reactant_state_tag']) == int:
        raise TypeError('reactant_state_tag')
    if not type(results['potential_tag']) == int:
        raise TypeError('potential_tag')
    if not type(results['force_calls']) == int:
        raise TypeError('force_calls')
    if not type(results['force_calls_saddle_point_concave']) == int:
        raise TypeError('force_calls_saddle_point_concave')
    if not type(results['force_calls_saddle_point_convex']) == int:
        raise TypeError('force_calls_saddle_point_convex')
    if not type(results['force_calls_prefactors']) == int:
        raise TypeError('force_calls_prefactors')
    if not type(results['prefactor_reactant_to_product']) == float:
        raise TypeError('prefactor_reactant_to_product')
    if not type(results['prefactor_product_to_reactant']) == float:
        raise TypeError('prefactor_product_to_reactant')
    if not type(results['potential_energy_saddle']) == float:
        raise TypeError('potential_energy_saddle')
    if not type(results['potential_energy_reactant']) == float:
        raise TypeError('potential_energy_reactant')
    if not type(results['potential_energy_product']) == float:
        raise TypeError('potential_energy_product')
    return results


def loadposcar(filename):
    '''
    Load the POSCAR file named filename
    Returns an atoms object
    '''
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    line_index = 0
    while line_index < len(lines):
        # Line 1: Atom types
        AtomTypes = lines[line_index].split() 
        # Line 2: scaling of coordinates
        line_index += 1
        scale = float(lines[line_index]) 
        # Lines 3-5: the box
        line_index += 1
        box = numpy.zeros((3, 3))
        for i in range(3):
            line = lines[line_index].split()
            line_index += 1
            box[i] = numpy.array([float(line[0]), float(line[1]), float(line[2])]) * scale
        # Line 6: number of atoms of each type. 
        line = lines[line_index].split()
        line_index += 1
        NumAtomsPerType = []
        for l in line:
            NumAtomsPerType.append(int(l))
        # Now have enough info to make the atoms object.
        num_atoms = sum(NumAtomsPerType)
        p = atoms.Atoms(num_atoms)
        # Fill in the box.
        p.box = box
        # Line 7: selective or cartesian
        sel = lines[line_index][0]
        line_index += 1
        selective_flag = (sel == 's' or sel == 'S')
        if not selective_flag: 
            car = sel
        else: 
            car = lines[line_index][0]
            line_index += 1
        direct_flag = not (car == 'c' or car == 'C' or car == 'k' or car == 'K')
        atom_index = 0
        for i in range(len(NumAtomsPerType)):
            for j in range(NumAtomsPerType[i]):
                p.names[atom_index] = AtomTypes[i]
                line = lines[line_index].split()
                line_index += 1
                if(selective_flag): 
                    assert len(line) >= 6
                else: 
                    assert len(line) >= 3
                pos = line[0:3]
                if selective_flag:
                    sel = line[3:7]
                    if sel[0] == 'T' or sel[0] == 't': 
                        p.free[atom_index] = 1
                    elif sel[0] == 'F' or sel[0] == 'f':
                        p.free[atom_index] = 0
                p.r[atom_index] = numpy.array([float(f) for f in pos])
                if direct_flag:
                    p.r[atom_index] = numpy.dot(p.r[atom_index], p.box)
                else: 
                    p.r[atom_index] *= scale
                atom_index += 1
    return p    


def saveposcar(filename, p, w='w', direct = False):
    '''
    Save a POSCAR
        filename - name to save it under
        point    - atoms object to save
        w        - write/append flag
    '''
    atom_types = []
    num_each_type = {}
    for name in p.names:
        if not name in atom_types:
            atom_types.append(name)
            num_each_type[name] = 1
        else:
            num_each_type[name] += 1
    poscar = open(filename, w)
    print >> poscar, " ".join(atom_types) #Line 1: Atom types
    print >> poscar, 1.0 #Line 2: scaling
    for i in range(3):
        print >> poscar, " ".join(['%20.14f' % s for s in p.box[i]]) #lines 3-5: box
    print >> poscar, " ".join(['%s' % num_each_type[key] for key in atom_types])
    print >> poscar, 'Selective Dynamics' #line 6: selective dynamics
    if direct:
        print >> poscar, 'Direct' #line 7 cartesian coordinates
        ibox = numpy.linalg.inv(numpy.array(p.box))
        p.r = numpy.dot(p.r, ibox)
    else:
        print >> poscar, 'Cartesian' #line 7 cartesian coordinates
    for i in range(len(p)):
            posline = " ".join(['%20.14f' % s for s in p.r[i]]) + " "
            for j in range(3):
                if(p.free[i]): 
                    posline+='   T'
                else: 
                    posline+='   F'
            print >> poscar, posline


if __name__=='__main__':
    import sys
    a = loadposcar(sys.argv[1])
    saveposcar("FOOCAR", a)
