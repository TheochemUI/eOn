#!/usr/bin/env python
import re, os, sys
sys.path.insert(0, '/Users/berthet/eon/')
from superbasin import *
from io import loadcon, savecon
from atoms import identical
from StringIO import StringIO

refloat='[\d\+\-\.e]+'

class Path:
    def __init__(self, file, state):
        """Constructor.
        A Path object is expected to be listed in a State obbject.
        file: open file from where to read the data.
        state: State object to which self is to belong.."""
        self.state=state
        self.energyPass=float(re.search('\|\s*energyPass_ (%s)' % refloat, file.readline()).group(1))
        self.framePass=int(re.search('framePass_ (\d+)', file.readline()).group(1))
        self.energyNewMinimum=float(re.search('energyNewMinimum_ (%s)' % refloat, file.readline()).group(1))
        self.frameNewMinimum=int(re.search('frameNewMinimum_ (\d+)', file.readline()).group(1))
        self.numberNewMinimum=int(re.search('numberNewMinimum_ (\d+)', file.readline()).group(1))
        file.readline()
    
    def __repr__(self):
        out=StringIO()
        print >>out, '|',
        print >>out, '\tenergyPass', self.energyPass
        print >>out, '\tframePass', self.framePass
        print >>out, '\tenergyNewMinimum', self.energyNewMinimum
        print >>out, '\tframeNewMinimum', self.frameNewMinimum
        print >>out, '\tnumberNewMinimum', self.numberNewMinimum
        return out.getvalue()
    
    def loadNewMinimum(self):
        """Load atomic configuration of the minimum the path leads to."""
        return self.state.load(self.frameNewMinimum)

    def loadPass(self):
        """Load atomic configuration at the saddle point."""
        return self.state.load(self.framePass)
    
    def rate(self, temperature):
        """Reaction rate of the process."""
        result=exp((self.energyNewMinimum - self.energyPass)/temperature)
        return result

class State:
    def __init__(self, number, temperature=None, path=''):
        """Container for a state defined as the minimum plus all the paths out of the minimum.
        number: state's number is used to determine the path to files from where to read the data.
        temperature: required to generate the rates.
        path: path to directory where the data are."""
        self.path=path
        self.temperature=temperature
        file = open(path+'%06d/transient_state.txt' % number, 'r')
        self.energy = float(re.search('energy_ (%s)' % refloat, file.readline()).group(1))
        self.number = int(re.search('stateNumber_ (\d+)', file.readline()).group(1))
        self.state_number=self.number
        n=int(re.search('paths_ (\d+)', file.readline()).group(1))
        self.paths=list()
        for i in range(0, n):
            self.paths.append(Path(file, self))
    
    def __repr__(self):
        out=StringIO()
        print >>out, 'energy', self.energy
        print >>out, 'number', self.number
        print >>out, 'nPath', len(self.paths)
        for path in self.paths:
            out.write(path)
        return out.getvalue()
    
    def get_directory(self):
        """Path to where the data are saved."""
        return self.path+'%06d/' % self.state_number
    
    def get_process_table(self):
        """List or processes.
        The function returns a list of dictionnaries. The list has one entry for per process (see self.paths). Each entry (i.e. dictionnary) has two entries: product: number of the state the process leads to; rate: rate of the process."""
        result=list()
        for path in self.paths:
            p=dict()
            p['product']=path.numberNewMinimum
            p['rate']=path.rate(self.temperature)
            result.append(p)
        return result

    def load(self, frame=0):
        """Load con file. In frame 0 (default) is the minimum, then odd frames are saddle points and even frames are the minima the saddle point in the preceding frame is leading to."""
        return loadcon(self.path+'%06d/%06d.con' % (self.state_number, frame))

class SuperBasinTestError:
    pass

def run_as_client():
    """Emulates the eOn client by reading the results from files."""
    a=loadcon('reactant_passed.con')
    path_states='/Users/berthet/eon/tests/cube/'
    for i in range(8):
        state=State(i, path=path_states)
        b=state.load()
        file=state.get_directory()+'super_basin_test.txt'
        if os.path.exists(file):
            index=int(open(file, 'r').read())
            index+=1
        else:
            index=0
        index%=6     
        open(file, 'w').write(str(index))
        path=state.paths[index]
        if identical(a, b, 0.1):
            savecon('reactant.con', b)
            b=path.loadPass()
            savecon('saddle.con', b)
            b=path.loadNewMinimum()
            savecon('product.con', b)
            file=open('results.dat', 'w')
            file.write("""0 termination_reason
            0 random_seed
            0 reactant_state_tag
            0 potential_tag
            0 force_calls
            0 force_calls_saddle_point_concave
            0 force_calls_saddle_point_convex
            %f potential_energy_saddle
            %f potential_energy_reactant
            %f potential_energy_product
            %f barrier_reactant_to_product
            %f barrier_product_to_reactant
            0.0 displacement_saddle_distance
            0 force_calls_prefactors
            1.0 prefactor_reactant_to_product
            1.0 prefactor_product_to_reactant""" % (path.energyPass, state.energy, \
                path.energyNewMinimum, path.energyPass - state.energy, path.energyPass - path.energyNewMinimum))
            file.close()
    
    open('mode.dat', 'w')

run_as_client()
