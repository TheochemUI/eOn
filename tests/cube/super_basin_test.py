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
        return self.state.load(self.frameNewMinimum)

    def loadPass(self):
        return self.state.load(self.framePass)
    
    def rate(self, temperature):
        result=exp((self.energyNewMinimum - self.energyPass)/temperature)
        return result

class State:
    def __init__(self, number, temperature=None, path=''):
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
        return self.path+'%06d/' % self.state_number
    
    def get_process_table(self):
        result=list()
        for path in self.paths:
            p=dict()
            p['product']=path.numberNewMinimum
            p['rate']=path.rate(self.temperature)
            result.append(p)
        return result

    def load(self, frame=0):
        return loadcon(self.path+'%06d/%06d.con' % (self.state_number, frame))

class SuperBasinTestError:
    pass

def wrap(states):
    temperature=0.031066746727980595
    temperature=0.010857362047581296343
    #temperature=0.0072382413650541981
    #temperature=0.0054286810237906486
    #temperature=0.0001
    n=len(states)
    result=list()
    for i in states:
        result.append(State(i, temperature))
    return result

def isClose(d1, d2):
    if abs((d1-d2)/d2) > 1e-8:
        print 'Error'
        raise SuperBasinTestError()

def areClose(a1, a2):
    if isinstance(a1, float):
        isClose(a1, a2)
    else:
        if len(a1) != len(a2):
            print 'Error'
            raise SuperBasinTestError()
        for a, b in zip(a1, a2):
            areClose(a, b)

def test_super_basin():       
    fundfile='super_basin_test.txt'
    if os.path.exists(fundfile): os.remove(fundfile)
    
    # Test number one
    states=[0]
    states=wrap(states)
    sb=Superbasin(fundfile, states)
    isClose(sb.mean_residence_times[0], 49.5000495)
    isClose(sb.probability_matrix[0][0], 1.0)
    sb.delete()
    
    # Test number two
    states=[0, 1]
    states=wrap(states)
    sb=Superbasin(fundfile, states)
    mean_residence_times=array([4950.4950495, 4950.4950495])
    areClose(sb.mean_residence_times, mean_residence_times)
    probability_matrix=array([[0.50251231, 0.49748769], [0.49748769, 0.50251231]])
    areClose(sb.probability_matrix, probability_matrix)
    sb.delete()
    
    # Test number three
    states=[0, 1, 2]
    states=wrap(states)
    sb=Superbasin(fundfile, states)
    mean_residence_times= array([ 4999.25752388,  4999.75002487,    98.99762424])
    probability_matrix=array([[  5.04950255e-01,   4.99950252e-01,   4.94951245e-03],
           [  4.95000250e-03,   4.99999752e-03,   4.95000250e-05],
           [  4.90099743e-01,   4.95049750e-01,   9.95000988e-01]])
    areClose(sb.mean_residence_times, mean_residence_times)
    areClose(sb.probability_matrix, probability_matrix)
    sb.delete()
    
    # Test number four
    states=[0, 1, 5]
    states=wrap(states)
    sb=Superbasin(fundfile, states)
    mean_residence_times= array([4950.74128214, 4950.74376909, 49.99017362])
    probability_matrix=\
    array([[ 5.02512555e-01, 4.97487932e-01, 4.92513545e-05],
           [ 4.92562309e-01, 4.97537188e-01, 4.92562309e-05],
           [ 4.92513545e-03, 4.97487932e-03, 9.99901492e-01]])
    areClose(sb.mean_residence_times, mean_residence_times)
    areClose(sb.probability_matrix, probability_matrix)
    #sb.delete()
    n=len(states)
    p=[0]*n
    tries=10
    for i in range(tries):
        time, i=sb.pick_exit_state(states[2])
        p[i]+=1
    for i in range(n):
        p[i]=float(p[i])/tries
    print time, p

def run_as_client():
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