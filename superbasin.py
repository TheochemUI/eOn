import os
import math
import numpy
from numpy import array #XXX: For repr-eval
import logging
logger = logging.getLogger('superbasin')


class Superbasin:
 
    def __init__(self, path, kT, state_list = None, get_state = None):
        #FIXME: self.states is literally a list of states, while in the superbasinscheme
        # self.states is a StateList object. Some renaming should happen.
        if state_list is None and get_state is None:
            raise ValueError('Superbasin must either have a list of states or a reference to get_state of a StateList')
        self.kT = kT
        self.path = path
        if not os.path.isfile(self.path):
            self.states = state_list
            self.state_numbers = [i.number for i in self.states]
            self._calculate_stuff()
            self.write_data()
        else:
            self.read_data(get_state)
      
    def pick_exit_state(self, entry_state):        
        #Find the exit state
        
        for i in range(len(self.state_numbers)):
            if entry_state.number == self.state_numbers[i]:
                entry_state_index = i
                break
        else:
            raise ValueError('Passed entry state is not in this superbasin')

        probability_vector = self.probability_matrix.transpose()[entry_state_index]
        if abs(1.0-numpy.sum(probability_vector)) < 1e-3:
            logger.warning("the probability vector isn't close to 1.0")
        probability_vector /= numpy.sum(probability_vector)
        
        u = numpy.random.random_sample()
        p = 0.0
        for i in range(len(self.states)):
            p += probability_vector[i]
            if p>u:
                exit_state_index = i 
                break
        else:
            logger.warning("Warning: failed to select exit state. p = " + str(p))
        time = self.mean_residence_times[entry_state_index]
        
         
        return time, exit_state_index

    def step(self, entry_state, get_product_state):
        time, exit_state_index=self.pick_exit_state(entry_state)
        exit_state=self.states[exit_state_index]
        #make a rate table for the exit state
        rate_table = []
        ratesum = 0.0
        process_table=exit_state.get_process_table()
        for process in process_table.values():
            if process['product'] not in self.state_numbers:
                rate_table.append(process['rate'])
                ratesum += process['rate']
       
        p = 0.0
        u = numpy.random.random_sample()
        for i in range(len(rate_table)):
            p += rate_table[i]/ratesum
            if p>u:
                exit_process_index = i 
                break
        else:
            logger.warning("Warning: failed to select rate. p = " + str(p))
        
        #akmc expects a state, not a process id
        exit_state = get_product_state(entry_state.number, exit_process_index)
        return time, exit_state

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self): 
        #Check detailed balance
        eq_prob = []
        for i in range(len(self.states)-1):
            for j in range(i+1, len(self.states)):
                proc_tab_i = self.states[i].get_process_table()
                proc_tab_j = self.states[j].get_process_table()
                
                rate_sum = 0.0
                pij = 0.0
                for id in proc_tab_i:
                    rate_sum += proc_tab_i[id]['rate']
                    if proc_tab_i[id]['product'] == self.state_numbers[j]:
                        pij += proc_tab_i[id]['rate']
                pij/=rate_sum

                rate_sum = 0.0
                pji = 0.0
                for id in proc_tab_j:
                    rate_sum += proc_tab_j[id]['rate']
                    if proc_tab_j[id]['product'] == self.state_numbers[i]:
                        pji += proc_tab_j[id]['rate']
                pji/=rate_sum

                #XXX: Hard coded comparison criteria
                if abs(math.exp((self.states[j].get_energy() - self.states[i].get_energy())/self.kT)*pij - pji) > .1:
                    logger.warning('states %d and %d do not satisfy detailed balance' % (self.state_numbers[i], self.state_numbers[j]))


        
        
        recurrent_vector = numpy.zeros(len(self.states))
        transient_matrix= numpy.zeros((len(self.states), len(self.states)))
        print 'remove'
        sum=0.0
        for i in range(len(self.states)):
            proc_table = self.states[i].get_process_table()
            for process in proc_table.values():
                sum+=process['rate']
                if process['product']==-1 or process['product'] not in self.state_numbers: 
                    recurrent_vector[i] += process['rate']
                else:
                    j = self.state_numbers.index(process['product'])
                    transient_matrix[i][j] += process['rate']
                transient_matrix[i][i] -= process['rate']
        print 'sum', sum
        print 'transient_matrix\n', transient_matrix
        print 'recurrent_vector', recurrent_vector

        
        #What is this stuff?
        #inverse free method
        #entry_state_vector = numpy.zeros(len(self.states))
        #entry_state_vector[self.state_numbers.index(entry_state.number)] = 1.0
        #fundamental_vector = numpy.linalg.solve(transient_matrix, entry_state_vector)
        # probabilitity to be the exit state
        #self.probability_vector = fundamental_vector/recurrent_vector
        #Calculate mean residence time
        #self.mean_residence_time = numpy.sum(fundamental_vector)
        
        fundamental_matrix = numpy.linalg.inv(transient_matrix)
        print 'fundamental_matrix\n', fundamental_matrix
        n=len(fundamental_matrix)
        id=numpy.zeros((n, n))
        self.mean_residence_times = numpy.zeros(len(self.states))
        self.probability_matrix = numpy.zeros((len(self.states), len(self.states)))
        
        for i in range(len(self.states)):
            for j in range(len(self.states)):
                self.mean_residence_times[j] -= fundamental_matrix[i][j]
                self.probability_matrix[i][j] = -recurrent_vector[i]*fundamental_matrix[i][j]
        print 'mean_residence_times=', repr(self.mean_residence_times)
        print 'probability_matrix=\\\n', repr(self.probability_matrix)
        print 'probability_matrix.sum', self.probability_matrix.sum(0)

    def write_data(self):
        logger.debug('saving data to %s' %self.path)
        f = open(self.path, 'w')
        print >> f, repr(self.state_numbers)
        print >> f, repr(self.mean_residence_times)
        print >> f, repr(self.probability_matrix)
        f.close()

    def read_data(self, get_state):
        logger.debug('reading data from %s' % self.path)
        print 'path', self.path
        f = open(self.path, 'r')
        self.state_numbers = eval(f.readline())
        print 'state_numbers', self.state_numbers
        self.mean_residence_times=eval(f.readline())
        print 'mean_residence_times', self.mean_residence_times
        self.probability_matrix=eval(f.read())
        print 'probability_matrix', self.probability_matrix

    def delete(self):
        logger.debug('deleting %s' % self.path)
        os.remove(self.path)
        self.states = None
        self.probability_matrix = None
        self.mean_residence_times = None

          



            
