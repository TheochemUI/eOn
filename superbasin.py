import os
import numpy
from numpy import array #XXX: For repr-eval
import logging
logger = logging.getLogger('superbasin')


class Superbasin:
    
    def __init__(self, path, state_list = None, get_state = None):
        #FIXME: self.states is literally a list of states, while in the superbasinscheme
        # self.states is a StateList object. Some renaming should happen.
        if state_list is None and get_state is None:
            raise ValueError('Superbasin must either have a list of states or a reference to get_state of a StateList')
        
        self.path = path
        if not os.path.isfile(self.path):
            self.states = state_list
            self.state_numbers = [i.number for i in self.states]
            self._calculate_stuff()
            self.write_data()
        else:
            self.read_data(get_state)

        
    
    def step(self, entry_state, get_product_state):        
        #Find the exit state
        
        for i in range(len(self.state_numbers)):
            if entry_state.number == self.state_numbers[i]:
                entry_state_index = i
                break
        else:
            raise ValueError('Passed entry state is not in this superbasin')

        probability_vector = self.probability_matrix[entry_state_index]
        probability_vector /= numpy.sum(probability_vector)
        if abs(1.0-numpy.sum(probability_vector)) < 1e-3:
            logger.warning("the probability vector isn't close to 1.0")
        
        u = numpy.random.random_sample()
        p = 0.0
        for i in range(len(self.states)):
            p += probability_vector[i]
            if p>u:
                exit_state_index = i 
                break
        else:
            logger.warning("Warning: failed to select exit state. p = " + str(p))
        
        #XXX: Do i even need to be keeping track of the exit_state_index?
        exit_state = self.states[exit_state_index]

        #make a rate table for the exit state
        rate_table = []
        ratesum = 0.0
        for key in exit_state.get_process_table(): 
            process = exit_state.get_process_table()[key]
            if process['product']==-1 or process['product'] not in self.state_numbers:
                rate_table.append([key, process['rate']])
                ratesum += process['rate']
        
        u = numpy.random.random_sample()
        for i in range(len(rate_table)):
            p += rate_table[i][1]/ratesum
            if p>u:
                nsid = i 
                break
        else:
            logger.warning("Warning: failed to select rate. p = " + str(p))

        absorbing_state = get_product_state(exit_state.number, rate_table[nsid][0])
        time = self.mean_residence_times[entry_state_index]

        if absorbing_state in self.states:
            #if we discovered a new connection between states in this basin
            #we need to update everything

            #TODO: need to ensure detailed balance for the subspace
            self._calculate_stuff()
            self.write_data()
        return time, absorbing_state

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self):
       
        recurrent_vector = numpy.zeros(len(self.states))
        transient_matrix= numpy.zeros((len(self.states), len(self.states)))
        for i in range(len(self.states)):
            proc_table = self.states[i].get_process_table()
            for key in proc_table:
                process = proc_table[key]

                if process['product']==-1 or process['product'] not in self.state_numbers: 
                    recurrent_vector[i] += process['rate']
                #elif process['product'] == self.state_numbers[i]:
                else:
                    j = self.state_numbers.index(process['product'])
                    transient_matrix[i][j] += process['rate']
                transient_matrix[i][i] -= process['rate']
        
        #What is this stuff?
        #inverse free method
        #entry_state_vector = numpy.zeros(len(self.states))
        #entry_state_vector[self.state_numbers.index(entry_state.number)] = 1.0
        #fundamental_vector = numpy.linalg.solve(transient_matrix, entry_state_vector)
        # probabilitity to be the exit state
        #self.probability_vector = fundamental_vector/recurrent_vector
        #Calculate mean residence time
        #self.mean_residence_time = numpy.sum(fundamental_vector)
        
        print recurrent_vector
        print transient_matrix
        fundamental_matrix = numpy.linalg.inv(transient_matrix)
        
        self.mean_residence_times = numpy.zeros(len(self.states))
        self.probability_matrix = numpy.zeros((len(self.states), len(self.states)))
        
        for i in range(len(self.states)):
            for j in range(len(self.states)):
                self.mean_residence_times[j] += fundamental_matrix[i][j]
                self.probability_matrix[j][i] = fundamental_matrix[i][j]/recurrent_vector[i]


   

    def write_data(self):
        logger.debug('saving data to %s' %self.path)
        f = open(self.path, 'w')
        print >> f, repr(self.state_numbers)
        print >> f, repr(self.mean_residence_times)
        print >> f, repr(self.probability_matrix)
        f.close()

    def read_data(self, get_state):
        logger.debug('reading data from %s' % self.path)
        f = open(self.path, 'r')
        self.state_numbers = eval(f.readline())
        self.states = [get_state(i) for i in self.state_numbers]
        
        matstr = ''
        for i in f.readlines():
            matstr += i
        self.probability_matrix = eval(matstr)

    def delete(self):
        logger.debug('deleting %s' % self.path)
        os.remove(self.path)
        self.states = None
        self.probability_matrix = None
        self.mean_residence_times = None

          



            
