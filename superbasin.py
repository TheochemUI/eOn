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


        exit_prob = self.prod_matrix[:,entry_state_index]
        u = numpy.random.random_sample()
        p = 0.0
        for i in range(len(exit_prob)):
            p += exit_prob[i]
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
        time = self.mean_residence_time

        if absorbing_state in self.states:
            #if we discovered a new connection between states in this basin
            #we need to update everything
            self._calculate_stuff()
            self.write_data()
        return time, absorbing_state

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self):
        '''easy to rip out later if necessary'''
       
        #Calculate the recurrent matrix
        recurrent_matrix = numpy.zeros((len(self.states), len(self.states)))
        for i in range(len(self.states)):
            proc_table = self.states[i].get_process_table()
            for key in proc_table:
                process = proc_table[key]
                if process['product']==-1 or process['product'] not in self.state_numbers: 
                    recurrent_matrix[i][i] += process['rate']
        self.recurrent_matrix = recurrent_matrix

        #Calculate the fundamental matrix
        fundamental_matrix = numpy.zeros((len(self.states), len(self.states)))
        for i in range(len(self.states)):
            for j in range(len(self.states)):
                if i == j:
                    #sum the rate table
                    for k in self.states[i].get_ratetable():
                        fundamental_matrix[i,i] += k[1]
                else:
                    #sum the rates from state j to state i
                    proc_table = self.states[j].get_process_table()
                    for key in proc_table:#Can I simply iterate over values?
                        proc = proc_table[key] 
                        if proc['product'] == self.state_numbers[i]:
                            fundamental_matrix[i,j] -= proc['rate']
        #XXX: Matrix inversion
        self.fundamental_matrix = numpy.linalg.inv(fundamental_matrix)
        
        self.prod_matrix = numpy.dot(self.recurrent_matrix, self.fundamental_matrix)
        
        #Calculate mean residence time
        self.mean_residence_time = 0.0
        for i in range(len(self.fundamental_matrix)):
            self.mean_residence_time += self.fundamental_matrix[i,i]

    
    def write_data(self):
        logger.debug('saving data to %s' %self.path)
        f = open(self.path, 'w')
        print >> f, repr(self.state_numbers)
        print >> f, repr(self.mean_residence_time)
        print >> f, repr(self.prod_matrix)
        f.close()

    def read_data(self, get_state):
        logger.debug('reading data from %s' % self.path)
        f = open(self.path, 'r')
        self.state_numbers = eval(f.readline())
        self.states = [get_state(i) for i in self.state_numbers]
        self.mean_residence_time = eval(f.readline())
        
        matstr = ''
        for i in f.readlines():
            matstr += i
        self.prod_matrix = eval(matstr)

    def delete(self):
        logger.debug('deleting %s' % self.path)
        os.remove(self.path)
        self.states = None
        self.prod_matrix = None
        self.mean_residence_time = None

          



            
