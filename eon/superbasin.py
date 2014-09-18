##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------


import os
import numpy
from mcamc import mcamc
import logging
logger = logging.getLogger('superbasin')


class Superbasin:

    def __init__(self, path, id, state_list = None, get_state = None):
        if state_list is None and get_state is None:
            raise ValueError('Superbasin must either have a list of states or a reference to get_state of a StateList')
        self.id = int(id)
        self.path = os.path.join(path, str(self.id))
        if not os.path.isfile(self.path):
            self.states = state_list
            self.state_numbers = [state.number for state in state_list]
            self.write_data()
        else:
            self.read_data(get_state)
        self.state_dict = {}
        for state in self.states:
            self.state_dict[state.number] = state

    def step(self, entry_state, get_product_state):
        # c_i (forming vector c) is the inverse of the sum of the rates for each transient state i
        # Q is the transient matrix of the canonical markov matrix. 
        # R is the recurrent matrix of the canonical markov matrix.
    
        # Build a mapping between transient states and row/column indices. Used in Q and R.
        st2i = {}
        i2st = {}
        index = 0
        for number in self.state_numbers:
            st2i[number] = index
            i2st[index] = number
            index += 1
            
        # Build a mapping between recurrent state identifiers [a (state number, process id) tuple] and 
        # column indices. Used only in R.
        st2col = {}
        col2st = {}
        index = 0
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in procs.iteritems():
                if proc['product'] not in self.state_numbers:
                    st2col[(number, id)] = index
                    col2st[index] = (number, id)
                    index += 1
        
        # Build c.
        c = numpy.zeros(len(self.state_numbers))
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in procs.iteritems():
                c[st2i[number]] += proc['rate']

        # Build Q and R.
        Q = numpy.zeros((len(self.state_numbers), len(self.state_numbers)))
        R = numpy.zeros((len(self.state_numbers), len(col2st)))
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in procs.iteritems():
                if proc['product'] in self.state_numbers:
                    Q[st2i[number], st2i[proc['product']]] += proc['rate']
                else:
                    R[st2i[number], st2col[(number, id)]] += proc['rate']
        
        # import pdb; pdb.set_trace()
        
        t, B, residual = mcamc(Q, R, c)
        logger.debug("residual %e" % residual)
        
        b = B[st2i[entry_state.number],:]
        p = 0.0
        u = numpy.random.sample()
        for i in range(len(b)):
            p += b[i]
            if p >= u:
                exit_state_number, exit_proc_id = col2st[i]
                break
        else:
            raise ValueError("Failed to select exit process.")
            
        exit_state = self.state_dict[exit_state_number]
        product_state = get_product_state(exit_state_number, exit_proc_id)        
        mean_time = t[st2i[entry_state.number]]
        return mean_time, exit_state, product_state, exit_proc_id, self.id

    def contains_state(self, state):
        return state in self.state_dict.values()

    def write_data(self):
        logger.debug('saving data to %s' %self.path)
        f = open(self.path, 'w')
        try:
            for number in self.state_numbers:
                f.write("%d " % number)
        finally:
            f.close()

    def read_data(self, get_state):
        logger.debug('reading data from %s' % self.path)
        f = open(self.path)
        try:
            self.state_numbers = [int(number)
                                  for number in f.read().split()]
        finally:
            f.close()
        self.states = [get_state(number) for number in self.state_numbers]

    def delete(self, storage=None):
        if storage is None:
            logger.debug('deleting %s' % self.path)
            os.remove(self.path)
        else:
            logger.debug('storing %s' % self.path)
            path_storage = os.path.join(storage, str(self.id))
            os.rename(self.path, path_storage)
        self.states = None

