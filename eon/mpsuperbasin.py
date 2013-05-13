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
import sys
import numpy
import logging
import mpmath
from mpmath import matrix
import pickle
import config
logger = logging.getLogger('mpsuperbasin')


class Superbasin:
    """Class to manage super basins: calculate the mean residence time, exit probabilities, and perform Monte Carlo transitions out of the basin, """\
    """based on Novotny's Absorbing Markov Chain algorithm."""

    def __init__(self, path, id, state_list = None, get_state = None):

        # Set the precision.
        mpmath.mp.dps = config.sb_arbitrary_precision

        #FIXME: self.states is literally a list of states, while in the superbasinscheme
        # self.states is a StateList object. Some renaming should happen.
        if state_list is None and get_state is None:
            raise ValueError('Superbasin must either have a list of states or a reference to get_state of a StateList')

        self.id = int(id)
        self.path = path+str(self.id)

        if not os.path.isfile(self.path):
            self.states = state_list
            self.state_numbers = [i.number for i in self.states]
            self._calculate_stuff()
            self.write_data()
        else:
            self.read_data(get_state)


    def pick_exit_state(self, entry_state):
        """Choose an exit state (state of the basin from which we will be leaving) using absorbing Markov chain theory."""
        for i in range(len(self.state_numbers)):
            if entry_state.number == self.state_numbers[i]:
                entry_state_index = i
                break
        else:
            raise ValueError('Passed entry state is not in this superbasin')

        probability_vector = self.probability_matrix.T[entry_state_index,0:]
        if abs(1.0-sum(probability_vector)) > 1e-3:
            logger.warning("MCAMC probability vector is not 1.0")
            logger.warning('Probability vector ' + str(probability_vector) + " " + str(sum(probability_vector)))
        probability_vector /= sum(probability_vector)

        u = numpy.random.random_sample()
        p = 0.0
        for i in range(len(self.states)):
            p += probability_vector[i]
            if p>u:
                exit_state_index = i 
                break
        else:
            logger.warning("Warning: Failed to select exit state; p = " + str(p))
        time = self.mean_residence_times[entry_state_index]
        return time, exit_state_index

    def test(self, entry_state):
        n=len(self.mean_residence_times)
        p=[0]*n
        m=10
        for i in range(m):
            time, exit_state_index=self.pick_exit_state(entry_state)
            p[exit_state_index]+=1
        for i in range(n):
            p[i]=float(p[i])/m
        print 'Pick up p:', p
        return time, exit_state_index

    def step(self, entry_state, get_product_state):
        """Perform a Monte Carlo transition: leave the basin."""\
        """The function returns a residence time as well as information to indenfity what saddle point was to leave the basin,"""\
        """from what state and to what state the system is moving to."""
        time, exit_state_index = self.pick_exit_state(entry_state)
        assert(time >= 0.0)
        exit_state = self.states[exit_state_index]

        # Make a rate table for all the exit state.  All processes are
        # needed as the might be a discrepancy in time scale
        # and it might be dangerous to weed out low rate events
        rate_table = []
        ratesum = 0.0
        process_table = exit_state.get_process_table()

        # Determine all process OUT of the superbasin (leaving from exit_state)
        for proc_id in process_table:
            process = process_table[proc_id]
            if process['product'] not in self.state_numbers:
                rate_table.append([proc_id, process['rate']])
                ratesum += process['rate']

        # picks the process to leave the superbasin
        p = 0.0
        u = numpy.random.random_sample()
        for i in range(len(rate_table)):
            p += rate_table[i][1]/ratesum
            if p>=u:
                exit_proc_id = rate_table[i][0]
                break
        else:
            logger.warning("Warning: Failed to select rate; p = " + str(p))
        
        # When requesting the product state the process
        # gets added to the tables of events for both the forward
        # and reverse process
        product_state = get_product_state(exit_state.number, exit_proc_id)

        return time, exit_state, product_state, exit_proc_id, self.id

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self):
        """Build the transient and recurrent matrices."""\
        """Calculate the fundamental matrix in order to be able to calculate the mean resisdence time"""\
        """and exit probablities any initial distribution."""

        # The i'th component of the recurrent vector contains the sum of all rates leaving state i (which 
        # is inside the composite) and entering a state which is not in the superbasin
        recurrent_vector = mpmath.matrix([[0]*len(self.states)])

        # The i'th diagonal component of the transient matrix contains minus the sum of _all_ rates 
        # from processes leaving state i to any other state (whether inside the composite or not).
        # the offdiagonal [i][j] components of the transient matrix contains the rate from state 
        # i (inside the superbasin) to state j (also in the superbasin)
        transient_matrix = mpmath.matrix(len(self.states))

        for i in range(len(self.states)):
            proc_table = self.states[i].get_process_table()
            for process in proc_table.values():

                # process is leaving the superbasin
                if process['product']==-1 or process['product'] not in self.state_numbers: 
                    recurrent_vector[i] += process['rate']

                # process remains in superbasin
                else:
                    j = self.state_numbers.index(process['product'])
                    # columns and rows interchanged as compared to theory?
                    transient_matrix[j,i] += process['rate']

                transient_matrix[i,i] -= process['rate']

        # Calculate mean residence time

        # Fundamental matrix is the inverse of the transient matrix T (not the inverse of (I-T) )
        fundamental_matrix = transient_matrix**-1


        # mean_residence_times contains the lifetime of state i in the composite state.
        self.mean_residence_times = mpmath.matrix([[0]*len(self.states)])
        # the probability matrix contains on the [i][j]'th position the probability of leaving 
        # the superbasin from state j, given that the system entered the superbasin from state i (or vice versa ;) )
        
        self.probability_matrix = mpmath.matrix(len(self.states))

        for i in range(len(self.states)):
            for j in range(len(self.states)):
                self.mean_residence_times[j] -= fundamental_matrix[i,j]
                self.probability_matrix[i,j] = -recurrent_vector[i] * fundamental_matrix[i,j]


        for i in range(self.probability_matrix.T.rows):
            row = self.probability_matrix.T[i,0:]
            if abs(1-sum(row)) > 1e-3:
                logger.debug('Probability matrix has row which does not add up to 1')
                logger.debug('Row: %s' % str(row))
                logger.debug('Transient matrix:\n%s' % str(transient_matrix))
                logger.debug('Recurrent vector:\n%s' % str(recurrent_vector))
                logger.debug('Fundamental matrix:\n%s' % str(fundamental_matrix))


    def write_data(self):
        logger.debug('Saving data to %s' %self.path)
        pickle.dump([repr(self.state_numbers), repr(self.mean_residence_times), repr(self.probability_matrix)], open(self.path, 'w'))

    def read_data(self, get_state):
        logger.debug('Reading data from %s' % self.path)
        data = [eval(d) for d in pickle.load(open(self.path, 'r'))]
        self.state_numbers = data[0]
        self.states = [get_state(i) for i in self.state_numbers]
        self.mean_residence_times = data[1]
        self.probability_matrix = data[2]

    def delete(self, storage=None):
        if storage is None:
            logger.debug('Deleting %s' % self.path)
            os.remove(self.path)
        else:
            logger.debug('Storing %s' % self.path)
            path_storage = storage+str(self.id)
            os.rename(self.path, path_storage)

        self.states = None
        self.probability_matrix = None
        self.mean_residence_times = None

