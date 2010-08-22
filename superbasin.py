import os
import math
import numpy
import logging
logger = logging.getLogger('superbasin')


class Superbasin:
 
    def __init__(self, path, kT, state_list = None, get_state = None):
        assert(isinstance(kT, float))
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
        for i in range(len(self.state_numbers)):
            if entry_state.number == self.state_numbers[i]:
                entry_state_index = i
                break
        else:
            raise ValueError('Passed entry state is not in this superbasin')

        probability_vector = self.probability_matrix.transpose()[entry_state_index]
        if abs(1.0-numpy.sum(probability_vector)) > 1e-3:
            logger.warning("the probability vector isn't close to 1.0")
            logger.warning('probability_vector ' + str(probability_vector) + " " + str(numpy.sum(probability_vector)))
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

    def test(self, entry_state):
        n=len(self.mean_residence_times)
        p=[0]*n
        m=10
        for i in range(m):
            time, exit_state_index=self.pick_exit_state(entry_state)
            p[exit_state_index]+=1
        for i in range(n):
            p[i]=float(p[i])/m
        print 'pick up p:', p
        return time, exit_state_index

    def step(self, entry_state, get_product_state):
        time, exit_state_index=self.pick_exit_state(entry_state)
        exit_state=self.states[exit_state_index]
        #make a rate table for the exit state
        rate_table = []
        ratesum = 0.0
        
        #XXX: Only process in the thermally accesbile window should be on this rate table
        process_table=exit_state.get_process_table()
        for proc_id in process_table:
            process = process_table[proc_id]
            if process['product'] not in self.state_numbers:
                rate_table.append([proc_id, process['rate']])
                ratesum += process['rate']
        
        p = 0.0
        u = numpy.random.random_sample()
        for i in range(len(rate_table)):
            p += rate_table[i][1]/ratesum
            if p>=u:
                exit_process_index = rate_table[i][0]
                break
        else:
            logger.warning("Warning: failed to select rate. p = " + str(p))
        absorbing_state = get_product_state(exit_state.number, exit_process_index)
        assert(time >= 0.0)
        return time, absorbing_state

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self): 
        recurrent_vector = numpy.zeros(len(self.states))
        transient_matrix= numpy.zeros((len(self.states), len(self.states)))
        sum=0.0
        for i in range(len(self.states)):
            proc_table = self.states[i].get_process_table()
            for process in proc_table.values():
                sum+=process['rate']
                if process['product']==-1 or process['product'] not in self.state_numbers: 
                    recurrent_vector[i] += process['rate']
                else:
                    j = self.state_numbers.index(process['product'])
                    transient_matrix[j][i] += process['rate']
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
        
        fundamental_matrix = numpy.linalg.inv(transient_matrix)
        n=len(fundamental_matrix)
        id=numpy.zeros((n, n))
        self.mean_residence_times = numpy.zeros(len(self.states))
        self.probability_matrix = numpy.zeros((len(self.states), len(self.states)))
        
        for i in range(len(self.states)):
            for j in range(len(self.states)):
                self.mean_residence_times[j] -= fundamental_matrix[i][j]
                self.probability_matrix[i][j] = -recurrent_vector[i]*fundamental_matrix[i][j]

        for i in self.probability_matrix.transpose():
            if abs(1-i.sum()) > 1e-3:
                logger.debug('Probability matrix has row which does not add up to 1')
                logger.debug('Row: %s' % str(i))
                logger.debug('Transient matrix:\n%s' % str(transient_matrix))
                logger.debug('Recurrent vector:\n%s' % str(recurrent_vector))
                logger.debug('Fundamental matrix:\n%s' % str(fundamental_matrix))
            


    def write_data(self):
        logger.debug('saving data to %s' %self.path)
        f = open(self.path, 'w')
        for i in [self.state_numbers, self.mean_residence_times, self.probability_matrix.ravel()]:
            for j in i:
                print >> f, repr(j), 
            print >> f
        f.close()


    def read_data(self, get_state):
        logger.debug('reading data from %s' % self.path)
        f = open(self.path, 'r')
        
        self.state_numbers = []
        for i in f.readline().rstrip().split():
            self.state_numbers.append(int(i)) 
        self.states = [get_state(i) for i in self.state_numbers]
        self.mean_residence_times = []
        for i in f.readline().rstrip().split():
            self.mean_residence_times.append(numpy.float64(i)) 
        pmat = [] 
        for i in f.readline().rstrip().split():
            pmat.append(numpy.float64(i))
        self.probability_matrix = numpy.array(pmat).reshape((len(self.states), len(self.states)))
        f.close()        

    def delete(self):
        logger.debug('deleting %s' % self.path)
        os.remove(self.path)
        self.states = None
        self.probability_matrix = None
        self.mean_residence_times = None
