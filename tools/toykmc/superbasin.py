import numpy
from state import State 

class Superbasin:
    """Class to manage super basin: calculate the mean residence time, exit probabilities, and perform Monte Carlo transitions out of the basin, """\
    """based on Novotny's Absorbing Markov Chain algorithm."""
    def __init__(self, path, id, statelist):
        #TODO: reinstate statelist!!!!
        self.nstates = len(State.states)
        self._calculate_stuff()


    def pick_exit_state(self, entry_state):
        """Chosse an exit state (state of the basin from which we will be leaving) using absorbing Markov chain theory."""
       entry_state_index = 
       if 
       raise ValueError('Passed entry state is not in this superbasin')

       probability_vector = self.probability_matrix.transpose()[entry_state_index]
       if abs(1.0-numpy.sum(probability_vector)) > 1e-3:
           print "the probability vector isn't close to 1.0"
           print 'probability_vector ' + str(probability_vector) + " " + str(numpy.sum(probability_vector))
       probability_vector /= numpy.sum(probability_vector)
       
       u = numpy.random.random_sample()
       p = 0.0
       for i in range(len(self.states)):
           p += probability_vector[i]
           if p>u:
               exit_state_index = i 
               break
       else:
           print "Warning: failed to select exit state. p = " + str(p)
       time = self.mean_residence_times[entry_state_index]
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

        # Determine all process OUT of the superbasin
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
            print "Warning: failed to select rate. p = " + str(p)
        
        # When requesting the product state the process 
        # gets added to the tables of events for both the forward 
        # and reverse process
        product_state = get_product_state(exit_state.number, exit_proc_id)

        return time, exit_state, product_state, exit_proc_id

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self):
        """Build the transient and recurrent matrices."""\
        """Calculate the fundamental matrix in order to be able to calculate the mean resisdence time"""\
        """and exit probablities any initial distribution."""

        recurrent_vector = numpy.zeros(self.nstates)
        transient_matrix= numpy.zeros((self.nstates, self.nstates))
        sum=0.0
        for i, item in enumerate(State.states.values()):
            proc_table = item.get_process_table()
            for process in proc_table.values():
                sum+=process['rate']
                #XXX: should be a hash check!!!
                if process['product'] not in State.states.values(): 
                    recurrent_vector[i] += process['rate']
                else:
                    #ouch that is complicated
                    j = State.states.index(State.gridhash(process['product'].grid))
                    transient_matrix[j][i] += process['rate']
                transient_matrix[i][i] -= process['rate']
        
        
        fundamental_matrix = numpy.linalg.inv(transient_matrix)
        self.mean_residence_times = numpy.zeros(len(self.states))
        self.probability_matrix = numpy.zeros((len(self.states), len(self.states)))
        
        for i in range(self.nstates):
            for j in range(self.nstates):
                self.mean_residence_times[j] -= fundamental_matrix[i][j]
                self.probability_matrix[i][j] = -recurrent_vector[i]*fundamental_matrix[i][j]

        for i in self.probability_matrix.transpose():
            if abs(1-i.sum()) > 1e-3:
                print "WARNING: Probability vector does not add up to 1" 


    
