import numpy

class Superbasin:
    """Class to manage super basin: calculate the mean residence time, exit probabilities, and perform Monte Carlo transitions out of the basin, """\
    """based on Novotny's Absorbing Markov Chain algorithm."""
    def __init__(self, statelist):
        #TODO: reinstate statelist!!!!
        self.nstates = len(statelist)
        self.states = statelist
        self._calculate_stuff()


    def pick_exit_state(self, entry_state):
       """Chosse an exit state (state of the basin from which we will be leaving) using absorbing Markov chain theory."""
       entry_state_index = self.states.index(entry_state)
       if entry_state_index is None:
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
       exit_state = self.states[exit_state_index]
       return time, exit_state


    def step(self, entry_state):
        """Perform a Monte Carlo transition: leave the basin."""\
        """The function returns a residence time as well as information to indenfity what saddle point was to leave the basin,"""\
        """from what state and to what state the system is moving to."""
        time, exit_state = self.pick_exit_state(entry_state)
        assert(time >= 0.0)

        # Make a rate table for all the exit state.  All processes are
        # needed as the might be a discrepancy in time scale
        # and it might be dangerous to weed out low rate events
        rate_table = []
        process_table = exit_state.get_rate_table()

        # Determine all process OUT of the superbasin
        for proc in process_table:
            if proc['product'] not in self.states:
                rate_table.append(proc)
        return rate_table, time, exit_state

    def contains_state(self, state):
        return state in self.states

    def _calculate_stuff(self):
        """Build the transient and recurrent matrices."""\
        """Calculate the fundamental matrix in order to be able to calculate the mean resisdence time"""\
        """and exit probablities any initial distribution."""

        recurrent_vector = numpy.zeros(self.nstates)
        transient_matrix= numpy.zeros((self.nstates, self.nstates))
        sum=0.0
        for i, item in enumerate(self.states):
            proc_table = item.get_rate_table()
            for process in proc_table:
                sum+=process['rate']
                if process['product'] not in self.states:
                    recurrent_vector[i] += process['rate']
                else:
                    #ouch that is complicated
                    j = self.states.index(process['product'])
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
