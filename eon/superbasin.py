
import os
import numpy
from eon.mcamc import mcamc
import logging
logger = logging.getLogger('superbasin')
from eon.config import config


class Superbasin:

    def __init__(self, path, id, state_list=None, get_state=None):
        """Initialize superbasin.

        Must pass either state_list or get_state:
        * state_list: Create a new state (must not exist, yet)
        * get_state: Read an existing state from path (must exist)

        """
        if (state_list is None) == (get_state is None):
            raise ValueError('Superbasin must either have a list of states or a reference to get_state of a StateList')
        self.id = int(id)
        self.path = os.path.join(path, str(self.id))
        # Get the states.
        if state_list is not None:
            if os.path.isfile(self.path):
                raise IOError("Superbasin file '%s' already exists!" % self.path)
            self.states = state_list
            self.state_numbers = [state.number for state in state_list]
            self.write_data()
        else: # get_state is not None
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
            for id, proc in list(procs.items()):
                if proc['product'] not in self.state_numbers:
                    st2col[(number, id)] = index
                    col2st[index] = (number, id)
                    index += 1

        # Build c.
        c = numpy.zeros(len(self.state_numbers))
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in list(procs.items()):
                c[st2i[number]] += proc['rate']

        # Build Q and R.
        Q = numpy.zeros((len(self.state_numbers), len(self.state_numbers)))
        R = numpy.zeros((len(self.state_numbers), len(col2st)))
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in list(procs.items()):
                if proc['product'] in self.state_numbers:
                    Q[st2i[number], st2i[proc['product']]] += proc['rate']
                else:
                    R[st2i[number], st2col[(number, id)]] += proc['rate']

        #lei debug
        print("################")
        print(str(self.id)+" c is: ")
        print(c)
        print(str(self.id)+" Q is: ")
        print(Q)
        print(str(self.id)+" R is: ")
        print(R)
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
        try:
            return state.number in self.state_dict
        except AttributeError:
            return False

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

    def _get_filtered_states(self):
        """Filter out states, which have no processes leading out of the superbasin
        and which have spent at least twice as much time in dynamics
        search as the longest time of all other states which do have
        processes leading out.

        Returns an iterator over states that should still be
        considered.

        """
        # Do not filter states if superbasin_confidence is disabled.
        if not config.sb_superbasin_confidence:
            for state in self.states:
                yield state
            return
        # The code for an enabled superbasin_confidence feature follows.
        if config.saddle_method == "dynamics":
            # Use time spent in dynamics.
            try:
                max_time = max(state.get_time()
                               for state in self.states
                               if state.get_ratetable(self))
            except ValueError:
                # There is no state with an exit process, yet.
                max_time = 0.0
            if max_time < 100.0:
                # If no state was searched for at least 100 fs, do not
                # filter yet.
                for state in self.states:
                    yield state
            for state in self.states:
                if state.get_ratetable(self) or state.get_time() < 2*max_time:
                    yield state
        else:
            # Use number of searches.
            try:
                max_searches = max(state.get_number_of_searches()
                                   for state in self.states
                                   if state.get_ratetable(self))
            except ValueError:
                # There is no state with an exit process, yet.
                max_searches = 0
            if max_searches < 1:
                # If no state was searched at least once, do not
                # filter yet.
                for state in self.states:
                    yield state
            for state in self.states:
                if state.get_ratetable(self) or state.get_number_of_searches() < 2*max_searches:
                    yield state

    def get_confidence(self):
        return min(state.get_confidence(self)
                   for state in self._get_filtered_states()) #self.states)

    def get_lowest_confidence_state(self):
        # If there are several states with equal confidence, order by
        # time spent searching.  This is needed in the case where we
        # know no exit at all from the superbasin. All states have
        # confidence 0.0 and we keep searching the same state. If that
        # state happens to have no exit but another one has, we will
        # never find anything. When we sort by time spent searching,
        # though, we will cycle the search through all states.
        if config.saddle_method == "dynamics":
            return sorted(self._get_filtered_states(),
                          key=lambda state: (state.get_confidence(self), state.get_time())
                          )[0]
        else:
            return sorted(self._get_filtered_states(),
                          key=lambda state: (state.get_confidence(self), state.get_number_of_searches())
                          )[0]
