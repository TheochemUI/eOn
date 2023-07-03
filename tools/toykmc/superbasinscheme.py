import superbasin

class SuperbasinScheme:
    ''' This poorly-named class handles keeping track of which states belong
    to which superbasins, the SuperBasin object of those superbasins, and
    superbasining criteria. It also expands and merges superbasins'''

    def __init__(self):
        self.superbasins = []

    def get_containing_superbasin(self, state):
        for i in self.superbasins:
            if i.contains_state(state):
                return i
        return None

    def make_basin(self, merge_states):
        print "ohdeargod make a basin"
        new_sb_states = []
        for i in merge_states:
            sb = self.get_containing_superbasin(i)
            if sb is None:
                if i not in new_sb_states:
                    new_sb_states.append(i)
            else:
                for j in sb.states:
                    if j not in new_sb_states:
                        new_sb_states.append(j)
                self.superbasins.remove(sb)

        #self.states.connect_states(new_sb_states) #XXX:This should ensure detailed balance
        #However, it will likely be very slow. We should be able to do without it.
        #Also, if confidence is changed and new processes are found, the superbasin
        #will ignore these new processes.

        self.superbasins.append(superbasin.Superbasin(new_sb_states))

        print "Created superbasin with states " #+ str([i.number for i in new_sb_states])

    def register_transition(self, start_state, end_state):
        raise NotImplementedError()

    def write_data(self):
        raise NotImplementedError()

    def read_data(self):
        raise NotImplementedError()

class TransitionCounting(SuperbasinScheme):
    ''' Implements the transition counting scheme for superbasining '''

    def __init__(self, num_transitions):
        self.num_transitions = num_transitions
        SuperbasinScheme.__init__(self)
        self.count = {}

    def register_transition(self, start_state, end_state):
        if start_state == end_state:
            return

        start_count = self.get_count(start_state)
        if end_state not in start_count:
            start_count[end_state] = 0
        start_count[end_state] += 1

        if start_count[end_state] >= self.num_transitions:
            self.make_basin([start_state, end_state])

    def get_count(self, state):
        try:
            return self.count[state]
        except:
            self.count[state] = {}
            return self.count[state]
