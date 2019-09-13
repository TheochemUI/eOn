
import os
import logging
logger = logging.getLogger('superbasinscheme')

import sys

import config
import superbasin

class SuperbasinScheme:
    ''' This poorly-named class handles keeping track of which states belong
        to which superbasins, the SuperBasin object of those superbasins, and
        superbasining criteria. It also expands and merges superbasins'''

    def __init__(self, superbasin_path, states, kT):

        self.path = superbasin_path
        self.path_storage = os.path.join(superbasin_path, "storage")

        self.states = states
        self.kT = kT

        if not os.path.isdir(self.path):
            logger.warning('Superbasin path does not exist, creating %s' % self.path)
            os.makedirs(self.path)
            os.makedirs(self.path_storage)

        self.superbasins = []
        self.next_sb_num = 0
        for i in os.listdir(self.path):
            if i == 'storage':
                continue
            self.next_sb_num = max(self.next_sb_num, int(i))
            self.superbasins.append(
                superbasin.Superbasin(self.path, i,
                                      get_state=states.get_state)
            )

        self.next_sb_num += 1
        self.read_data()

    def get_containing_superbasin(self, state):
        for sb in self.superbasins:
            if sb.contains_state(state):
                return sb
        return None

    def make_basin_from_sets(self, start_state, end_state):
        # print "start_state: ",start_state
        # print "end_state: ",end_state

        start_states = set()
        sb = self.get_containing_superbasin(start_state)
        if sb is None:
            start_states.add(start_state)
        else:
            start_states.update(sb.states)

        end_states = set()
        sb = self.get_containing_superbasin(end_state)
        if sb is None:
            end_states.add(end_state)
        else:
            end_states.update(sb.states)

        # Is there an upper limit for the size of a superbasin?
        if config.sb_max_size:
            numstates = len(start_states) + len(end_states)
            # if number of states in the new superbasin will be larger
            # than the maximum size, do not proceed
            if numstates > config.sb_max_size:
                return

        # Now start creating the new superbasin.
        merge_states = [start_state, end_state]
        new_sb_states = set()
        for i in merge_states:
            sb = self.get_containing_superbasin(i)
            if sb is None:
                new_sb_states.add(i)
            else:
                new_sb_states.update(sb.states)
                # keep basins to analyze data
                if True:
                    sb.delete(self.path_storage)
                else:
                    sb.delete()
                self.superbasins.remove(sb)
        new_sb_states = list(new_sb_states)

        # print "connect_state: before"
        # self.states.connect_states(new_sb_states) #XXX:This should ensure detailed balance
        # However, it will likely be very slow. We should be able to do without it.
        # Also, if confidence is changed and new processes are found, the superbasin
        # will ignore these new processes.
        self.states.connect_state_sets(start_states, end_states)
        # print "connect_state: after"

        # print "superbasins.append: before"
        self.superbasins.append(
            superbasin.Superbasin(self.path, self.next_sb_num,
                                  state_list=new_sb_states)
        )
        # print "superbasins.append: after"
        logger.info("Created superbasin with states " + str([i.number for i in new_sb_states]))
        self.next_sb_num += 1

    def make_basin(self, merge_states):
        # Is there an upper limit for the size of a superbasin?
        if config.sb_max_size:
            # first determine how many states will be in the new superbasin
            numstates = 0
            for i in merge_states:
                sb = self.get_containing_superbasin(i)
                if sb is None:
                    numstates += 1
                else:
                    numstates += len(sb.states)
            # if number of states in the new superbasin will be larger
            # than the maximum size, do not proceed
            if numstates > config.sb_max_size:
                return
        # Now start creating the new superbasin.
        new_sb_states = set()
        for i in merge_states:
            sb = self.get_containing_superbasin(i)
            if sb is None:
                new_sb_states.add(i)
            else:
                new_sb_states.update(sb.states)
                # keep basins to analyze data
                if True:
                    sb.delete(self.path_storage)
                else:
                    sb.delete()
                self.superbasins.remove(sb)
        new_sb_states = list(new_sb_states)

        self.states.connect_states(new_sb_states) #XXX:This should ensure detailed balance
        # However, it will likely be very slow. We should be able to do without it.
        # Also, if confidence is changed and new processes are found, the superbasin
        # #will ignore these new processes.

        self.superbasins.append(
            superbasin.Superbasin(self.path, self.next_sb_num,
                                  state_list=new_sb_states)
        )
        logger.info("Created superbasin with states " + str([i.number for i in new_sb_states]))
        self.next_sb_num += 1

    def register_transition(self, start_state, end_state):
        raise NotImplementedError()

    def write_data(self):
        raise NotImplementedError()

    def read_data(self):
        raise NotImplementedError()

    def __del__(self):
        self.write_data()
    #def get_superbasins(self):
    #    if self.superbasins is not None:
    #        return self.superbasins
    #    else:
    #        dirs = os.listdir(self.path)
    #        for i in dirs:
    #            path = os.path.join(self.path, i)
    #            self.superbasins[int(i)] = superbasin.Superbasin(path)


class TransitionCounting(SuperbasinScheme):
    ''' Implements the transition counting scheme for superbasin detection '''

    def __init__(self, superbasin_path, states, kT, num_transitions):
        self.num_transitions = num_transitions
        SuperbasinScheme.__init__(self,superbasin_path, states, kT)

    def register_transition(self, start_state, end_state):
        logger.debug('Registering transitions')

        if start_state == end_state and not config.comp_use_identical:
            return

        start_count = self.get_count(start_state)
        if end_state not in start_count:
            start_count[end_state] = 0
        start_count[end_state] += 1

        if start_count[end_state] >= self.num_transitions:
            logger.debug( "Making basin ....")
            # self.make_basin([start_state, end_state])
            self.make_basin_from_sets(start_state, end_state)

    def write_data(self):
        logger.debug('writing')
        for start_state in self.count:
            data_path = os.path.join(start_state.path, config.sb_state_file)
            f = open(data_path, 'w')
            for end_state in self.count[start_state]:
                print(end_state.number, self.count[start_state][end_state], f)
            f.close()

    def read_data(self):
        self.count = {}

    def get_count(self, state):
        try:
            return self.count[state]
        except KeyError:
            data_path = os.path.join(state.path, config.sb_state_file)
            self.count[state] = {}
            if os.path.isfile(data_path):
                f = open(data_path, 'r')
                for i in f:
                    i = i.strip().split()
                    self.count[state][self.states.get_state(int(i[0]))] = int(i[1]) 
                f.close()
            return self.count[state]


class EnergyLevel(SuperbasinScheme):

    def __init__(self, superbasin_path, states, kT, energy_increment):
        self.energy_increment = energy_increment
        self.levels = {}
        SuperbasinScheme.__init__(self,superbasin_path, states, kT)

    def get_energy_increment(self, e_min, e_global_min, e_saddle):
        return self.energy_increment * (e_saddle - e_min) / (e_saddle - e_global_min)

    def get_statelist(self, state):
        sb = self.get_containing_superbasin(state)
        if sb:
            return sb.states
        else:
            return [ state ]

    # start_state and end_state are the non-sb ids. 
    def register_transition(self, start_state, end_state):
        '''Increments the energy level of the end state or sets it equal to the energy
           of the end_state if it hasn't been visited before.'''

        # error
        if start_state == end_state and not config.comp_use_identical:
        #if start_state == end_state:
            return

        # if the start state does not have an energy level yet, we set it to the energy of the state.
        if start_state not in self.levels:
            self.levels[start_state] = start_state.get_energy()

        # if the end state does not have an energy level yet, set it to the energy of the state.
        if end_state not in self.levels:
            self.levels[end_state] = end_state.get_energy()

        # determine wheter or not the start and end states belong to a superbasin.
        # if they belong to a superbasin, we need to know the ids of the states inside
        # that superbasin to determine the minimum energy state inside the superbasin.

        start_statelist = self.get_statelist(start_state)
        end_statelist   = self.get_statelist(end_state)

        e_min        = min ( state.get_energy() for state in end_statelist )
        e_global_min = min ( e_min , ( state.get_energy() for state in start_statelist ) )

        # determine the barrier
        # (LJ: What happens if there is more then one path between start_state and end_state???)
        barrier = 1e200
        proc_tab = start_state.get_process_table()
        for key in proc_tab:
            if proc_tab[key]['product'] == end_state.number:
                barrier = min(proc_tab[key]['barrier'], barrier)

        if barrier > 1e199:
            logger.warning("Start and end state have no direct connection")
            return
        saddle_energy = barrier + start_state.get_energy()

        # increment the level of the end_state. 
        for i in end_statelist:
            self.levels[i] += self.get_energy_increment(e_min, e_global_min, saddle_energy)

        # merge states if one of the two levels is higher than the saddle energy
        largest_level = max(self.levels[start_state], self.levels[end_state])

        if largest_level > saddle_energy:
            self.levels[start_state] = largest_level 
            self.levels[end_state] = largest_level 
            for i in end_statelist: #this may also be just [ end_state ]
                self.levels[i] = largest_level
            self.make_basin([start_state, end_state])

    def read_data(self):
        logger.debug('reading')
        for i in range(self.states.get_num_states()):
            state = self.states.get_state(i)
            data_path = os.path.join(state.path, config.sb_state_file)
            if os.path.isfile(data_path):
                f = open(data_path, 'r')
                self.levels[self.states.get_state(i)] = float(f.read().strip())
                f.close()

    def write_data(self):
        logger.debug('writing')
        for i in self.levels:
            data_path = os.path.join(i.path, config.sb_state_file)
            f = open(data_path, 'w')
            print("%f\n" % self.levels[i], f)
            f.close()

