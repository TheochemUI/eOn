import os
import logging
logger = logging.getLogger('superbasinscheme')

import superbasin

class SuperbasinScheme:
    ''' This poorly-named class handles keeping track of which states belong
    to which superbasins, the SuperBasin object of those superbasins, and 
    superbasining criteria. It also expands and merges superbasins'''

    def __init__(self, superbasin_path, states, kT):
        self.path = superbasin_path
        self.states = states
        self.kT = kT

        if not os.path.isdir(self.path):
            logger.warning('superbasin path does not exist, creating %s' % self.path)
            os.makedirs(self.path)
        
        self.superbasins = []
        dirs = os.listdir(self.path)
        self.next_sb_num = 0
        for i in dirs:
            self.next_sb_num = max(self.next_sb_num, int(i))
            path = os.path.join(self.path, i)
            self.superbasins.append(superbasin.Superbasin(path, self.kT, get_state = states.get_state))
        self.next_sb_num += 1
        self.read_data()

    def get_containing_superbasin(self, state):
        for i in self.superbasins:
            if i.contains_state(state):
                return i
        return None

    def make_basin(self, merge_states):
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
                sb.delete()
                self.superbasins.remove(sb)
        
        
        self.states.connect_states(new_sb_states) #XXX:This should ensure detailed balance
        #However, it will likely be very slow. We should be able to do without it.
        #Also, if confidence is changed and new processes are found, the superbasin
        #will ignore these new processes.

        new_sb_path = os.path.join(self.path, str(self.next_sb_num))
        self.superbasins.append(superbasin.Superbasin(new_sb_path, self.kT, state_list = new_sb_states)) 
        
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
    ''' Implements the transition counting scheme for superbasining '''
    
    def __init__(self, superbasin_path, states, kT, num_transitions):
        self.num_transitions = num_transitions
        SuperbasinScheme.__init__(self,superbasin_path, states, kT)

    def register_transition(self, start_state, end_state):
        logger.debug('Registering transitions')
        
        if start_state == end_state:
            return
        if start_state not in self.count:
            self.count[start_state] = {}
            table_path = os.path.join(start_state.path)
            if os.path.isfile(table_path):
                f = open(table_path, 'r')
                for i in f:
                    i = i.strip().split()
                    self.count[start_state][self.states.get_state(int(i[0]))] = i[1]
                f.close()
        if end_state not in self.count[start_state]:
            self.count[start_state][end_state] = 0
        self.count[start_state][end_state] += 1
        
        if self.count[start_state][end_state] >= self.num_transitions:
            logger.debug( "Making basin....")
            self.make_basin([start_state, end_state])

    

    def write_data(self):
        logger.debug('writing')
        for start_state in self.count:
            data_path = os.path.join(start_state.path, 'superbasin_tc')
            f = open(data_path, 'w')
            for end_state in self.count[start_state]:
                print >> f, end_state.number, self.count[start_state][end_state]
            f.close()

    def read_data(self):
        self.count = {}
        logger.debug('reading')
        for i in range(self.states.get_num_states()):
            #TODO: Should each scheme have a unique file?
            state = self.states.get_state(i)
            data_path = os.path.join(state.path, 'superbasin_tc')
            if os.path.isfile(data_path):
                f = open(data_path, 'r')
                self.count[state] = {}
                for i in f:
                    i = i.strip().split()
                    self.count[state][self.states.get_state(int(i[0]))] = int(i[1]) 
                f.close()


class EnergyLevel(SuperbasinScheme):

    def __init__(self, superbasin_path, states, kT, energy_increment):
        self.energy_increment = energy_increment
        SuperbasinScheme.__init__(self,superbasin_path, states, kT)
        self.levels = {}

    def get_energy_increment(self):
        # implement variable increment according to JC paper
        return self.energy_increment

    def register_transition(self, start_state, end_state):
        '''Increments the energy level of the end state or sets it equal to the energy 
           of the end_state if it hasn't been visited before.'''
        if start_state == end_state:
            return
        sb = self.get_containing_superbasin(end_state)
        if not sb:
            up_states = [end_state]
        else:
            up_states = sb.states
        for i in up_states:
            if i not in self.levels:
                self.levels[i] = i.get_energy()
            else:
                self.levels[i] += self.energy_increment()

        #saddle energy is the total energy of the saddle
        largest_level = max(self.levels[start_state], self.levels[end_state])
       
        barrier = 0.0
        proc_tab = start_state.get_process_table()
        for key in proc_tab:
            if proc_tab[key][product] == end_state.number:
                barrier = min(proc_tab[key]['barrier'], barrier)
        
        if barrier == 0.0:
            logger.warning("Start and end state have no connections")
            return
        
        saddle_energy = barrier + start_state.get_energy()
        if largest_level > saddle_energy:
            self.levels[start_state] = largest_level 
            self.levels[end_state] = largest_level 
            self.merge_basin([start_state, end_state])


    def read_data(self):
        logger.debug('reading')
        for i in range(self.states.get_num_states()):
            state = self.states.get_state(i)
            data_path = os.path.join(state.path, 'superbasin_el')
            if os.path.isfile(data_path):
                f = open(data_path, 'r')
                self.levels[self.states.get_state(i)] = float(f.read().strip())
                f.close()

    def write_data(self):
        logger.debug('writing')
        for i in self.levels:
            data_path = os.path.join(start_state.path, 'superbasin_el')
            f = open(data_path, 'w')
            print >> f, "%f\n" % self.levels[i]
            f.close()

    
        
        

