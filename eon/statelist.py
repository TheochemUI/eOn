
""" The statelist module. """

import logging
logger = logging.getLogger('statelist')
import os

from eon import atoms
from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing


class StateList:
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, StateClass, initial_state = None, config: ConfigClass = EON_CONFIG):
        ''' Check to see if state_path exists and that state zero exists.
            Initializes state zero when passed a initial_state only if state
            zero doesn't already exist. '''

        self.config = config
        self.path = self.config.path_states
        self.epsilon_e = self.config.comp_eps_e
        self.epsilon_r = self.config.comp_eps_r
        self.use_identical = self.config.comp_use_identical
        self.StateClass = StateClass

        # Paths
        self.state_table_path = os.path.join(self.path, "state_table")

        # Create the statelist directory if it does not exist.
        if not os.path.isdir(self.path):
            logger.warning("State list path does not exist; Creating: %s" % self.path)
            os.makedirs(self.path)
            open(self.state_table_path, 'w').close()

        # Create the zero state directory if it does not exist.
        if not os.path.isdir(os.path.join(self.path, "0")):
            if initial_state is None:
                raise IOError("Missing zeroth state directory and no reactant provided")
            self.StateClass(
                statepath=os.path.join(self.path, "0"),
                statenumber=0,
                statelist=self,
                previous_state_num=-1,
                reactant_path=initial_state,
                config=config,
            )

        # Other class variables.
        self.states = {}
        self.state_table = None

    def get_num_states(self):
        """ Returns the number of lines in the state_table file. """
        self.load_state_table()
        return len(self.state_table)

    def get_product_state(self, state_number, process_id):
        ''' Returns a State object referenced by state_number and process_id. '''
        #TODO: Compare self.configuration of product with existing states.

        # If the number of states in state_table is zero
        #we need to add the zero state and energy to the state table.
        if self.get_num_states() == 0:
            zst = self.get_state(0)
            self.append_state_table(zst.get_energy())

        # Load the state object containing the process we want the product for.
        st = self.get_state(state_number)
        st.load_process_table()

        # Get the state number for the product.
        newstnr = st.procs[process_id]['product']

        # If the product id is not initialized, make sure it is not a copy of an existing state.
        # Otherwise, create it, connect it to st, and return it.
        if newstnr == -1:

            # Make a list of states for which we need to compare self.configurations.
            enew = st.procs[process_id]['product_energy']
            energetically_close = []
            for id in range(self.get_num_states()):
                if abs(self.get_state(id).get_energy() - enew) < self.epsilon_e:
                    energetically_close.append(id)

            # Perform distance checks on the energetically close self.configurations.
            if len(energetically_close) > 0:
                pnew = st.get_process_product(process_id)
                for id in energetically_close:
                    p = self.get_state(id).get_reactant()
                    if atoms.match(p, pnew, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, True):
                        if id == state_number:
                            logging.warning("State %i process %i leads back to initial state",
                                            state_number, process_id)
                        self.register_process(st.number, id, process_id)
                        return self.get_state(id)

            # The id for the new state is the number of states.
            newstnr = self.get_num_states()

            # Create the new state object.
            newst = self.StateClass(
                statepath=self.state_path(newstnr),
                statenumber=newstnr,
                statelist=self,
                previous_state_num=state_number,
                reactant_path=st.proc_product_path(process_id),
                config=self.config,
            )
            self.register_process(st.number, newstnr, process_id)

            # Append the new state to the state table.
            self.append_state_table(st.procs[process_id]['product_energy'])

        # The product state already exists, so get it.
        else:
            newst = self.get_state(newstnr)

        # Return the product state.
        return newst

    def get_state(self, state_number):
        ''' Returns a state object. '''
        if state_number in self.states:
            return self.states[state_number]
        st = self.StateClass(
            statepath=os.path.join(self.path, str(state_number)),
            statenumber=state_number,
            statelist=self,
            config=self.config,
        )
        self.states[state_number] = st
        return st

    def load_state_table(self):
        if self.state_table is None:
            f = open(self.state_table_path, 'r')
            lines = f.readlines()
            f.close()
            self.state_table = []
            for l in lines:
                self.state_table.append(float(l.strip().split()[1]))

    def save_state_table(self):
        if self.state_table != None:
            f = open(self.state_table_path, 'w')
            for i in range(len(self.state_table)):
                f.write("% 7d %16.5f\n" % (i, self.state_table[i]))
            f.close()

    def append_state_table(self, energy):
        number = self.get_num_states()
        f = open(self.state_table_path, 'a')
        f.write("% 7d %16.5f\n" % (number, energy))
        f.close()
        if self.state_table != None:
            self.state_table.append(energy)

    def state_path(self, state_number):
        """ Utility function to return the compiled path of a state, whether it exists or not. """
        return os.path.join(self.path, str(state_number))

if __name__ == "__main__":
    pass
