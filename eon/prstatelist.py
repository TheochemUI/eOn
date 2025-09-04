
""" The statelist module. """

import logging
logger = logging.getLogger('statelist')

from eon import prstate
from eon import statelist


class PRStateList(statelist.StateList):
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, initial_state = None):
        statelist.StateList.__init__(self, prstate.PRState, initial_state)

    def register_process(self, reactant_number, product_number, process_id):
        # Get the reactant and product state objects.
        reactant = self.get_state(reactant_number)
        product = self.get_state(product_number)
        reactant.load_process_table()
        product.load_process_table()

        # Make the reactant process point to the product state number.
        reactant.procs[process_id]["product"] = product_number
        reactant.save_process_table()

        product.set_energy(reactant.procs[process_id]["product_energy"])
        # Update the reactant state to point at the new state id.
        reactant.procs[process_id]['product'] = product_number
        reactant.save_process_table()
