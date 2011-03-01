##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

""" The tables module. """

import os
import shutil
import math
import exceptions
import logging
logger = logging.getLogger('process_table')
akmclogger = logging.getLogger('process_table_akmc')

import atoms
import config
import io

class Table():
    def __init__(self, table_path):#, table_file_name):
        """ Create a new Table """
        
        self._table = None
        self._table_path = table_path
                
        self.initialize_table()
        self._changed = False

    def __del__(self):
        """ Saves the table to disk if it has been changed """
        if self._changed:
            self.save()
        return

    def __len__(self):
        length = 0
        try:
            length = len(self._table.keys())
        except:
            pass
        return length

    def initialize_table(self):
        """ If the process table is not loaded, load it and set the parameters, else nothing. """
        if self._table == None:
            self._table = {}
            # Only load data if file exsists
            if os.path.isfile(self._table_path):
                self.load()
        return
        
    def load(self):
        f = open(self._table_path)
        lines = f.readlines()
        f.close()
        data_cells = lines[0].strip().split()

        for l in lines[2:]:
            l = l.strip().split()
            # First element is the process nr
            nr = int(l[0])
            self._table[nr] = {}
            for j in xrange(1, len(data_cells)):
                try:
                    value = int(l[j])
                except exceptions.ValueError:
                    value = float(l[j])
                # Places the value in the table
                self._table[nr][data_cells[j]] = value
        return
        
    def save(self):
        """ If the processtable is present in memory, writes it to disk. """
        if self._table != None:
            format_int = "%10s"
            format_float = "%20s"

            f = open(self._table_path, 'w')
            
            key_table = self._table.keys()
            keys_data_cells = self._table[key_table[0]].keys()

            # Create header
            header = format_int % "nr"
            for key_data_cell in keys_data_cells:
                if type(self._table[key_table[0]][key_data_cell]) is int:
                    header += format_int % key_data_cell
                else:
                    header += format_float % key_data_cell
            header += "\n" + "-" * len(header) + "\n"
            f.write(header)
            
            # Write data
            for nr in key_table:
                data_line = format_int % nr
                for key_data_cell in keys_data_cells:
                    if type(self._table[nr][key_data_cell]) is int:
                        data_line += format_int % self._table[nr][key_data_cell]
                    else:
                        data_line += format_float % self._table[nr][key_data_cell]
                data_line += "\n"
                f.write(data_line)
                
            f.close()
        self._changed = False
        return

    # Pass data from this table        
    def get_table(self):
        return self._table

class StateTable(Table):
    def __init__(self, state_path):
        """ Create a new StateTable """
        Table.__init__(self, state_path)

    def add_state(self, nr, energy):
        self._table[nr] = {"energy" : energy}
        self._changed = True
        return 
    
    def get_energy_state_nr(self, nr):
        return self._table[nr]["energy"]
    
class ProcessTable(Table):
    def __init__(self, state_path, state_nr):
        """ Create a new ProcessTable """
        table_path = os.path.join(state_path, config.file_process_table)
        Table.__init__(self, table_path)
        
        self._state_nr = state_nr
#        self._procdata_path = state_path# os.path.join(state_path)
        self._procdata_path = os.path.join(state_path, "procdata")
        return
 
    # Pass data from this table                
    def get_process(self, nr):
        return self._table[nr]
        
    def get_process_nrs(self):
        return self._table.keys()
 
    # Utility functions for loading process con and mode files.
    def pass_process_reactant(self, nr):    
        return io.loadcon(self.proc_reactant_path(nr))
        
    def proc_reactant_path(self, nr):
        return os.path.join(self._procdata_path, "reactant_%d.con" % nr)
        
    def pass_process_product(self, nr):
        return io.loadcon(self.proc_product_path(nr))

    def proc_product_path(self, nr):
        return os.path.join(self._procdata_path, "product_%d.con" % nr)
        
        
class ProcessTableAKMC(ProcessTable):
    def __init__(self, state_path, state_nr):
        ProcessTable.__init__(self, state_path, state_nr)

    def register_good_process(self, process):
        """ Adds a process to the table and returns its nr, if it is a new process and if it is a new lowest barrier"""
        lowest_barrier = self.get_lowest_barrier()
        new_lowest_barrier = False
        new_process = False
                
        comment = None
        
        resultdata = process['results']

        barrier = resultdata["potential_energy_saddle"] - resultdata["potential_energy_reactant"]
        nr = self.known_process(process, barrier)
        
       # One more repetion
        if nr != -1:
            self._table[nr]['repeats'] += 1
            comment = "Process %i redetermined" % nr
        else:
            # This appears to be a unique process.
            # Check if the mode, reactant, saddle, and product are legit
            try:
                if 'mode' not in process:
                    io.load_mode(process['mode.dat'])
                if 'reactant' not in process:
                    io.loadcon(process['reactant.con'])
                if 'saddle' not in process:
                    io.loadcon(process['saddle.con'])
                if 'product' not in process:
                    io.loadcon(process['product.con'])
            except:
                akmclogger.exception("Mode, reactant, saddle, or product has incorrect format")
                return (None, comment, new_process, new_lowest_barrier)

            # The nr of the new process
            nr = len(self._table)
            
            if barrier < lowest_barrier:
                new_lowest_barrier = True
                akmclogger.info("found new lowest barrier %f for state %i", barrier, self._state_nr)
         
            rate = resultdata["prefactor_reactant_to_product"] * math.exp(-barrier / config.kT)

            # Append this mechanism to the process table.
            self._table[int(nr)] =  {"saddle_energy"     : resultdata["potential_energy_saddle"], 
                                     "prefactor"         : resultdata["prefactor_reactant_to_product"], 
                                     "product"           : -1, 
                                     "product_energy"    : resultdata["potential_energy_product"],
                                     "product_prefactor" : resultdata["prefactor_product_to_reactant"], 
                                     "barrier"           : barrier, 
                                     "rate"              : rate, 
                                     "repeats"           : 1}
            self._save_process(process, nr)
            new_process = True
            comment = "Process %i determined" % nr

        self._changed = True
        return nr, comment, new_process, new_lowest_barrier
 
    def register_reverse_process(self, forward_reactant, forward_process, forward_process_nr):
        lowest_barrier = self.get_lowest_barrier()
        new_lowest_barrier = False
        reverse_process_nr = len(self._table)

        reverse_barrier = forward_process["saddle_energy"] - forward_process["product_energy"]
        if reverse_barrier < lowest_barrier:
            new_lowest_barrier = True
            akmclogger.info("found new lowest barrier %f for state %i", reverse_barrier, self._state_nr)

        rate = forward_process["product_prefactor"] * math.exp(-reverse_barrier / config.kT)
        
        self._table[reverse_process_nr] = {"saddle_energy"     : forward_process["saddle_energy"], 
                                           "prefactor"         : forward_process["product_prefactor"], 
                                           "product"           : forward_reactant.get_nr(), 
                                           "product_energy"    : forward_reactant.get_energy(),
                                           "product_prefactor" : forward_process["prefactor"],
                                           "barrier"           : reverse_barrier, 
                                           "rate"              : rate, 
                                           "repeats"           : 1}

        self._copy_process(forward_reactant, forward_process_nr, reverse_process_nr)
        comment = "Process %i reverse" % reverse_process_nr
        self._changed = True
        
        return reverse_process_nr, comment, True, new_lowest_barrier

    def known_process(self, process, barrier):
        """Returns -1 if the process is unknown else the nr of the known process."""
        nr_same_sp = -1
        
        # Determine the number of processes in the process table that have a similar energy.
        energetically_close = []
        for nr in self._table.keys():
            if abs(self._table[nr]['barrier'] - barrier) < config.comp_eps_e:
                energetically_close.append(nr)

        # If the number of energetically similar saddles is > 0, we need to do distance checks on them.
        if len(energetically_close) > 0:
            #load the saddle
            process["saddle"] = io.loadcon(process["saddle.con"])
            p0 = process["saddle"]
            #ibox = numpy.linalg.inv(p0.box)
            for nr in energetically_close:
                p1 = self.pass_process_saddle(nr)
                if atoms.match(p1, p0, False):
                    nr_same_sp = nr
                    break
        return nr_same_sp

    def get_ratetable(self, energy_window=float("inf")):
        """ Generates a rate table according cointaining process within window. """
        ratetable = []
        for nr in self._table.keys():
            proc = self._table[nr]
            if proc['barrier'] < energy_window:
                ratetable.append((nr, proc['rate']))
        return ratetable

    def get_indices_processes_below(self, energy_window=float("inf")):
        """ Generates a rate table according cointaining process within window. """
        indices = []
        for nr in self._table.keys():
            proc = self._table[nr]
            if proc['barrier'] < energy_window:
                indices.append(nr)
        return indices

    # Pass data to this table
    def set_product_nr(self, proc_nr, prod_nr):
        self._table[proc_nr]['product'] = prod_nr
        self._changed = True
        return

    # Pass data from this table
    def get_lowest_barrier(self):
        return self.get_lowest_barrier_excluding_product_nrs([])[0]

    def get_lowest_barrier_excluding_product_nrs(self, exclude_nrs):
        """ Returns the lowest barrier and process nr for the processes excluding """
        """ the ones connected to the passed product numbers"""
        lowest_barrier = float("inf")
        lowest_barrier_nr = None
        lowest_barrier_product = None
        for key in self._table.keys():
            if not (self.get_product_nr(key) in exclude_nrs):
                barrier = self.get_barrier(key)
                if barrier < lowest_barrier:
                    lowest_barrier = barrier
                    lowest_barrier_nr = key
                    lowest_barrier_product = self.get_product_nr(key)
        return lowest_barrier, lowest_barrier_nr, lowest_barrier_product

    def get_saddle_energies_excluding_product_nrs(self, exclude_nrs):
        """ Returns the saddle energies for the processes excluding """
        """ the ones connected to the passed product numbers"""
        saddle_energies = []
        for key in self._table.keys():
            if not (self.get_product_nr(key) in exclude_nrs):
                saddle_energies.append(self.get_saddle_energy(key))
        return saddle_energies
        
    def get_unique_saddle_count(self):
        return len(self._table)

    def get_repeats(self, nr):
        return self._table[nr]['repeats']

    def get_product_nr(self, nr):
        return self._table[nr]['product']

    def get_product_energy(self, nr):
        return self._table[nr]['product_energy']
    
    def get_saddle_energy(self, nr):
        return self._table[nr]['saddle_energy']

    def get_barrier(self, nr):
        return self._table[nr]['barrier']
 
     # Utility functions for loading process con and mode files.
    def pass_process_saddle(self, nr):    
        return io.loadcon(self.proc_saddle_path(nr))

    def proc_saddle_path(self, nr):
        return os.path.join(self._procdata_path, "saddle_%d.con" % nr)

    def pass_process_mode(self, nr):
        return io.load_mode(self.proc_mode_path(nr))

    def proc_mode_path(self, nr):
        return os.path.join(self._procdata_path, "mode_%d.dat" % nr) 

    def proc_results_path(self, nr):
        return os.path.join(self._procdata_path, "results_%d.dat" % nr)        

    def proc_stdout_path(self, nr):
        return os.path.join(self._procdata_path, "stdout_%d.dat" % nr)        

    def _save_process(self, process, nr):
        # Save files to the procdata directory.
        open(self.proc_reactant_path(nr), 'w').writelines(process['reactant.con'].getvalue())
        open(self.proc_saddle_path(nr),   'w').writelines(process['saddle.con'].getvalue())
        open(self.proc_product_path(nr),  'w').writelines(process['product.con'].getvalue())
        open(self.proc_mode_path(nr),     'w').writelines(process['mode.dat'].getvalue())
        open(self.proc_results_path(nr),  'w').writelines(process['results.dat'].getvalue())

    def _copy_process(self, reactant, passed_nr, nr):
        # Copy files to the procdata directory.
        shutil.copy(reactant.proc_saddle_path(passed_nr),   self.proc_saddle_path(nr))
        shutil.copy(reactant.proc_reactant_path(passed_nr), self.proc_product_path(nr))
        shutil.copy(reactant.proc_product_path(passed_nr),  self.proc_reactant_path(nr))
        shutil.copy(reactant.proc_results_path(passed_nr),  self.proc_results_path(nr))
        shutil.copy(reactant.proc_mode_path(passed_nr),     self.proc_mode_path(nr))
        



def register_bad_process(process, store = False):
    """ Registers a bad saddle. """
    result_state_code = ["Good",
                         "Init",
                         "Saddle Search No Convex Region",
                         "Saddle Search Terminated Barrier",
                         "Saddle Search Terminated Total Iterations",
                         "Saddle Search Terminated Concave Iterations",
                         "Not Connected",
                         "Bad Prefactor",
                         "Bad Barrier",
                         "Minimum Not Converged"]
    if store:
        bad_procdata_path = os.path.join(config.path, "badprocesses")
        if not os.path.isdir(bad_procdata_path):
            os.mkdir(bad_procdata_path)
        open(os.path.join(bad_procdata_path, "reactant_%d.con" % process['wuid']), 'w').writelines(process['reactant.con'].getvalue())
        open(os.path.join(bad_procdata_path, "product_%d.con" % process['wuid']), 'w').writelines(process['product.con'].getvalue())
        open(os.path.join(bad_procdata_path, "mode_%d.dat" % process['wuid']), 'w').writelines(process['mode.dat'].getvalue())
        open(os.path.join(bad_procdata_path, "results_%d.dat" % process['wuid']), 'w').writelines(process['results.dat'].getvalue())
        open(os.path.join(bad_procdata_path, "saddle_%d.con" % process['wuid']), 'w').writelines(process['saddle.con'].getvalue())
    return result_state_code[process["results"]["termination_reason"]]

def append_search_result(process, state_nr, comment):
    f = None

    path_list_search_results = os.path.join(config.path_states, str(state_nr), "search_results.txt")

    if os.path.isfile(path_list_search_results):
        f = open(path_list_search_results, 'a')
    else:
        f = open(path_list_search_results, 'w')
        header = "%8s %10s %10s %10s %10s %10s %10s    %s\n" % ("wuid", "type",
             "barrier",
             "max-dist", "sad-fcs",
             "mins-fcs", "pref-fcs",
             "result             ")
        header += "-" * len(header) + '\n'
        f.write(header)
    try:
        r = process['results']
        f.write("%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (
            process["wuid"], process["type"], 
            r["potential_energy_saddle"] - r["potential_energy_reactant"],
            r["displacement_saddle_distance"], r["force_calls_saddle"] ,
            r["force_calls_minimization"], r["force_calls_prefactors"],
            comment))
    except:
        f.write("%8d %10s %10.5f %10.5f %10d %10d %10d    %s\n" % (
            0, "reverse", process["barrier"], 0, 0, 0, 0, comment))
    f.close()
    return