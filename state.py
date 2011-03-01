##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

""" The state module. """

import os
import shutil
import math
import numpy

import io
import parsers
import tables
import atoms
import config
#import state

import logging
amkc_state_logger = logging.getLogger('state')


class State:
    """ The State super class. """
#    def __init__(self, path, nr, state_list, previous_state_num = -1, 
    def __init__(self, path, nr, state_list, 
                 reactant_path = None, energy = None):                 
        """ Creates a new State, with lazily loaded data. """

        # The parent statelist.
        self._state_list = state_list
        self._nr = nr

        # Lazily loaded data.
        self._info = None
        self._processtable = None

        # Data stored on the disk
        self._path = path        
        self._reactant_path = os.path.join(self._path, "reactant.con")

        # If this state does not exist on disk, create it.
        if not os.path.isdir(self._path):
            if reactant_path == None:
                raise IOError("State needs a reactant_path when it is being instantiated to disk.")
            os.mkdir(self._path)
            shutil.copy(reactant_path, self._reactant_path)

            # Creates new info file for the state
            self._info = parsers.StateInfo(self._path, self._nr)
            self._info.new()
#            self._info.new(previous_state_num)
            if energy is not None:
                self._set_energy(energy)
                
            self._info.save()

        self._procdata_path = os.path.join(self._path, "procdata")
        self._tar_path = os.path.join(self._path, "procdata.tar")   
        self._procdata_tarred = False
        # Create storage directory for the processes
        if not os.path.isdir(self._procdata_path):
            os.mkdir(self._procdata_path)

    def __repr__(self):
        return "State #%i" % self._nr
        
    def get_nr(self):
        return self._nr

    def get_path(self):
        return self._path

    def get_reactant(self):
        """ Loads the reactant.con into a point and returns it. """
        return io.loadcon(self._reactant_path)

    def tar_procdata(self):
        if not self._procdata_tarred:
            tar = tarfile.TarFile(self._tar_path, 'w')
            for i in os.listdir(self._procdata_path):
                tar.add(i)
                os.unlink(i)
            tar.close()
        else:
            amkc_state_logger.warning("Attempted to tar an already tarred procdata")
        self._procdata_tarred = True
        return


    #------------------
    # Info file related
    def _initialize_info(self):
        if self._info == None:
            self._info = parsers.StateInfo(self._path, self._nr)
        return
        
    def _set_energy(self, energy):
        """ Loads the info file if it has not been loaded and sets the reactant energy 
            variable. """
        self._initialize_info()
        self._info.set_reactant_energy(energy)

    def get_energy(self):
        """ Loads the info file if it is not already loaded and returns the energy, or None
            if it is not there. """
        self._initialize_info()
        return self._info.get_reactant_energy()

#    def get_previous_state(self):
#        self._initialize_info()
#        return self._info.get_previous_state()

    #----------------------
    # Process table related
    def _initialize_process_table(self):
        if self._processtable == None:
            self._processtable = tables.ProcessTable(self._path, self._nr)    
        return        


class AKMCState(State):
    def __init__(self, statepath, statenumber, statelist,
                 reactant_path = None, energy = None):                 
        """ Creates a new State, with lazily loaded data. """
        State.__init__(self, statepath, statenumber, statelist, 
                             reactant_path, energy)
        return


    def append_search_result(self, process, comment):
        tables.append_search_result(result, self._nr, comment)
        return

    def get_ratetable(self):
        """ Loads the process table if it has not been loaded and generates a rate table 
            according to kT and thermal_window. """
            
        self._initialize_process_table()
        lowest = self._processtable.get_lowest_barrier()
        energy_window = lowest + (config.kT * config.akmc_thermal_window)
        return self._processtable.get_ratetable(energy_window)

    def table_complete(self):
        table_comp = False
        if config.akmc_sampling_criteria == "confidence":
            if self.get_confidence() >= config.akmc_confidence:
                table_comp = True
        return table_comp
  

    def get_confidence(self):
        """ The confidence is a function of the ratio Nf/Ns, where Nf is the number of unique 
            processes and Ns is the number of searches performed.  
            
            When Nf or Ns are zero, there is no confidence.  Zero is returned in this case. 
            
            As the ratio Nf/Ns decreases, the confidence increases from 0.0 to a limit of 1.0.
            This is roughly equivalent to the statement: I am more confident that I have found
            all relevant processes if I have found 10 processes after 1000 searches than if I 
            have found 900 processes after 1000 searches.
            
            Nf is calculated to be the number of processes in the rate table. Ns is calculated
            to be the number of searches that resulted in a process on the rate table.
            
            When using recycling or kdb, it is useful to ignore processes that occur outside
            the hole, or the region in which the last process took place.  Focusing on this
            region means you can do possibly far fewer searches to reach confidence. When using
            the hole to filter processes, Nf and Ns only take into account processes that 
            intersect the hole. """

        self._initialize_process_table()
        lowest = self._processtable.get_lowest_barrier()
        energy_window = lowest + (config.kT * config.akmc_thermal_window)
        
        indices_relevant_processes = self._processtable.get_indices_processes_below(energy_window)
        
        unique_processes = len(indices_relevant_processes)
        total_repetitions = 0

        for index in indices_relevant_processes:
            total_repetitions += self._processtable.get_repeats(index)

        if unique_processes < 1 or total_repetitions < 1:
            confidence = 0.0
        else:
#            nr_u = float(unique_processes) # Nf
#            nr_t = float(total_repetitions) # Ns
            ratio = float(unique_processes) / float(total_repetitions) # Nf/Ns
            confidence = 1.0 + ratio * lambertw(-math.exp(-1.0 / ratio) / ratio)
        return confidence
#
#        prc = self.get_proc_repeat_count()
#        rt = self.get_ratetable()
#        Nf = 0.0
#        Ns = 0.0
#        conf = 0.0
#        #GH print "\nIn confidence"
#        for r in rt:
#            if self._state_list.filter_hole:
#                #GH print "filter_hole"
#                if self.proc_in_hole(r[0]):
#                    Ns += prc[r[0]]
#                    Nf += 1
#            else:
#                Ns += prc[r[0]]
#                Nf += 1
#            #GH print " proc:",r[0]," Nf:",Nf," dgen:",prc[r[0]]," Ns:",Ns
#        if Nf < 1 or Ns < 1:
#            conf = 0.0
#            #GH return 0.0
#        else:
#            conf = 1.0 + (Nf/Ns) * lambertw(-math.exp(-1.0 / (Nf/Ns))/(Nf/Ns))
#        #GH return 1.0 + (Nf/Ns) * lambertw(-math.exp(-1.0 / (Nf/Ns))/(Nf/Ns))
#        #GH print " conf: ",conf
#        return conf


    def proc_in_hole(self, procid):
        """ Returns True if the given process intersects the previous state/current state hole,
            False if not. """
        if self.get_previous_state() == -1:
            return True
        if procid in self.get_procs_in_hole():
            return True
        if procid in self.get_procs_not_in_hole():
            return False
        ps = self._state_list.get_state(self.get_previous_state()).get_reactant()
        cs = self.get_reactant()
        pan = atoms.per_atom_norm(ps.r - cs.r, ps.box)
        hole_atoms = []
        for i in range(len(pan)):
            if pan[i] > 0.2:
                hole_atoms.append(i)
        pp = self.get_process_product(procid)
        pan = atoms.per_atom_norm(pp.r - cs.r, cs.box)
        for i in range(len(pan)):
            if pan[i] > 0.2:
                if i in hole_atoms:
                    pih = self.get_procs_in_hole()
                    pih.append(procid)
                    self.set_procs_in_hole(pih)
                    return True
        pnih = self.get_procs_not_in_hole()
        pnih.append(procid)
        self.set_procs_not_in_hole(pnih)
        return False

    #------------------
    # Info file related 
    def get_lowest_barrier(self):
        self._initialize_info()
        return self._info.get_lowest_barrier()

    def get_unique_saddle_count(self):
        self._initialize_info()
        return self._info.get_unique_saddle_count()

    def get_good_saddle_count(self):
        self._initialize_info()
        return self._info.get_good_saddle_count()

    def get_bad_saddle_count(self):
        self._initialize_info()
        return self._info.get_bad_saddle_count()
        
    def get_proc_repeat_count(self):
        self._initialize_info()
        return self._info.get_proc_repeat_count()

    def get_total_saddle_count(self):
        self._initialize_info()
        return self._info.get_total_saddle_count()
    
    #----------------------
    # Process table related
    def _initialize_process_table(self):
        if self._processtable == None:
            self._processtable = tables.ProcessTableAKMC(self._path, self._nr)
        return 

    def register_good_process(self, process):
        self._initialize_info()
        self._initialize_process_table()

        # Forces that the reactant energy is fixed 
        resultdata = process["results"]
        if self._info.get_reactant_energy() is None:
            self._info.set_reactant_energy(resultdata["potential_energy_reactant"])
        else:
            resultdata["potential_energy_reactant"] = self._info.get_reactant_energy()

        # Adds the process to the process table and info file
        process_nr, comment, new_process, new_lowest_barrier = self._processtable.register_good_process(process)
        self._info.register_process(process_nr)
        if new_lowest_barrier:
            self._info.set_lowest_barrier(self._processtable.get_lowest_barrier())

        # Append to the search report file for this state 
        tables.append_search_result(process, self._nr, comment)        
        return process_nr

    def register_reverse_process(self, forward_reactant, forward_process, forward_process_nr):
        # This must be a new state
        self._initialize_info()
        self._initialize_process_table()
        
        reverse_process_nr, comment, new_process, new_lowest_barrier =\
            self._processtable.register_reverse_process(forward_reactant, forward_process, forward_process_nr)

        self._info.register_process(reverse_process_nr)
        if new_lowest_barrier:
            self._info.set_lowest_barrier(self._processtable.get_lowest_barrier())

        # Append to the search report file for this state 
        tables.append_search_result(forward_process, self._nr, comment)        
        return reverse_process_nr        

    def register_bad_process(self, process, store):
        self._initialize_info()
        self._info.increase_bad_saddle_count()

        comment = tables.register_bad_process(process, store)

        # Append to the search report file for this state
        tables.append_search_result(process, self._nr, comment)
        return

    def set_process_product_nr(self, proc_nr, prod_nr):
        self._initialize_process_table()
        return self._processtable.set_product_nr(proc_nr, prod_nr)

    def get_process_product_nr(self, nr):
        self._initialize_process_table()
        return self._processtable.get_product_nr(nr)

    def get_process_product_energy(self, nr):
        self._initialize_process_table()
        return self._processtable.get_product_energy(nr)

    def get_process_saddle_energy(self, nr):
        self._initialize_process_table()
        return self._processtable.get_saddle_energy(nr)

    def get_process(self, id):
        self._initialize_process_table()
        return self._processtable.get_process(id)

#    def get_reverse_process(self, id):
#        self._initialize_process_table()
#        return self._processtable.get_reverse_process(id)

    def get_num_procs(self):
        """ Loads the process table if it is not already loaded and returns the length of it """
        self._initialize_process_table()
        return self._processtable.get_unique_saddle_count()

    def get_process_table(self):
        self._initialize_process_table()
        return self._processtable.get_table()
    
    # Utility functions for loading process .con and mode files.
    
    def get_process_reactant(self, nr):
        self._initialize_process_table()
        return self._processtable.pass_process_reactant(nr)

    def get_process_saddle(self, nr):
        self._initialize_process_table()
        return self._processtable.pass_process_saddle(nr)

    def get_process_product(self, nr):
        self._initialize_process_table()
        return self._processtable.pass_process_product(nr)

#    def get_process_product(self, nr):
#        self._initialize_process_table()
#        return self._processtable.pass_process_product(nr)

    def get_process_mode(self, nr):
        self._initialize_process_table()
        return self._processtable.pass_process_mode(nr)

    # Utility functions for compiling procdata paths, whether the files exist or not.  
    def proc_reactant_path(self, nr):
        self._initialize_process_table()
        return self._processtable.proc_reactant_path(nr)

    def proc_saddle_path(self, nr):
        self._initialize_process_table()
        return self._processtable.proc_saddle_path(nr)

    def proc_product_path(self, nr):
        self._initialize_process_table()
        return self._processtable.proc_product_path(nr)

    def proc_mode_path(self, nr):
        self._initialize_process_table()
        return self._processtable.proc_mode_path(nr)        

    def proc_results_path(self, nr):
        self._initialize_process_table()
        return self._processtable.proc_results_path(nr)
 
    def proc_stdout_path(self, nr):
        self._initialize_process_table()
        return self._processtable.proc_stdout_path(nr)

    def get_lowest_barrier_excluding_product_nrs(self, nrs):
        self._initialize_process_table()
        return self._processtable.get_lowest_barrier_excluding_product_nrs(nrs)

    def get_saddle_energies_excluding_product_nrs(self, nrs):
        self._initialize_process_table()
        return self._processtable.get_saddle_energies_excluding_product_nrs(nrs)

   





# Lambert-W function: http://keithbriggs.info/software/LambertW.py
def lambertw(z):
    """Lambert W function, principal branch"""
    eps = 1.0e-12
    em1 = 0.3678794411714423215955237701614608
    if z < -em1:
        amkc_state_logger.error("Tried to evaluate Lambert W function @ < -1/e")
        raise ValueError()
    if 0.0 == z: 
        return 0.0
    if z < -em1 + 1e-4:
        q = z + em1
        r = math.sqrt(q)
        q2 = q * q
        q3 = q2 * q
        return\
         -1.0\
         +2.331643981597124203363536062168 * r\
         -1.812187885639363490240191647568 * q\
         +1.936631114492359755363277457668 * r * q\
         -2.353551201881614516821543561516 * q2\
         +3.066858901050631912893148922704 * r * q2\
         -4.175335600258177138854984177460 * q3\
         +5.858023729874774148815053846119 * r * q3\
         -8.401032217523977370984161688514 * q3 * q
    if z < 1.0:
        p = math.sqrt(2.0 * (2.7182818284590452353602874713526625 * z + 1.0))
        w = -1.0 + p * (1.0 + p * (-0.333333333333333333333 + p * 0.152777777777777777777777))
    else:
        w = math.log(z)
    if z > 3.0: w-=math.log(w)
    for i in xrange(10):
        e = math.exp(w)
        t = w * e - z
        p = w + 1.0
        t /= e * p - 0.5 * (p + 1.0) * t / p
        w -= t
        if abs(t) < eps * (1.0 + abs(w)): 
            return w
    amkc_state_logger.error("Failed to converge Lambert W function.")
    raise ValueError()
        
        
        
        
        
        
        
        
        
        
        
