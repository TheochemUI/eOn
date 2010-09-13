
import subprocess
import os
import shutil
import signal
import glob
import numpy
import logging
logger = logging.getLogger('kdb')    

import io



class KDB:


    def __init__(self, kdbpath, querypath, addpath):
        self.path = kdbpath
        self.querypath = querypath
        self.addpath = addpath


    def add_process(self, state, process_id):
        logger.debug("Adding process.")
        try:
            os.makedirs(self.path)
        except OSError:
            pass
        p = state.get_process_reactant(process_id)
        io.saveposcar(os.path.join(self.path, "REACTANT"), p)
        p = state.get_process_saddle(process_id)
        io.saveposcar(os.path.join(self.path, "SADDLE"), p)
        p = state.get_process_product(process_id)
        io.saveposcar(os.path.join(self.path, "PRODUCT"), p)
        mode = state.get_process_mode(process_id)
        f = open(os.path.join(self.path, "MODE"), 'w')
        for m in mode:
            f.write("%f %f %f\n" % (m[0], m[1], m[2]))
        f.close()
        sp = subprocess.Popen([self.addpath, "REACTANT", "PRODUCT", "SADDLE", "MODE"], 
                              cwd = self.path, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        sp.wait()
        return sp.communicate()[0].strip()
                

    def query(self, state, rhsco = 1.0, wait = False):
        # If the path already exists, remove it and create a new one.
        if os.path.isdir(self.path):
            # See if there is a PID file for a possibly already running query process.
            if os.path.exists(os.path.join(self.path, "PID")):
                # If so, try to kill it.
                f = open(os.path.join(self.path, "PID") , 'r')
                pid = int(f.readline())
                f.close()
                try:
                    os.kill(pid, signal.SIGKILL)
                except OSError:
                    # Probably wasn't running.
                    pass
            # Delete the scratch path.
            shutil.rmtree(self.path)        
        # Create the scratch path and save the POSCAR
        os.makedirs(self.path)
        p = state.get_reactant()
        io.saveposcar(os.path.join(self.path, "POSCAR"), p)
        sp = subprocess.Popen([self.querypath, "POSCAR", "--rhsco=%f" % rhsco], cwd = self.path, 
                stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        # Save the PID of this running process.
        f = open(os.path.join(self.path, "PID"), 'w')
        f.write("%d" % sp.pid)
        f.close()
        if wait:
            sp.wait()


    def make_suggestion(self, keep_path = None):
        if os.path.isdir(os.path.join(self.path, "kdbmatches")):
            dones = glob.glob(os.path.join(self.path, "kdbmatches",".done_*"))
            if len(dones) > 0:
                number = dones[0].split("_")[1]
                displacement = io.loadposcar(os.path.join(self.path, "kdbmatches", "SADDLE_%s" % number))
                mode = [[float(i) for i in l.strip().split()] for l in open(os.path.join(self.path, "kdbmatches", "MODE_%s" % number), 'r').readlines()[:]]
                if keep_path is not None:
                    if not os.path.isdir(keep_path):
                        os.mkdir(keep_path)
                    shutil.move(os.path.join(self.path, "kdbmatches", ".done_%s" % number), os.path.join(keep_path, ".done_%s" % number))
                    shutil.move(os.path.join(self.path, "kdbmatches", "SADDLE_%s" % number), os.path.join(keep_path, "SADDLE_%s" % number))
                    shutil.move(os.path.join(self.path, "kdbmatches", "MODE_%s" % number), os.path.join(keep_path, "MODE_%s" % number))
                else:
                    os.remove(os.path.join(self.path, "kdbmatches", ".done_%s" % number))
                    os.remove(os.path.join(self.path, "kdbmatches", "SADDLE_%s" % number))
                    os.remove(os.path.join(self.path, "kdbmatches", "MODE_%s" % number))
                return displacement, mode
        return None, None
                
            
    
















