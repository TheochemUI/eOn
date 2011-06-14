#!/usr/bin/env python

import os
import math
import time

import gtk
import gtk.gdk as gdk
import gtk.glade as glade
import gobject

import numpy as np

import pathfix
import atomview
import atoms
import glob

class akmcgui(atomview.atomview):



    
    def __init__(self):
        # Glade
        gladetree = gtk.glade.XML(os.path.join(pathfix.path, "tools/akmcgui.glade"))
        gui = gladetree.get_widget("akmcgui")
        atomview.atomview.__init__(self, gui)
        self.stateRB = gladetree.get_widget("stateRB")
        self.processesRB = gladetree.get_widget("processesRB")
        self.stateScale = gladetree.get_widget("stateScale")
        self.processesScale = gladetree.get_widget("processesScale")
        self.interpolationCB = gladetree.get_widget("interpulationCB")
        self.interpolationSB = gladetree.get_widget("interpolationSB")
        #defaulting options
        self.processesScale.set_sensitive(False)
        self.interpolationCB.set_sensitive(False)
        self.interpolationSB.set_sensitive(False)
        self.state_scale_changed()
        interpolation = []
        #event handling
        self.stateRB.connect("clicked", self.RB_changed)
        self.processesRB.connect("clicked", self.RB_changed)
        self.processesRB.connect("clicked", self.interpolation) 
        self.processesRB.connect("clicked", self.processes_scale_changed)
        self.interpolationCB.connect("toggled", self.interpolation)
        self.interpolationSB.connect("value-changed", self.interpolation)
        self.interpolationCB.connect("toggled", self.interpolationCB_changed)
        
#
# StateScale--------------------------------------------------        
#
        states = glob.glob("./states/*")
        
        i = 0
        while ("./states/%d" %i) in states:
            i+=1    
        numStates = i-1
        
        self.stateScale.set_range(0, numStates)       
        self.stateScale.connect("value-changed", self.state_scale_changed)

#
# ProcessesScale---------------------------------------------- 
#
        processes = glob.glob("./states/%d/procdata/*" % self.stateScale.get_value())
        
        j = 0
        while ("./states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), j)) in processes:
            j+=1
        numProcesses = j-1
        
        self.processesScale.set_range(0, numProcesses)
        self.processesScale.connect("value-changed", self.processes_scale_changed)
        self.processesScale.connect("value-changed", self.interpolation)
        

#     
# Events------------------------------------------------------        
#
    
    def RB_changed(self, widget, data=None):
        if self.stateRB.get_active() == True:
            self.processesScale.set_value(0)
            self.processesScale.set_sensitive(False)
            self.interpolationCB.set_sensitive(False)
            self.interpolationSB.set_sensitive(False)
            self.interpolationCB.set_active(False)
            state = io.loadcons("states/%d/reactant.con" % self.stateScale.get_value())
            self.data_set(state)
            return True
        else:
            self.processesScale.set_sensitive(True)
            self.interpolationCB.set_sensitive(True)
            return True
            
     
    def state_scale_changed(self, *args):
        state = io.loadcons("states/%d/reactant.con" % self.stateScale.get_value())
        self.data_set(state)
        return True
        
    def processes_scale_changed(self, *args):
        saddle = io.loadcon("states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), self.processesScale.get_value()))
        reactant = io.loadcon("states/%d/procdata/reactant_%d.con" %(self.stateScale.get_value(), self.processesScale.get_value()))
        product = io.loadcon("states/%d/procdata/product_%d.con" %(self.stateScale.get_value(), self.processesScale.get_value()))
        movie = [reactant, saddle, product]
        self.data_set(movie)
        return True
        
    def interpolationCB_changed(self, widget, data=None):
        if self.interpolationCB.get_active() == True:
            self.interpolationSB.set_sensitive(True)
        else:  
            self.interpolationSB.set_sensitive(False)
        
    def interpolation(self, widget, data=None):
        saddle = io.loadcon("states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), self.processesScale.get_value()))
        reactant = io.loadcon("states/%d/procdata/reactant_%d.con" %(self.stateScale.get_value(), self.processesScale.get_value()))
        product = io.loadcon("states/%d/procdata/product_%d.con" %(self.stateScale.get_value(), self.processesScale.get_value()))
        movie = [reactant, saddle, product]
        if self.interpolationCB.get_active() == True:
            N = int (self.interpolationSB.get_text())
            p = [reactant, saddle, product]
            q = []
            c = 1
            for k in range(len(p) -1):
                v = atoms.pbc(p[k+1].r - p[k].r, p[k].box)
                d = np.linalg.norm(v)
                v /= d
                q.append(p[k].copy())
                for l in range(1,N):
                    temp = p[k].copy()
                    temp.r += l * (v * (d/(N)))
                    q.append(temp)
            q.append(p[2].copy())
            self.data_set(q)
        else: 
            self.data_set(movie)              
    
        
 #       
 # Main-------------------------------------------------------       
 #   

     
if __name__ == "__main__":
    print os.getcwd()
    pid = os.fork()
    if pid:
        os._exit(0)
    import io
    import sys
    q = akmcgui()
    if len(sys.argv) > 1:
        data = io.loadposcars(sys.argv[1])
        if len(data) < 1:
            data = io.loadcons(sys.argv[1])
        q.data_set(data)
    gtk.main()

        
        

