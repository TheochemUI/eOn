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
import config
import ConfigParser
import atomview
import atoms
import glob

import pylab as p
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as drawArea
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar


class akmcgui(atomview.atomview):



    
    def __init__(self):
        # imports from glade
        gladetree = gtk.glade.XML(os.path.join(pathfix.path, "tools/akmcgui.glade"))
        gui = gladetree.get_widget("akmcgui")
        atomview.atomview.__init__(self, gui)
        self.stateRB = gladetree.get_widget("stateRB")
        self.processesRB = gladetree.get_widget("processesRB")
        self.stateScale = gladetree.get_widget("stateScale")
        self.processesSB = gladetree.get_widget("processesSB")
        self.interpolationCB = gladetree.get_widget("interpulationCB")
        self.interpolationSB = gladetree.get_widget("interpolationSB")
        self.stateEnergy = gladetree.get_widget("stateEnergy")
        self.statePlayTB = gladetree.get_widget("statePlayTB")
        self.state_fpsSB = gladetree.get_widget("state_fpsSB")
        self.plotWindow = gladetree.get_widget("plotWindow")
        self.energy_plotButton = gladetree.get_widget("energy_plotButton")
        self.tempLabel = gladetree.get_widget("tempLabel")
        self.rateLabel = gladetree.get_widget("rateLabel")
        self.barrierLabel = gladetree.get_widget("barrierLabel")
        self.prefactorLabel = gladetree.get_widget("prefactorLabel")
        # defaults
        self.config = ConfigParser.SafeConfigParser()
        self.config.read(os.path.join(pathfix.path, "default_config.ini"))
        try:
            self.config.read("./config.ini")
        except:
            print "No config.ini found in local directory, using default values."
        self.processesSB.set_sensitive(False)
        self.interpolationCB.set_sensitive(False)
        self.interpolationSB.set_sensitive(False)
        self.changeImage()
        self.energy_changed()
        self.state_fpsSB.set_value(1)
        self.playImage = gtk.Image()
        self.playImage.set_from_stock(gtk.STOCK_MEDIA_PLAY, gtk.ICON_SIZE_BUTTON)
        self.pauseImage = gtk.Image()
        self.pauseImage.set_from_stock(gtk.STOCK_MEDIA_PAUSE, gtk.ICON_SIZE_BUTTON)
        self.statePlayTB.set_image(self.playImage)
        self.tempLabel.set_text(self.config.get("Main", "temperature"))
        self.processtable()
        self.plotImage = gtk.Image()
        self.plotImage.set_from_file("plot_icon.png")
        self.energy_plotButton.set_image(self.plotImage)
        # event handeling
        self.stateRB.connect("clicked", self.changeImage)
        self.processesRB.connect("clicked", self.changeImage)
        self.interpolationCB.connect("toggled", self.changeImage)
        self.interpolationSB.connect("value-changed", self.changeImage)
        self.stateScale.connect("value-changed", self.changeImage)
        self.processesSB.connect("value-changed", self.changeImage)
        self.stateRB.connect("clicked", self.RB_changed)
        self.processesRB.connect("clicked", self.RB_changed)
        self.interpolationCB.connect("clicked", self.interpolationCB_changed)
        self.stateScale.connect("value-changed", self.energy_changed)
        self.statePlayTB.connect("toggled", self.state_play)
        self.stateScale.connect("value-changed", self.processtable)
        self.state_fpsSB.connect("value-changed", self.state_play)
        self.energy_plotButton.connect("clicked", self.energy_plot)
        self.processesSB.connect("value-changed", self.processtable)
        
        
        

#
# Events-----------------------------------------------------
#

    def changeImage(self, *args):
        if self.stateRB.get_active() == True:
            states = glob.glob("./states/*")
            i = 0
            while ("./states/%d" %i) in states:
                i+=1    
            numStates = i-1
            self.stateScale.set_range(0, numStates)
            datapass = io.loadcons("states/%d/reactant.con" % self.stateScale.get_value())
        else: 
            if self.interpolationCB.get_active() == False:
                processes = glob.glob("./states/%d/procdata/*" % self.stateScale.get_value())
                j = 0
                while ("./states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), j)) in processes:
                    j+=1
                numProcesses = j-1
                if numProcesses == 0:
                    self.processesSB.set_value(0)
                    self.processesSB.set_sensitive(False)
                else:
                    self.processesSB.set_sensitive(True)            
                    self.processesSB.set_range(0, numProcesses)
                saddle = io.loadcon("states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), self.processesSB.get_value()))
                reactant = io.loadcon("states/%d/procdata/reactant_%d.con" %(self.stateScale.get_value(), self.processesSB.get_value()))
                product = io.loadcon("states/%d/procdata/product_%d.con" %(self.stateScale.get_value(), self.processesSB.get_value()))
                datapass = [reactant, saddle, product]
            else:
                processes = glob.glob("./states/%d/procdata/*" % self.stateScale.get_value())
                j = 0
                while ("./states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), j)) in processes:
                    j+=1
                numProcesses = j-1        
                if numProcesses == 0:
                    self.processesSB.set_value(0)
                    self.processesSB.set_sensitive(False)
                else:
                    self.processesSB.set_sensitive(True)            
                    self.processesSB.set_range(0, numProcesses)
                saddle = io.loadcon("states/%d/procdata/saddle_%d.con" %(self.stateScale.get_value(), self.processesSB.get_value()))
                reactant = io.loadcon("states/%d/procdata/reactant_%d.con" %(self.stateScale.get_value(), self.processesSB.get_value()))
                product = io.loadcon("states/%d/procdata/product_%d.con" %(self.stateScale.get_value(), self.processesSB.get_value()))
                datapass = [reactant, saddle, product]
                N = int (self.interpolationSB.get_text())+1
                p = [reactant, saddle, product]
                q = []
                for k in range(len(p)-1):
                    v = atoms.pbc(p[k+1].r - p[k].r, p[k].box)
                    d = np.linalg.norm(v)
                    v /=d
                    q.append(p[k].copy())
                    for l in range(1,N):
                        temp = p[k].copy()
                        temp.r += l * (v * (d/(N)))
                        q.append(temp)
                q.append(p[2].copy())
                datapass = q
        self.data_set(datapass)
        
        
    def RB_changed(self, widget, data=None):
        if self.stateRB.get_active() == True:
            self.processesSB.set_value(0)
            self.processesSB.set_sensitive(False)
            self.interpolationCB.set_sensitive(False)
            self.interpolationSB.set_sensitive(False)
            self.statePlayTB.set_sensitive(True)
            self.state_fpsSB.set_sensitive(True)
            return True
        else:
            self.processesSB.set_sensitive(True)
            self.interpolationCB.set_sensitive(True)
            self.statePlayTB.set_sensitive(False)
            self.state_fpsSB.set_sensitive(False)
            self.statePlayTB.set_active(False)
            if self.interpolationCB.get_active() == True:
                self.interpolationSB.set_sensitive(True)
            return True
           
            
    def interpolationCB_changed(self, widget, data=None):
        if self.interpolationCB.get_active() == True:
            self.interpolationSB.set_sensitive(True)
        else:  
            self.interpolationSB.set_sensitive(False)
        
        
    def energy_changed(self, *args):
        energy = open("states/state_table", "r")
        energyNumber = energy.readlines()
        a = energyNumber[int (self.stateScale.get_value())].split()[1]
        self.stateEnergy.set_markup("%f<b>eV</b>" % float (a))
        
    def processtable(self, *args):
        a = open("states/%s/processtable" % int(self.stateScale.get_value()) , "r" )
        a = a.readlines()
        b = a[int (self.processesSB.get_value()) +1].split()
        self.rateLabel.set_text(b[7])
        self.barrierLabel.set_text(b[6])
        self.prefactorLabel.set_text(b[2])
      
        
    def state_play(self, *args): 
        try:
            gobject.source_remove(self.timer_id)
        except:
            pass 
        if self.statePlayTB.get_active() == False:
            self.statePlayTB.set_image(self.playImage)
        else:
            self.statePlayTB.set_image(self.pauseImage)
            states = glob.glob("./states/*")
            i = 0
            while ("./states/%d" %i) in states:
                i+=1    
            numStates = i-1     
            def loop():    
                self.stateScale.set_value(int (self.stateScale.get_value())+1)
                if self.stateScale.get_value() >= numStates:
                    self.stateScale.set_value(0)
                return True 
            self.timer_id = gobject.timeout_add(1000/(int (self.state_fpsSB.get_value())), loop)
           
            
    def energy_plot(self, *args):
        a = open("states/state_table", "r")
        b = a.readlines()
        x = []
        y = []
        for i in range(len(b)):
            y.append(float (b[i].split()[1]))     
        for i in range(len(b)):
            x.append(i)
        self.plotWindow.set_default_size(600,600)    
        fig = p.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(x, y)
        p.xlabel("State")
        p.ylabel("eV")
        container = gtk.VBox()
        self.plotWindow.add(container)
        graph = drawArea(fig)
        toolbar = NavigationToolbar(graph, self.plotWindow)
        container.pack_start(graph)
        container.pack_start(toolbar, False, False)
        self.plotWindow.show_all()
                     

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
        
