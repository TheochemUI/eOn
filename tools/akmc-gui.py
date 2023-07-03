#!/usr/bin/env python
##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

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
import configparser
import atomview
import atoms
import glob
import pylab as p
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as drawArea
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import importlib


class akmcgui(atomview.atomview):


    def __init__(self):
        # imports from glade
        gladetree = gtk.glade.XML(os.path.join(pathfix.path, "tools/akmc-gui.glade"))
        gui = gladetree.get_widget("akmcgui")
        atomview.atomview.__init__(self, gui)
        self.stateScale = gladetree.get_widget("stateScale")
        self.processesSB = gladetree.get_widget("processesSB")
        self.interpolationCB = gladetree.get_widget("interpulationCB")
        self.interpolationSB = gladetree.get_widget("interpolationSB")
        self.stateEnergy = gladetree.get_widget("stateEnergy")
        self.statePlayTB = gladetree.get_widget("statePlayTB")
        self.plotWindow = gladetree.get_widget("plotWindow")
        self.energy_plotButton = gladetree.get_widget("energy_plotButton")
        self.hbox = gladetree.get_widget("hbox5")
        self.tempLabel = gladetree.get_widget("simtemp")
        self.fileButton = gladetree.get_widget("fileButton")
        self.fileEntry = gladetree.get_widget("locationEntry")
        self.stateSB = gladetree.get_widget("stateSB")
        # defaults
        button_size = (20, 20)
        self.open = gtk.Image()
        self.open.set_from_stock(gtk.STOCK_FILE, gtk.ICON_SIZE_BUTTON)
        self.fileButton.set_image(self.open)
        self.stateScale.set_size_request(400,-1)
        self.interpolationSB.set_sensitive(False)
        self.playImage = gtk.Image()
        self.playImage.set_from_stock(gtk.STOCK_MEDIA_PLAY, gtk.ICON_SIZE_BUTTON)
        self.pauseImage = gtk.Image()
        self.pauseImage.set_from_stock(gtk.STOCK_MEDIA_PAUSE, gtk.ICON_SIZE_BUTTON)
        self.statePlayTB.set_image(self.playImage)
        self.plotImage = gtk.Image()
        image = gdk.pixbuf_new_from_file_at_size(os.path.join(pathfix.path, "tools/plot_icon.png"), 15, 15)
        self.plotImage.set_from_pixbuf(image)
        self.energy_plotButton.set_image(self.plotImage)
        self.directory = "./config.ini"
        self.fileEntry.set_text(self.directory)
        self.fileButton.set_sensitive(False)
        self.fileEntry.set_sensitive(False)
        self.startup()
        # event handeling
        self.interpolationCB.connect("toggled", self.changeImage)
        self.interpolationSB.connect("value-changed", self.changeImage)
        self.stateScale.connect("value-changed", self.changeImage)
        self.processesSB.connect("value-changed", self.changeImage)
        self.stateSB.connect("value-changed", self.stateSB_changed)
        self.interpolationCB.connect("clicked", self.interpolationCB_changed)
        self.stateScale.connect("value-changed", self.energy_changed)
        self.stateScale.connect("value-changed", self.state_play)
        self.statePlayTB.connect("toggled", self.state_play)
        self.stateScale.connect("value-changed", self.processtable)
        self.energy_plotButton.connect("clicked", self.energy_plot)
        self.stateScale.connect("value-changed", self.state_changed)
        self.fileButton.connect("clicked", self.filechanged)
        self.playbutton.connect("clicked", self.movieplaying)
        self.pausebutton.connect("clicked", self.movieplaying)

#
# Events-----------------------------------------------------
#

#used to reload table and states if config.ini is changed
    def startup(self, *args):
        importlib.reload(config)
        config.init(self.directory)
        config.path_root = os.path.dirname(self.directory)
        #changing path_states could cause errors
        config.path_states = "%s/states/" % config.path_root
        self.tempLabel.set_text(str(config.main_temperature))
        self.changeImage()
        self.energy_changed()
        self.processtable()


#display area
    def changeImage(self, *args):
        # without interpolation
        if self.interpolationCB.get_active() == False:
            states = glob.glob("%s*" %config.path_states)
            i = 0
            while ("%s%d" %(config.path_states,i)) in states:
                i+=1
            numStates = i-1
            self.stateScale.set_range(0, numStates)
            self.stateSB.set_range(0, numStates)
            processes = glob.glob("%s%d/procdata/*" %(config.path_states, self.stateScale.get_value()))
            j = 0
            while ("%s%d/procdata/saddle_%d.con" %(config.path_states,self.stateScale.get_value(), j)) in processes:
                j+=1
            numProcesses = j-1
            if numProcesses == 0:
                self.processesSB.set_value(0)
                self.processesSB.set_sensitive(False)
            else:
                self.processesSB.set_sensitive(True)
                self.processesSB.set_range(0, numProcesses)
            saddle = io.loadcon("%s%d/procdata/saddle_%d.con" %(config.path_states,self.stateScale.get_value(), self.processesSB.get_value()))
            reactant = io.loadcon("%s%d/procdata/reactant_%d.con" %(config.path_states,self.stateScale.get_value(), self.processesSB.get_value()))
            product = io.loadcon("%s%d/procdata/product_%d.con" %(config.path_states,self.stateScale.get_value(), self.processesSB.get_value()))
            datapass = [reactant, saddle, product]
        # with interpolation
        else:
            processes = glob.glob("%s%d/procdata/*" %(config.path_states, self.stateScale.get_value()))
            j = 0
            while ("%s%d/procdata/saddle_%d.con" %(config.path_states,self.stateScale.get_value(), j)) in processes:
                j+=1
            numProcesses = j-1
            if numProcesses == 0:
                self.processesSB.set_value(0)
                self.processesSB.set_sensitive(False)
            else:
                self.processesSB.set_sensitive(True)
                self.processesSB.set_range(0, numProcesses)
            saddle = io.loadcon("%s%d/procdata/saddle_%d.con" %(config.path_states,self.stateScale.get_value(), self.processesSB.get_value()))
            reactant = io.loadcon("%s%d/procdata/reactant_%d.con" %(config.path_states,self.stateScale.get_value(), self.processesSB.get_value()))
            product = io.loadcon("%s%d/procdata/product_%d.con" %(config.path_states,self.stateScale.get_value(), self.processesSB.get_value()))
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


    # creates process table
    def processtable(self, *args):
        #allows user to re-order table
        def sortable(a, b, c, d):
            x = self.store.get_value(b, d)
            y = self.store.get_value(c, d)
            if x>y:
                return 1
            if x==y:
                return 0
            if x<y:
                return -1

        # re creates table if a table is already there
        try:
            self.scrollview.destroy()
            self.view.destroy()
        except:
            pass
        self.store = gtk.ListStore(int, float, str, int, float, str, float, str, int)
        self.view = gtk.TreeView(model=self.store)
        self.store.set_sort_func(2, sortable, 2)
        self.store.set_sort_func(5, sortable, 5)
        self.store.set_sort_func(7, sortable, 7)

        self.proc_column = gtk.TreeViewColumn('proc #')
        self.proc_column.set_sort_column_id(0)
        self.saddleenergy_column = gtk.TreeViewColumn('saddle energy')
        self.saddleenergy_column.set_sort_column_id(1)
        self.prefactor_column = gtk.TreeViewColumn('prefactor')
        self.prefactor_column.set_sort_column_id(2)
        self.product_column = gtk.TreeViewColumn('product')
        self.product_column.set_sort_column_id(3)
        self.productenergy_column = gtk.TreeViewColumn('product energy')
        self.productenergy_column.set_sort_column_id(4)
        self.productprefactor_column = gtk.TreeViewColumn('product prefactor')
        self.productprefactor_column.set_sort_column_id(5)
        self.barrier_column = gtk.TreeViewColumn('barrier')
        self.barrier_column.set_sort_column_id(6)
        self.rate_column = gtk.TreeViewColumn('rate')
        self.rate_column.set_sort_column_id(7)
        self.repeat_column = gtk.TreeViewColumn('repeats')
        self.repeat_column.set_sort_column_id(8)

        self.view.append_column(self.proc_column)
        self.view.append_column(self.saddleenergy_column)
        self.view.append_column(self.prefactor_column)
        self.view.append_column(self.product_column)
        self.view.append_column(self.productenergy_column)
        self.view.append_column(self.productprefactor_column)
        self.view.append_column(self.barrier_column)
        self.view.append_column(self.rate_column)
        self.view.append_column(self.repeat_column)

        proc_cell = gtk.CellRendererText()
        saddleenergy_cell = gtk.CellRendererText()
        prefactor_cell = gtk.CellRendererText()
        product_cell = gtk.CellRendererText()
        productenergy_cell = gtk.CellRendererText()
        productprefactor_cell = gtk.CellRendererText()
        barrier_cell = gtk.CellRendererText()
        rate_cell = gtk.CellRendererText()
        repeat_cell = gtk.CellRendererText()

        self.proc_column.pack_start(proc_cell)
        self.saddleenergy_column.pack_start(saddleenergy_cell)
        self.prefactor_column.pack_start(prefactor_cell)
        self.product_column.pack_start(product_cell)
        self.productenergy_column.pack_start(productenergy_cell)
        self.productprefactor_column.pack_start(productprefactor_cell)
        self.barrier_column.pack_start(barrier_cell)
        self.rate_column.pack_start(rate_cell)
        self.repeat_column.pack_start(repeat_cell)

        self.proc_column.add_attribute(proc_cell, 'text', 0)
        self.saddleenergy_column.add_attribute(saddleenergy_cell, 'text', 1)
        self.prefactor_column.add_attribute(prefactor_cell, 'text', 2)
        self.product_column.add_attribute(product_cell, 'text', 3)
        self.productenergy_column.add_attribute(productenergy_cell, 'text', 4)
        self.productprefactor_column.add_attribute(productprefactor_cell, 'text', 5)
        self.barrier_column.add_attribute(barrier_cell, 'text', 6)
        self.rate_column.add_attribute(rate_cell, 'text', 7)
        self.repeat_column.add_attribute(repeat_cell, 'text', 8)
        # populates table with values from processtable
        a = open("%s%s/processtable" %(config.path_states, int(self.stateScale.get_value())) , "r")
        lines = a.readlines()
        rows = {}

        for i in range(1,len(lines)):
            rows[i] = lines[i].split()
            rows[i] = [float(j) for j in rows[i]]
            self.store.append(rows[i])
        # connects, packs, and shows table
        self.view.connect("row-activated", self.rowactivated)
        self.scrollview = gtk.ScrolledWindow()
        self.scrollview.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        self.scrollview.set_size_request(-1, 150)
        self.scrollview.add(self.view)
        self.view.set_headers_clickable(True)
        self.hbox.pack_start(self.scrollview)
        self.scrollview.show_all()
        self.view.show_all()


    # lets user click on a table row to set process number
    def rowactivated(self, *args):
        selection = self.view.get_selection().get_selected()
        proc = self.store.get(selection[1], 0)
        self.processesSB.set_value(proc[0])


    # Dialog window to change directory or config file
    def filechanged(self, *args):
        dialog = gtk.FileChooserDialog(action=gtk.FILE_CHOOSER_ACTION_OPEN, buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))
        response = dialog.run()
        if response == gtk.RESPONSE_OK:
            self.fileEntry.set_text(dialog.get_filename())
            self.directory = dialog.get_filename()
            dialog.destroy()
        elif response == gtk.RESPONSE_CANCEL:
            dialog.destroy()
        self.startup()


    # allows use of interpolation Spin Button
    def interpolationCB_changed(self, widget, data=None):
        if self.interpolationCB.get_active() == True:
            self.interpolationSB.set_sensitive(True)
        else:
            self.interpolationSB.set_sensitive(False)


    # allows stateSlider to change with stateSB
    def stateSB_changed(self, *args):
        self.stateScale.set_value(self.stateSB.get_value())


    # allows stateSB to change with stateSlider
    def state_changed(self, *args):
        self.stateSB.set_value(self.stateScale.get_value())


    def movieplaying(self, *args):
        if self.playing == True:
            self.statePlayTB.set_sensitive(False)
        else:
            self.statePlayTB.set_sensitive(True)


    # allows clicking play button near statesScale to play through states
    def state_play(self, *args):
        try:
            gobject.source_remove(self.timer_id)
        except:
            pass
        #pause
        if self.statePlayTB.get_active() == False:
            self.statePlayTB.set_image(self.playImage)
            self.interpolationCB.set_sensitive(True)
            self.processesSB.set_sensitive(True)
            self.playbutton.set_sensitive(True)
            self.pausebutton.set_sensitive(True)
            self.changeImage()
        #play
        else:
            self.statePlayTB.set_image(self.pauseImage)
            self.interpolationCB.set_sensitive(False)
            self.processesSB.set_sensitive(False)
            self.playbutton.set_sensitive(False)
            self.pausebutton.set_sensitive(False)
            states = glob.glob("%s*" %config.path_states)
            i = 0
            while ("%s%d" %(config.path_states,i)) in states:
                i+=1
            numStates = i-1
            def loop():
                self.stateScale.set_value(int (self.stateScale.get_value())+1)
                if self.stateScale.get_value() >= numStates:
                    self.stateScale.set_value(0)
                return True
            self.timer_id = gobject.timeout_add(1000/int(self.fps.get_value()), loop)
            #self.timer_id = gobject.timeout_add(1000/(int (self.fps.get_value())), loop)


    # Displays current state's energy
    def energy_changed(self, *args):
        energy = open("%sstate_table" %config.path_states, "r")
        energyNumber = energy.readlines()
        a = energyNumber[int (self.stateScale.get_value())].split()[1]
        self.stateEnergy.set_markup("%.3f<b>eV</b>" % float (a))


    # creates plot graph when plot button is pressed
    def energy_plot(self, *args):
        a = open("%sstate_table" %config.path_states, "r")
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
        ax.plot(x, y, '-' 'o')
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
    print(os.getcwd())
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
