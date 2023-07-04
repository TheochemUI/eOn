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
# Author: Ian Johnson

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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as drawArea
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
import importlib

class BHgui(atomview.atomview):
    def __init__(self):
        # Glade imports
        gladetree = gtk.glade.XML(os.path.join(pathfix.path, "tools/bh-gui.glade"))
        gui = gladetree.get_widget("BHgui")
        atomview.atomview.__init__(self,gui)
        self.stateSB = gladetree.get_widget("stateSB")
        self.stateScale = gladetree.get_widget("stateScale")
        self.scrollwindow = gladetree.get_widget("scrollwindow")
        self.energyplotButton = gladetree.get_widget("energyplotButton")
        self.repeatplotButton = gladetree.get_widget("repeatplotButton")
        self.eplotWindow = gladetree.get_widget("eplotWindow")
        self.rplotWindow = gladetree.get_widget("rplotWindow")
        self.acceptanceratio = gladetree.get_widget("acceptanceratio")
        # Image setting
        self.eplotImage = gtk.Image()
        image = gdk.pixbuf_new_from_file_at_size(os.path.join(pathfix.path, "tools/plot_icon.png"), 15, 15)
        self.eplotImage.set_from_pixbuf(image)
        self.energyplotButton.set_image(self.eplotImage)
        self.rplotImage = gtk.Image()
        self.rplotImage.set_from_pixbuf(image)
        self.repeatplotButton.set_image(self.rplotImage)
        # Config.ini imports
        self.directory = "./config.ini"
        importlib.reload(config)
        config.init(self.directory)
        config.path_root = os.path.dirname(self.directory)
        config.path_states = "%s/states/" % config.path_root
        # Signal connectors
        self.stateScale.connect("value-changed", self.changeImage)
        self.stateSB.connect("value-changed", self.stateSB_changed)
        self.stateScale.connect("value-changed", self.state_changed)
        self.energyplotButton.connect("clicked", self.energyplot)
        self.repeatplotButton.connect("clicked", self.repeatplot)

        self.startup()

    # startup
    def startup(self, *args):
        self.changeImage()
        self.statetable()
        self.acceptanceratio_average()


    # minimum.con display
    def changeImage(self, *args):
        states = glob.glob("%s*" %config.path_states)
        i = 0
        while ("%s%d" %(config.path_states,i)) in states:
            i+=1
        self.numstates = i-1
        self.stateScale.set_range(0,self.numstates)
        self.stateSB.set_range(0,self.numstates)
        datapass = io.loadcon("%s%d/minimum.con" %(config.path_states,self.stateScale.get_value()))
        self.data_set([datapass])


    # changes state by clicking table rows
    def rowactivated(self, *args):
        selection = self.view.get_selection().get_selected()
        proc = self.store.get(selection[1], 0)
        self.stateScale.set_value(proc[0])


    # table
    def statetable(self, *args):
        self.store = gtk.ListStore(int, float, int)
        self.view = gtk.TreeView(model=self.store)
        self.view.connect("row-activated", self.rowactivated)

        self.statecolumn = gtk.TreeViewColumn('state')
        self.statecolumn.set_sort_column_id(0)
        self.energycolumn = gtk.TreeViewColumn('energy')
        self.energycolumn.set_sort_column_id(1)
        self.repeatcolumn = gtk.TreeViewColumn('repeats')
        self.repeatcolumn.set_sort_column_id(2)

        self.view.append_column(self.statecolumn)
        self.view.append_column(self.energycolumn)
        self.view.append_column(self.repeatcolumn)

        statecell = gtk.CellRendererText()
        energycell = gtk.CellRendererText()
        repeatcell = gtk.CellRendererText()

        self.statecolumn.pack_start(statecell)
        self.energycolumn.pack_start(energycell)
        self.repeatcolumn.pack_start(repeatcell)

        self.statecolumn.add_attribute(statecell,'text',0)
        self.energycolumn.add_attribute(energycell,'text',1)
        self.repeatcolumn.add_attribute(repeatcell,'text',2)

        table = open('%s/state_table' %config.path_states,'r')
        lines = table.readlines()
        rows = {}

        for i in range(2,len(lines)):
            rows[i] = lines[i].split()
            rows[i] = [float(j) for j in rows[i]]
            self.store.append(rows[i])

        self.scrollwindow.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
        self.scrollwindow.set_size_request(-1, 100)
        self.scrollwindow.add(self.view)
        self.view.set_headers_clickable(True)
        self.scrollwindow.show_all()
        self.view.show_all()


    # keeps spinbutton and scale the same
    def stateSB_changed(self, *args):
        self.stateScale.set_value(self.stateSB.get_value())


    # keeps scale and spinbutton the same
    def state_changed(self, *args):
        self.stateSB.set_value(self.stateScale.get_value())


    # energy vs state plot
    def energyplot(self, *args):
        table = open("%sstate_table" %config.path_states, "r")
        lines = table.readlines()
        x = []
        y = []
        grouplabels = []
        for i in range(2,len(lines)):
            x.append(i-1)
            y.append(float (lines[i].split()[1]))
            grouplabels.append(lines[i].split()[0])
        y.reverse()
        grouplabels.reverse()
        self.eplotWindow.set_default_size(425,425)
        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y, 'bo' '-')
        p.xlabel("State")
        p.ylabel("eV")
        ax.set_xticks(x)
        ax.set_xticklabels(grouplabels)
        container = gtk.VBox()
        self.eplotWindow.add(container)
        graph = drawArea(fig)
        toolbar = NavigationToolbar(graph, self.eplotWindow)
        container.pack_start(graph)
        container.pack_start(toolbar, False, False)
        self.eplotWindow.show_all()


    # repeat vs state plot
    def repeatplot(self, *args):
        table = open("%sstate_table" %config.path_states, "r")
        lines = table.readlines()
        x = []
        y = []
        unsortedlines = []
        sortedlines = []
        grouplabels = []
        for i in range(2,len(lines)):
            x.append(i-2)
            unsortedlines.append((lines[i].split()[0],lines[i].split()[2]))
        sortedlines = sorted(unsortedlines, key  = lambda tuble: int (tuble[1]))
        for i in sortedlines:
            y.append(i[1])
            grouplabels.append(i[0])
        fig = p.figure()
        ax = fig.add_subplot(111)
        p.plot(x,y, 'bo' '-')
        p.ylabel("number of repeats")
        p.xlabel("states")
        lnr = float (max(y))
        fp = lnr*.05
        p.ylim([0-fp,lnr+fp])
        p.xlim([0-fp,max(x)+fp])
        ax.set_xticks(x)
        ax.set_xticklabels(grouplabels)
        self.rplotWindow.set_default_size(425,425)
        container = gtk.VBox()
        self.rplotWindow.add(container)
        graph = drawArea(fig)
        toolbar = NavigationToolbar(graph, self.rplotWindow)
        container.pack_start(graph)
        container.pack_start(toolbar, False, False)
        self.rplotWindow.show_all()


    # average acceptanceratio
    def acceptanceratio_average(self, *args):
        ratios = []
        for i in range(self.numstates+1):
            results = open("%s%d/results.dat" %(config.path_states, i), "r")
            lines = results.readlines()
            for j in lines:
                if "acceptance_ratio" in j:
                    ratios.append(float (j.split()[0]))
        if len(ratios) > 0:
            average = sum(ratios) / len(ratios)
        else:
            average = 0
        self.acceptanceratio.set_text(str (average))


if __name__ == "__main__":
    print(os.getcwd())
    pid = os.fork()
    if pid:
        os._exit(0)
    import io
    import sys
    q = BHgui()
    if len(sys.argv) > 1:
        data = io.loadposcars(sys.argv[1])
        if len(data) < 1:
            data = io.loadcons(sys.argv[1])
        q.data_set(data)
    gtk.main()
