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

class eongui(atomview.atomview):

    def __init__(self):
        gladetree = gtk.glade.XML("eon-gui.glade")
        gui = gladetree.get_widget("eongui")
        atomview.atomview.__init__(self, gui)

if __name__ == "__main__":
    pid = os.fork()
    if pid:
        os._exit(0)
    import io
    import sys
    q = eongui()
    if len(sys.argv) > 1:
        data = io.loadposcars(sys.argv[1])
        if len(data) < 1:
            data = io.loadcons(sys.argv[1])
        q.data_set(data)
    gtk.main()
