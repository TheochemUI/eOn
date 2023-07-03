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

# Authors: Ian Johnson

import pathfix
import config
import configparser
import gtk
import gtk.gdk as gdk
import pygtk
import os


class cfggui():


    def delete_event(self, widget, data=None):
        gtk.main_quit()
        return False


    # bolds options that will be saved
    def defaultchanged(self, widget, list=None):
        name = list[0]
        option = name.split(', ')[1]
        i = list[1]
        j = list[2]

        if config.format[i].keys[j].kind == "string" and len(config.format[i].keys[j].values) > 0:
            if self.buttons[name].get_active_text() != config.format[i].keys[j].default:
                self.nameLabels[name].set_markup("<b>%s:</b>" %option)
            else:
                self.nameLabels[name].set_markup("%s:" %option)
        if (config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0) or config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float':
            if self.buttons[name].get_text() != str (config.format[i].keys[j].default):
                self.nameLabels[name].set_markup("<b>%s:</b>" %option)
            else:
                self.nameLabels[name].set_markup("%s:" %option)
        if config.format[i].keys[j].kind == 'boolean':
            if config.format[i].keys[j].default == True and self.buttons[name].get_active() == False:
                self.nameLabels[name].set_markup("<b>%s:</b>" %option)
            elif config.format[i].keys[j].default == False and self.buttons[name].get_active() == True:
                self.nameLabels[name].set_markup("<b>%s:</b>" %option)
            else:
                self.nameLabels[name].set_markup("%s:" %option)


    # activates save button
    def saveCheck(self, widget, data=None):
        self.saveButton.set_sensitive(True)
        self.closeButton.set_label("Cancel")


    # allows only changed RB options to save
    def buttonChanged(self, widget, name=None):
        split = name.split(', ')
        if name not in self.changedbuttons:
            self.changedbuttons.append(name)
        if split[0] not in self.config.sections():
            self.config.add_section(split[0])


    # changes options to config default
    def refreshoption(self, widget, event, refreshlist=None):
        name = refreshlist[0]
        key = refreshlist[1]

        if key.kind == "int" or key.kind == "float" or key.kind == "string" and len(key.values) == 0:
            self.buttons[name].set_text(str (key.default))
        if key.kind == "string" and len(key.values) != 0:
            counter = 0
            for i in key.values:
                if i.name == key.default:
                    self.buttons[name].set_active(counter)
                counter +=1
        if key.kind == "boolean":
            if key.default == True:
                self.buttons[name].set_active(True)
            else:
                self.RB2[name].set_active(True)


#saves changes to config.ini
    def save(self, widget, data=None):
        for i in range(len(config.format)):
            for j in range(len(config.format[i].keys)):
                #strings with values
                if len(config.format[i].keys[j].values) != 0:
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    if name in self.changedbuttons:
                        self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, '%s' %self.buttons[name].get_active_text())

                #string without values, ints, and floats
                if (config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0) or config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float':
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    if name in self.changedbuttons:
                        self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, '%s' %self.buttons[name].get_text())
                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    if name in self.changedbuttons:
                        if self.buttons[name].get_active() == True:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, 'True')
                        else:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, 'False')

        f = open("config.ini", 'w')
        self.config.write(f)
        f.close()
        self.saveButton.set_sensitive(False)
        self.closeButton.set_label("Close")


#display
    def __init__(self):
        button_width = 100
        self.config = configparser.SafeConfigParser()
        try:
            self.config.read("./config.ini")
        except:
            print("No config.ini found in local directory, using default values.")
        #window, table, and notebook
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_border_width(10)
        self.window.set_geometry_hints(max_width = self.window.allocation.width)
        self.window.set_title("Eon Config")
        self.table = gtk.Table(2,6,False)
        self.table.set_row_spacings(10)
        self.notebook = gtk.Notebook()
        self.notebook.set_tab_pos(gtk.POS_LEFT)
        self.HbuttonBox = gtk.HBox()
        self.buttonBox = gtk.Table(1,2, True)
        self.buttonBox.set_col_spacings(10)
        self.HbuttonBox.pack_end(self.buttonBox, False, False)
        self.window.add(self.table)
        self.table.attach(self.HbuttonBox, 0,6,1,2)
        self.table.attach(self.notebook, 0,6,0,1)
        #save button and close button
        self.closeButton = gtk.Button()
        self.closeButton.set_label("Close")
        self.closeButton.connect("clicked", self.delete_event)
        self.saveButton = gtk.Button()
        self.saveButton.set_label("Save")
        self.saveButton.connect("clicked", self.save)
        saveButtonTT = gtk.Tooltips()
        saveButtonTT.set_tip(self.saveButton, "saves changed options to current directory's config.ini")
        self.saveButton.set_size_request(button_width, -1)
        self.closeButton.set_size_request(button_width, -1)
        self.buttonBox.attach(self.saveButton,0,1,0,1, False | gtk.FILL, False)
        self.buttonBox.attach(self.closeButton,1,2,0,1, False | gtk.FILL, False)
        #dictionaries & lists
        self.buttons = {}
        self.descLabels = {}
        self.changedbuttons = []
        self.nameLabels = {}
        self.RB2 = {}


#Sections
        for i in range(len(config.format)):
            scrollWindow = gtk.ScrolledWindow()
            label = gtk.Label("%s" %config.format[i].name)
            Htable = gtk.Table((len(config.format[i].keys)),2, True)
            tabBox = gtk.VBox()
            tabBox.set_border_width(10)
            frame = gtk.Frame()
            viewport = gtk.Viewport()
            viewport.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("#F7F6F6"))
            viewport.set_shadow_type(gtk.SHADOW_NONE)
            viewport.add(tabBox)
            scrollWindow.add(viewport)
            scrollWindow.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
            scrollWindow.set_shadow_type(gtk.SHADOW_NONE)
            self.notebook.append_page(scrollWindow, label)
            tabBox.pack_start(Htable, False, False)


#keys
            for j in range(len(config.format[i].keys)):
                refreshButton = gtk.EventBox()
                refreshImage = gtk.Image()
                tooltip2 = gtk.Tooltips()
                tooltip2.set_tip(refreshImage, "resets value to default from config.py")
                refreshImage.set_from_stock(gtk.STOCK_REFRESH, gtk.ICON_SIZE_BUTTON)
                refreshButton.add(refreshImage)


                #strings with values
                if len(config.format[i].keys[j].values) != 0:
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    self.buttons[name] = gtk.combo_box_new_text()
                    self.nameLabels[name] = gtk.Label("%s:" %config.format[i].keys[j].name)
                    for k in range(len(config.format[i].keys[j].values)):
                        self.buttons[name].append_text(str (config.format[i].keys[j].values[k]))
                        if config.format[i].keys[j].values[k] == config.format[i].keys[j].default:
                            default = k
                    try:
                        for k in range(len(config.format[i].keys[j].values)):
                            if config.format[i].keys[j].values[k].name == self.config.get(config.format[i].name,config.format[i].keys[j].name):
                                default = k
                        self.changedbuttons.append(name)
                    except:
                        pass
                    self.buttons[name].connect("changed", self.defaultchanged, [name,i,j])
                    try:
                        self.buttons[name].set_active(default)
                    except:
                        pass
                    hbox = gtk.HBox(spacing = 5)
                    hbox.pack_start(self.nameLabels[name], False, False)
                    buttonbox = gtk.HBox(spacing =5)
                    buttonbox.pack_start(self.buttons[name],True,True)
                    buttonbox.pack_start(refreshButton, False, True)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(buttonbox,1,2,j,j+1, False | gtk.FILL, False, ypadding = 5)
                    self.buttons[name].connect("changed", self.buttonChanged, name)
                    self.buttons[name].connect("changed", self.saveCheck)


                #strings without values
                if config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0:
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    self.nameLabels[name] = gtk.Label("%s:" %config.format[i].keys[j].name)
                    self.buttons[name] = gtk.Entry()
                    try:
                        self.buttons[name].set_text(config.format[i].keys[j].default)
                    except:
                        pass
                    self.buttons[name].connect("changed", self.defaultchanged, [name,i,j])
                    try:
                        self.buttons[name].set_text(self.config.get(config.format[i].name,config.format[i].keys[j].name))
                        self.changedbuttons.append(name)
                    except:
                        pass
                    hbox = gtk.HBox(spacing = 5)
                    hbox.pack_start(self.nameLabels[name], False, False)
                    buttonbox = gtk.HBox(spacing =5)
                    buttonbox.pack_start(self.buttons[name],True,True)
                    buttonbox.pack_start(refreshButton, False, True)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(buttonbox,1,2,j,j+1, False | gtk.FILL, False, ypadding=5)
                    self.buttons[name].connect("changed", self.buttonChanged, name)
                    self.buttons[name].connect("changed", self.saveCheck)


                #ints & floats
                if config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float' :
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    self.nameLabels[name] = gtk.Label("%s:" %config.format[i].keys[j].name)
                    self.buttons[name] = gtk.Entry()
                    try:
                        self.buttons[name].set_text(str (config.format[i].keys[j].default))
                    except:
                        pass
                    self.buttons[name].connect("changed", self.defaultchanged, [name,i,j])
                    try:
                        self.buttons[name].set_text(str (self.config.get(config.format[i].name,config.format[i].keys[j].name)))
                        self.changedbuttons.append(name)
                    except:
                        pass
                    hbox = gtk.HBox(spacing =5)
                    hbox.pack_start(self.nameLabels[name], False, False)
                    buttonbox = gtk.HBox(spacing =5)
                    buttonbox.pack_start(self.buttons[name],True,True)
                    buttonbox.pack_start(refreshButton, False, True)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(buttonbox,1,2,j,j+1, False | gtk.FILL, False, ypadding=5)
                    self.buttons[name].connect("changed", self.buttonChanged, name)
                    self.buttons[name].connect("changed", self.saveCheck)
                    self.buttons[name].connect("changed", self.buttonChanged, name)


                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    name = config.format[i].name + ", " + config.format[i].keys[j].name
                    self.nameLabels[name] = gtk.Label("%s:" %config.format[i].keys[j].name)
                    self.buttons[name] = gtk.RadioButton(label="True")
                    self.RB2[name] = gtk.RadioButton(label="False", group = self.buttons[name])
                    self.buttons[name].connect("clicked", self.saveCheck)
                    try:
                        if config.format[i].keys[j].default == True:
                            self.buttons[name].set_active(True)
                        if config.format[i].keys[j].default == False:
                            self.RB2[name].set_active(True)
                    except:
                        pass
                    self.buttons[name].connect("toggled", self.defaultchanged, [name,i,j])
                    self.RB2[name].connect("toggled", self.defaultchanged, [name,i,j])
                    try:
                        if self.config.get(config.format[i].name,config.format[i].keys[j].name) == 'True':
                            self.buttons[name].set_active(True)
                        if self.config.get(config.format[i].name,config.format[i].keys[j].name) == 'False':
                            self.RB2[name].set_active(True)
                        self.changedbuttons.append(name)
                    except:
                        pass
                    alignment = gtk.Alignment(yalign = .5)
                    RBhbox = gtk.HBox()
                    separator = gtk.SeparatorToolItem()
                    RBhbox.pack_start(self.buttons[name], False, False)
                    RBhbox.pack_start(separator)
                    RBhbox.pack_start(self.RB2[name], False, False)
                    alignment.add(RBhbox)
                    boolhbox = gtk.HBox()
                    boolhbox.pack_start(alignment, False, False)
                    boolhbox.pack_end(refreshButton, False, True)
                    hbox = gtk.HBox(spacing = 5)
                    hbox.pack_start(self.nameLabels[name], False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(boolhbox,1,2,j,j+1, False| gtk.FILL, False,ypadding=5)
                    self.buttons[name].connect("toggled", self.buttonChanged, name)
                    self.RB2[name].connect("toggled", self.buttonChanged, name)

                refreshlist = [name, config.format[i].keys[j]]
                refreshButton.connect("button_press_event", self.refreshoption, refreshlist)

        self.saveButton.set_sensitive(False)
        self.closeButton.set_label("Close")
        self.window.connect("delete_event", self.delete_event)
        self.window.show_all()


if __name__ == "__main__":
    gui = cfggui()
    gtk.main()
