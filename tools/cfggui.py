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

#-------------------------------------- Instructions ------------------------------------------------------
#   To add options to the window edit the eon/config.py file.
#
# adding sections: fadd("your_section_name", description = "your_description")
# adding keys:     fadd("your_section_name", "your_key_name", kind ="key_type(string, int, ect)", description = "your_key_description")
# adding values:   fadd("your_section_name", "your_key_name", "your_value", description = "your_value_description")


import pathfix
import config
import ConfigParser
import gtk
import gtk.gdk as gdk
import pygtk
import os


class cfggui():  

    def delete_event(self, widget, data=None):
        gtk.main_quit()
        return False
            
    #creates info windows        
    def info(self, widget, event, d=None): 
        self.infoWindow = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "%s" %d)
        self.infoWindow.show_all()
        result = self.infoWindow.run()
        self.infoWindow.hide()
        
    #activates save button    
    def saveCheck(self, widget, data=None):
        self.saveButton.set_sensitive(True)
        
    #not implemented ignore for now    
    def expose_event(self, descLabel, allocation, data=None):
        for i in range(len(self.descLabels)):
            self.descLabels[config.format[i].name].set_size_request(allocation.width, -1)
            
    
#saves changes to config.ini
    def save(self, widget, data=None):
        for i in range(len(config.format)):
            for j in range(len(config.format[i].keys)):
                #strings with values
                if len(config.format[i].keys[j].values) != 0:
                    name = config.format[i].name + "," + config.format[i].keys[j].name
                    try:
                        self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, '%s' %self.CBbuttons[name].get_active_text())
                    except:
                        pass
                #string without values, ints, and floats
                if (config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0) or config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float':
                    name = config.format[i].name + "," + config.format[i].keys[j].name 
                    if self.TEbuttons[name].get_text() != '':
                        try:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, '%s' %self.TEbuttons[name].get_text())     
                        except:
                            pass
                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    name = config.format[i].name + "," + config.format[i].keys[j].name
                    try:
                        if self.RBbuttons[name].get_active() == True:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, 'True')
                        else:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, 'False')
                    except:
                        pass
        f = open("config.ini", 'w')
        self.config.write(f)
        f.close()
        self.saveButton.set_sensitive(False)
            

    def __init__(self):
#display
        button_width = 100
        self.config = ConfigParser.SafeConfigParser()
        self.config.read(os.path.join(pathfix.path, "default_config.ini"))
        try:
            self.config.read("./config.ini")
        except:
            print "No config.ini found in local directory, using default values."      
        #window, table, and notebook      
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.set_border_width(10)
        self.window.set_geometry_hints(max_width = self.window.allocation.width)
        self.window.set_title("Eon Config")           
        self.table = gtk.Table(2,6,False)
        self.table.set_row_spacings(10)
        self.window.add(self.table)
        self.notebook = gtk.Notebook()
        self.table.attach(self.notebook, 0,6,0,1)
        self.notebook.set_tab_pos(gtk.POS_LEFT)
        self.HbuttonBox = gtk.HBox()
        self.buttonBox = gtk.Table(1,2, True)
        self.buttonBox.set_col_spacings(10)
        self.HbuttonBox.pack_end(self.buttonBox, False, False)
        self.table.attach(self.HbuttonBox, 0,6,1,2)
        #save button and close button
        closeButton = gtk.Button()
        closeButton.set_label("Cancel")
        closeButton.connect("clicked", self.delete_event)
        self.saveButton = gtk.Button()
        self.saveButton.set_label("Save")
        self.saveButton.connect("clicked", self.save)
        self.saveButton.set_sensitive(False)
        self.saveButton.set_size_request(button_width, -1)
        closeButton.set_size_request(button_width, -1)
        self.buttonBox.attach(self.saveButton,0,1,0,1, False | gtk.FILL, False)
        self.buttonBox.attach(closeButton,1,2,0,1, False | gtk.FILL, False)
        #dictionaries
        self.CBbuttons = {}
        self.TEbuttons = {}
        self.RBbuttons = {}
        self.descLabels = {}
        
        
#Sections
        for i in range(len(config.format)):
            scrollWindow = gtk.ScrolledWindow()
            label = gtk.Label("%s" %config.format[i].name)
            Htable = gtk.Table((len(config.format[i].keys)),2, True)
            tabBox = gtk.VBox()
            tabBox.set_border_width(10)
            frame = gtk.Frame()
            descSeparator = gtk.HSeparator()
            viewport = gtk.Viewport()
            viewport.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse("#F7F6F6"))
            viewport.set_shadow_type(gtk.SHADOW_NONE)
            viewport.add(tabBox)
            scrollWindow.add(viewport)
            scrollWindow.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
            scrollWindow.set_shadow_type(gtk.SHADOW_NONE)
            self.notebook.append_page(scrollWindow, label)
            sectionLabel = gtk.Label()
            sectionLabel.set_markup("<b>%s: </b> " %config.format[i].name)
            sectionLabel.set_alignment(0,0)
            descname = config.format[i].name
            self.descLabels[descname] = gtk.Label("%s" %config.format[i].description)
            self.descLabels[descname].set_line_wrap(True)
            self.descLabels[descname].set_alignment(0,0)
            self.descLabels[descname].set_justify(gtk.JUSTIFY_FILL) 
            self.descLabels[descname].connect("size-allocate", self.expose_event)
            tabBox.pack_start(sectionLabel, False, False)
            tabBox.pack_start(self.descLabels[descname], False, False)
            tabBox.pack_start(descSeparator, False, False, 6)
            tabBox.pack_start(Htable, False, False)
#keys
            for j in range(len(config.format[i].keys)):
                infoButton = gtk.EventBox()
                infoImage = gtk.Image()
                infoImage.set_from_stock(gtk.STOCK_INFO, gtk.ICON_SIZE_BUTTON)
                infoButton.add(infoImage)
                infoButton.connect("button_press_event", self.info, config.format[i].keys[j].description)
                
                #strings with values
                if len(config.format[i].keys[j].values) != 0:
                    name = config.format[i].name + "," + config.format[i].keys[j].name
                    self.CBbuttons[name] = gtk.combo_box_new_text()
                    for k in range(len(config.format[i].keys[j].values)):
                        self.CBbuttons[name].append_text(str (config.format[i].keys[j].values[k].name))
                        try:
                            if config.format[i].keys[j].values[k].name == self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name):
                                self.CBbuttons[name].set_active(k)
                        except:
                            pass                                   
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.CBbuttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding = 5) 
                    self.CBbuttons[name].connect("changed", self.saveCheck)   
                           
                #strings without values
                if config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0:
                    name = config.format[i].name + "," + config.format[i].keys[j].name
                    self.TEbuttons[name] = gtk.Entry()
                    self.TEbuttons[name].connect("changed", self.saveCheck)
                    try:
                        self.TEbuttons[name].set_text(self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name))
                    except:
                        pass
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox() 
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.TEbuttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding=5)    
                    
                #ints & floats
                if config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float' :
                    name = config.format[i].name + "," + config.format[i].keys[j].name
                    self.TEbuttons[name] = gtk.Entry()
                    self.TEbuttons[name].connect("changed", self.saveCheck)
                    try:
                        self.TEbuttons[name].set_text(self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name))
                    except:
                        pass
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.TEbuttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding=5)   
                    
                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    name = config.format[i].name + "," + config.format[i].keys[j].name
                    self.RBbuttons[name] = gtk.RadioButton(label="True")
                    self.RBbuttons[name].connect("clicked", self.saveCheck)
                    RB2 = gtk.RadioButton(label ="False", group=self.RBbuttons[name])
                    try:
                        if self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name) == 'True':
                            self.RBbuttons[name].set_active(True)
                        else:
                            RB2.set_active(True)
                    except:
                        pass   
                    alignment = gtk.Alignment(yalign = .5)
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    RBhbox = gtk.HBox()
                    separator = gtk.SeparatorToolItem()
                    RBhbox.pack_start(self.RBbuttons[name], False, False)
                    RBhbox.pack_start(separator)
                    RBhbox.pack_start(RB2, False, False)
                    alignment.add(RBhbox)
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(alignment,1,2,j,j+1, False | gtk.FILL, False,ypadding=5)
        self.saveButton.set_sensitive(False)                  
        self.window.connect("delete_event", self.delete_event)
        self.window.show_all()
        
                  
if __name__ == "__main__":
    gui = cfggui()
    gtk.main()
