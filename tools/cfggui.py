#!/usr/bin/env python

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
            
    def info(self, widget, event, d=None): 
        self.infoWindow = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "%s" %d)
        self.infoWindow.show_all()
        result = self.infoWindow.run()
        self.infoWindow.hide()
        
    def saveCheck(self, widget, data=None):
        self.saveButton.set_sensitive(True)
        
    def expose_event(self, descLabel, allocation, data=None):
        sections = config.format .keys()
        for i in range(len(self.descLabels)):
            self.descLabels[sections[i]].set_size_request(allocation.width, -1)
            
    
            
    
    #saves changes to config.ini
    def save(self, widget, data=None):
        sections = config.format.keys()
        for i in range(len(sections)):
            keys = config.format[sections[i]]['keys'].keys()
            for j in range(len(keys)):
                #strings with values
                if len(config.format[sections[i]]['keys'][keys[j]]['values']) != 0:
                    try:
                        self.config.set('%s' %sections[i], '%s' %keys[j], '%s' %self.buttons[keys[j]].get_active_text())
                    except:
                        pass
                #string without values, ints, and floats
                if (config.format[sections[i]]['keys'][keys[j]]['kind'] == 'string' and len(config.format[sections[i]]['keys'][keys[j]]['values']) == 0) or config.format[sections[i]]['keys'][keys[j]]['kind'] == 'int' or config.format[sections[i]]['keys'][keys[j]]['kind'] == 'float':
                    try:
                        self.config.set('%s' %sections[i], '%s' %keys[j], '%s' %self.buttons[keys[j]].get_text())
                    except:
                        pass
                #booleans
                if config.format[sections[i]]['keys'][keys[j]]['kind'] == 'boolean':
                    try:
                        if self.buttons[keys[j]].get_active() == True:
                            self.config.set('%s' %sections[i], '%s' %keys[j], 'True')
                        else:
                            self.config.set('%s' %sections[i], '%s' %keys[j], 'False')
                    except:
                        pass
        f = open("config.ini", 'w')
        self.config.write(f)
        f.close()
            

    def __init__(self):
        
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
        self.table = gtk.Table(3,6,False)
        self.window.add(self.table)
        self.notebook = gtk.Notebook()
        self.table.attach(self.notebook, 0,6,1,2)
        self.notebook.set_tab_pos(gtk.POS_LEFT)
        self.VbuttonBox = gtk.VBox()
        tempLabel =gtk.Label(' ')
        self.buttonBox = gtk.HBox()
        self.VbuttonBox.pack_start(self.buttonBox, False, False)
        self.VbuttonBox.pack_start(tempLabel, False, False)
        self.table.attach(self.VbuttonBox, 0,6,2,3)
        
        #save button and close button
        closeButton = gtk.Button()
        closeButton.set_label("Close")
        closeButton.connect("clicked", self.delete_event)
        self.saveButton = gtk.Button()
        self.saveButton.set_label("Save")
        self.saveButton.connect("clicked", self.save)
        self.saveButton.set_sensitive(False)
        self.buttonBox.pack_start(self.saveButton, False, False)
        self.buttonBox.pack_end(closeButton, False, False)
        
        
        #sections
        self.buttons = {}
        self.descLabels = {}
        sections = config.format.keys()
        sections.remove("Main")
        sections.insert(0, "Main")
        for i in range(len(sections)):
            scrollWindow = gtk.ScrolledWindow()
            label = gtk.Label("%s" %sections[i])
            Htable = gtk.Table((len(config.format[sections[i]]['keys'].keys())),2, True)
            tabBox = gtk.VBox()
            tabBox.set_border_width(10)
            frame = gtk.Frame()
            descSeparator = gtk.HSeparator()
            viewport = gtk.Viewport()
            viewport.set_shadow_type(gtk.SHADOW_NONE)
            viewport.add(tabBox)
            scrollWindow.add(viewport)
            scrollWindow.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
            scrollWindow.set_shadow_type(gtk.SHADOW_NONE)
            self.notebook.append_page(scrollWindow, label)
            sectionLabel = gtk.Label()
            sectionLabel.set_markup("<b>%s: </b> " %sections[i])
            sectionLabel.set_alignment(0,0)
            self.descLabels[sections[i]] = gtk.Label("%s" %config.format[sections[i]]['description'])
            self.descLabels[sections[i]].set_line_wrap(True)
            self.descLabels[sections[i]].set_alignment(0,0)
            self.descLabels[sections[i]].set_justify(gtk.JUSTIFY_FILL) 
            self.descLabels[sections[i]].connect("size-allocate", self.expose_event)
            tabBox.pack_start(sectionLabel, False, False)
            tabBox.pack_start(self.descLabels[sections[i]], False, False)
            tabBox.pack_start(descSeparator, False, False, 6)
            tabBox.pack_start(Htable, False, False)
            #keys
            for j in range(len(config.format[sections[i]]['keys'].keys())):
                keys = config.format[sections[i]]['keys'].keys()
                infoButton = gtk.EventBox()
                infoImage = gtk.Image()
                infoImage.set_from_stock(gtk.STOCK_INFO, gtk.ICON_SIZE_BUTTON)
                infoButton.add(infoImage)
                infoButton.connect("button_press_event", self.info, config.format[sections[i]]['keys'][keys[j]]['description'])
                
                #strings with values
                if len(config.format[sections[i]]['keys'][keys[j]]['values']) != 0:
                    name = keys[j]
                    self.buttons[name] = gtk.combo_box_new_text()
                    values = config.format[sections[i]]['keys'][keys[j]]['values'].keys()
                    for k in range(len(values)):
                        self.buttons[name].append_text(values[k])
                        try:
                            if values[k] == self.config.get('%s' %sections[i], '%s' %keys[j]):
                                self.buttons[name].set_active(k)
                        except:
                            pass                            
                            
                    nameLabel = gtk.Label(" %s:   " %keys[j])
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.buttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding = 5)    
                           
                #strings without values
                if config.format[sections[i]]['keys'][keys[j]]['kind'] == 'string' and len(config.format[sections[i]]['keys'][keys[j]]['values']) == 0:
                    name = keys[j]
                    self.buttons[name] = gtk.Entry()
                    self.buttons[name].connect("changed", self.saveCheck)
                    try:
                        self.buttons[name].set_text(self.config.get('%s' %sections[i], '%s' %keys[j]))
                    except:
                        pass
                    nameLabel = gtk.Label(" %s:   " %keys[j])
                    hbox = gtk.HBox() 
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.buttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding=5)    
                    
                #ints & floats
                if config.format[sections[i]]['keys'][keys[j]]['kind'] == 'int' or config.format[sections[i]]['keys'][keys[j]]['kind'] == 'float' :
                    name = keys[j]
                    self.buttons[name] = gtk.Entry()
                    self.buttons[name].connect("changed", self.saveCheck)
                    try:
                        self.buttons[name].set_text(self.config.get('%s' %sections[i], '%s' %keys[j]))
                    except:
                        pass
                    nameLabel = gtk.Label(" %s:   " %keys[j])
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.buttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding=5)   
                    
                #booleans
                if config.format[sections[i]]['keys'][keys[j]]['kind'] == 'boolean':
                    name = keys[j]
                    self.buttons[name] = gtk.RadioButton(label="True")
                    self.buttons[name].connect("clicked", self.saveCheck)
                    RB2 = gtk.RadioButton(label ="False", group=self.buttons[name])
                    try:
                        if self.config.get('%s' %sections[i], '%s' %keys[j]) == 'True':
                            self.buttons[name].set_active(True)
                        else:
                            RB2.set_active(True)
                    except:
                        pass   
                    alignment = gtk.Alignment(yalign = .5)
                    nameLabel = gtk.Label(" %s:   " %keys[j])
                    RBhbox = gtk.HBox()
                    separator = gtk.SeparatorToolItem()
                    RBhbox.pack_start(self.buttons[name], False, False)
                    RBhbox.pack_start(separator)
                    RBhbox.pack_start(RB2, False, False)
                    alignment.add(RBhbox)
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(alignment,1,2,j,j+1, False | gtk.FILL, False,ypadding=5)
                          
        self.window.connect("delete_event", self.delete_event)
        self.window.show_all()
        
    
        
        
        
        
              
if __name__ == "__main__":
    gui = cfggui()
    gtk.main()
