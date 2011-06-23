#!/usr/bin/env python

import pathfix
import config
import ConfigParser
import gtk
import gtk.gdk as gdk
import pygtk
import os


class cfggui():

    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False
            
    def info(self, widget, event, d=None): 
        infoWindow = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "%s" %d)
        infoWindow.show_all()
        result = infoWindow.run()
        infoWindow.hide()  
            
    
            
    
    #saves changes to config.ini
    def save(self, widget, data=None):
        for i in range(len(config.format)):
            for j in range(len(config.format[i].keys)):
                #strings with values
                if len(config.format[i].keys[j].values) != 0:
                    try:
                        self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, '%s' %self.buttons[config.format[i].keys[j].name].get_active_text())
                    except:
                        pass
                #string without values, ints, and floats
                if (config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0) or config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float':
                    try:
                        self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, '%s' %self.buttons[config.format[i].keys[j].name].get_text())
                    except:
                        pass
                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    try:
                        if self.buttons[config.format[i].keys[j].name].get_active() == True:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, 'True')
                        else:
                            self.config.set('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name, 'False')
                    except:
                        pass
        self.saveLabel.set_text("Changes saved")
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
        self.window.set_title("Eon Config")           
        self.table = gtk.Table(3,6,False)
        self.window.add(self.table)
        self.notebook = gtk.Notebook()
        self.table.attach(self.notebook, 0,6,1,2)
        self.notebook.set_tab_pos(gtk.POS_LEFT)
        self.show_tabs = True
        self.show_boader = True
        
        #save button and text
        menu = gtk.HBox()
        self.table.attach(menu, 0,1,0,1)
        saveImage = gtk.Image()
        saveImage.set_from_stock(gtk.STOCK_FLOPPY, gtk.ICON_SIZE_LARGE_TOOLBAR)
        self.saveLabel = gtk.Label("Please save Changes   ")
        saveButton = gtk.EventBox()
        saveButton.add(saveImage)
        menu.pack_start(self.saveLabel, False, False)
        menu.pack_start(saveButton, False, False)
        saveButton.connect("button_press_event", self.save)
        
        
        
        self.buttons = {}
        #sections
        for i in range(len(config.format)):
            label = gtk.Label("%s" %config.format[i].name)
            vbox = gtk.VBox()
            self.notebook.append_page(vbox, label)
            descLabel = gtk.Label("%s" %config.format[i].description)
            descLabel.set_line_wrap(True)
            descLabel.set_alignment(0,0)
            vbox.pack_start(descLabel, False, True) 
            #keys
            for j in range(len(config.format[i].keys)):
                infoButton = gtk.EventBox()
                infoImage = gtk.Image()
                infoImage.set_from_stock(gtk.STOCK_INFO, gtk.ICON_SIZE_BUTTON)
                infoButton.add(infoImage)
                infoButton.connect("button_press_event", self.info, config.format[i].keys[j].description)
                
                #strings with values
                if len(config.format[i].keys[j].values) != 0:
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.combo_box_new_text()
                    for k in range(len(config.format[i].keys[j].values)):
                        self.buttons[name].append_text(config.format[i].keys[j].values[k].name)
                        try:
                            if config.format[i].keys[j].values[k].name == self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name):
                                self.buttons[name].set_active(k)
                        except:
                            pass                            
                        
                        
                    frame = gtk.Frame()    
                    alignment = gtk.Alignment(yalign = .5)
                    alignment.add(self.buttons[name])
                    nameLabel = gtk.Label("%s:   " %config.format[i].keys[j].name)    
                    hbox = gtk.HBox()
                    vbox.pack_start(frame)
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    frame.add(hbox)
                    hbox.pack_start(alignment, False, False)
                           
                #strings without values
                if config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0 :
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.Entry()
                    try:
                        self.buttons[name].set_text(self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name))
                    except:
                        pass
                    frame = gtk.Frame()
                    nameLabel = gtk.Label("%s:   " %config.format[i].keys[j].name)    
                    hbox = gtk.HBox()
                    vbox.pack_start(frame)
                    frame.add(hbox)
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    hbox.pack_start(self.buttons[name], False, False)
                    
                #ints & floats
                if config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float':
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.Entry()
                    try:
                        self.buttons[name].set_text(self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name))
                    except:
                        pass
                    
                    frame = gtk.Frame()
                    nameLabel = gtk.Label("%s:   " %config.format[i].keys[j].name)    
                    hbox = gtk.HBox()
                    frame.add(hbox)
                    vbox.pack_start(frame)
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    hbox.pack_start(self.buttons[name], False, False)
                    
                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.RadioButton(label="True")
                    RB2 = gtk.RadioButton(label ="False", group=self.buttons[name])
                    try:
                        if self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name) == 'True':
                            self.buttons[name].set_active(True)
                        else:
                            RB2.set_active(True)
                    except:
                        pass   
                    
                    frame = gtk.Frame()
                    alignment = gtk.Alignment(yalign = .5)
                    nameLabel = gtk.Label("%s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox()
                    frame.add(hbox)
                    vbox.pack_start(frame)
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    RBvbox = gtk.VBox()
                    RBvbox.pack_start(self.buttons[name], False, False)
                    RBvbox.pack_start(RB2, False, False)
                    alignment.add(RBvbox)
                    hbox.pack_start(alignment, False, False)
                    
              
        self.window.connect("delete_event", self.delete_event)
        self.window.show_all()
        
    
        
        
        
        
              
if __name__ == "__main__":
    gui = cfggui()
    gtk.main()
