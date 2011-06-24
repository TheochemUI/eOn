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
        
    def saveCheck(self, widget, data=None):
        self.saveLabel.set_text("Please save changes")  
            
    
            
    
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
        self.window.set_default_size(700, 1000)
        self.window.set_title("Eon Config")           
        self.table = gtk.Table(3,6,False)
        self.window.add(self.table)
        self.notebook = gtk.Notebook()
        self.table.attach(self.notebook, 0,6,1,2)
        self.notebook.set_tab_pos(gtk.POS_LEFT)
        self.show_tabs = True
        self.show_boader = True
        
        #save button and save text
        menu = gtk.HBox()
        self.table.attach(menu, 0,1,0,1, False, False)
        saveImage = gtk.Image()
        saveImage.set_from_stock(gtk.STOCK_FLOPPY, gtk.ICON_SIZE_LARGE_TOOLBAR)
        self.saveLabel = gtk.Label("Please save changes   ")
        saveButton = gtk.EventBox()
        saveButton.add(saveImage)
        menu.pack_start(self.saveLabel, False, False)
        menu.pack_start(saveButton, False, False)
        saveButton.connect("button_press_event", self.save)
        
        
        #sections
        self.buttons = {}
        for i in range(len(config.format)):
            scrollWindow = gtk.ScrolledWindow()
            label = gtk.Label("%s" %config.format[i].name)
            Htable = gtk.Table((len(config.format[i].keys)),2, True)
            tabBox = gtk.VBox()
            descSeparator = gtk.HSeparator()
            scrollWindow.add_with_viewport(tabBox)
            self.notebook.append_page(scrollWindow, label)
            sectionLabel = gtk.Label()
            sectionLabel.set_markup("<b> %s </b> " %config.format[i].name)
            descLabel = gtk.Label("%s" %config.format[i].description)
            descLabel.set_line_wrap(True)
            descLabel.set_alignment(0,0) 
            tabBox.pack_start(sectionLabel, False, False)
            tabBox.pack_start(descLabel, False, False)
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
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.combo_box_new_text()
                    self.buttons[name].connect("changed", self.saveCheck)
                    for k in range(len(config.format[i].keys[j].values)):
                        self.buttons[name].append_text(config.format[i].keys[j].values[k].name)
                        try:
                            if config.format[i].keys[j].values[k].name == self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name):
                                self.buttons[name].set_active(k)
                        except:
                            pass                            
                            
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.buttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding = 5)    
                           
                #strings without values
                if config.format[i].keys[j].kind == 'string' and len(config.format[i].keys[j].values) == 0 :
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.Entry()
                    self.buttons[name].connect("changed", self.saveCheck)
                    try:
                        self.buttons[name].set_text(self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name))
                    except:
                        pass
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox() 
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.buttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding=5)    
                    
                #ints & floats
                if config.format[i].keys[j].kind == 'int' or config.format[i].keys[j].kind == 'float':
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.Entry()
                    self.buttons[name].connect("changed", self.saveCheck)
                    try:
                        self.buttons[name].set_text(self.config.get('%s' %config.format[i].name, '%s' %config.format[i].keys[j].name))
                    except:
                        pass
                    nameLabel = gtk.Label(" %s:   " %config.format[i].keys[j].name)
                    hbox = gtk.HBox()
                    hbox.pack_start(infoButton, False, False)
                    hbox.pack_start(nameLabel, False, False)
                    Htable.attach(hbox,0,1,j,j+1, False | gtk.FILL, False)
                    Htable.attach(self.buttons[name],1,2,j,j+1, False | gtk.FILL, False, ypadding=5)   
                    
                #booleans
                if config.format[i].keys[j].kind == 'boolean':
                    name = config.format[i].keys[j].name
                    self.buttons[name] = gtk.RadioButton(label="True")
                    self.buttons[name].connect("clicked", self.saveCheck)
                    RB2 = gtk.RadioButton(label ="False", group=self.buttons[name])
                    try:
                        if self.config.get(' %s' %config.format[i].name, '%s' %config.format[i].keys[j].name) == 'True':
                            self.buttons[name].set_active(True)
                        else:
                            RB2.set_active(True)
                    except:
                        pass   
                    alignment = gtk.Alignment(yalign = .5)
                    nameLabel = gtk.Label("%s:   " %config.format[i].keys[j].name)
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
