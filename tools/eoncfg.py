#!/usr/bin/env python

import os
import pygtk
import gtk
import gtk.glade
import ConfigParser
import pathfix



class eoncfg(object):
    def __init__(self):
        #creates a builder, adds glade file, creates an object for main window
        self.builder = gtk.Builder()
        self.builder.add_from_file(os.path.join(pathfix.path, "tools/eoncfg.glade"))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object("window1")
        #adds dialog windows
        self.jobDialog = self.builder.get_object("jobDialog")
        self.potentialDialog = self.builder.get_object("potentialDialog")
        self.tempDialog = self.builder.get_object("tempDialog")
        self.random_seedDialog = self.builder.get_object("random_seedDialog")
        self.energyDialog = self.builder.get_object("energyDialog")
        self.distanceDialog = self.builder.get_object("distanceDialog")
        self.ind_atomsDialog = self.builder.get_object("ind_atomsDialog")
        self.check_rotationDialog = self.builder.get_object("check_rotationDialog")
        self.neighborDialog = self.builder.get_object("neighborDialog")
        self.use_covalentDialog = self.builder.get_object("use_covalentDialog")
        self.covalent_scaleDialog = self.builder.get_object("covalent_scaleDialog")
        self.brute_neighborsDialog = self.builder.get_object("brute_neighborsDialog")
        self.typeDialog = self.builder.get_object("typeDialog")
        self.num_jobsDialog = self.builder.get_object("num_jobsDialog")
        self.jobs_per_bundleDialog = self.builder.get_object("jobs_per_bundleDialog")
        self.client_pathDialog = self.builder.get_object("client_pathDialog")
        self.number_of_cpusDialog = self.builder.get_object("number_of_cpusDialog")
        self.script_pathDialog = self.builder.get_object("script_pathDialog")
        self.name_prefixDialog = self.builder.get_object("name_prefixDialog")
        self.queued_jobsDialog = self.builder.get_object("queued_jobsDialog")
        self.submit_jobDialog = self.builder.get_object("submit_jobDialog")
        self.cancel_jobDialog = self.builder.get_object("cancel_jobDialog")
        self.project_dirDialog = self.builder.get_object("project_dirDialog")
        self.wu_template_pathDialog = self.builder.get_object("wu_template_pathDialog")
        self.re_template_pathDialog = self.builder.get_object("re_template_pathDialog")
        self.appnameDialog = self.builder.get_object("appnameDialog")
        self.results_pathDialog = self.builder.get_object("results_pathDialog")

        #adds info buttons
        self.jobInfoButton = self.builder.get_object("jobInfoButton")
        self.potentialInfoButton = self.builder.get_object("potentialInfoButton")
        self.tempInfoButton = self.builder.get_object("tempInfoButton")
        self.random_seedInfoButton = self.builder.get_object("random_seedInfoButton")
        self.energyInfoButton = self.builder.get_object("energyInfoButton")
        self.distanceInfoButton = self.builder.get_object("distanceInfoButton")
        self.ind_atomsInfoButton = self.builder.get_object("ind_atomsInfoButton")
        self.check_rotationInfoButton = self.builder.get_object("check_rotationInfoButton")
        self.neighborInfoButton = self.builder.get_object("neighborInfoButton")
        self.use_covalentInfoButton = self.builder.get_object("use_covalentInfoButton")
        self.covalent_scaleInfoButton = self.builder.get_object("covalent_scaleInfoButton")
        self.brute_neighborsInfoButton = self.builder.get_object("brute_neighborsInfoButton")
        self.typeInfoButton = self.builder.get_object("typeInfoButton")
        self.num_jobsInfoButton = self.builder.get_object("num_jobsInfoButton")
        self.jobs_per_bundleInfoButton = self.builder.get_object("jobs_per_bundleInfoButton")
        self.client_pathInfoButton = self.builder.get_object("client_pathInfoButton")
        self.number_of_cpusInfoButton = self.builder.get_object("number_of_cpusInfoButton")
        self.script_pathInfoButton = self.builder.get_object("script_pathInfoButton")
        self.name_prefixInfoButton = self.builder.get_object("name_prefixInfoButton")
        self.queued_jobsInfoButton = self.builder.get_object("queued_jobsInfoButton")
        self.submit_jobInfoButton = self.builder.get_object("submit_jobInfoButton")
        self.cancel_jobInfoButton = self.builder.get_object("cancel_jobInfoButton")
        self.project_dirInfoButton = self.builder.get_object("project_dirInfoButton")
        self.wu_template_pathInfoButton = self.builder.get_object("wu_template_pathInfoButton")
        self.re_template_pathInfoButton = self.builder.get_object("re_template_pathInfoButton")
        self.appnameInfoButton = self.builder.get_object("appnameInfoButton")
        self.results_pathInfoButton = self.builder.get_object("results_pathInfoButton")
        

        #function to run dialog window
        def rundialog(self, widget, dialog=None):
            dialog.show_all()
            result = dialog.run()
            dialog.hide()
        #attaches dialog buttons to rundialog function
        self.jobInfoButton.connect("button_press_event", rundialog, self.jobDialog)
        self.potentialInfoButton.connect("button_press_event", rundialog, self.potentialDialog)
        self.tempInfoButton.connect("button_press_event", rundialog, self.tempDialog)
        self.random_seedInfoButton.connect("button_press_event", rundialog, self.random_seedDialog)
        self.energyInfoButton.connect("button_press_event", rundialog, self.energyDialog)
        self.distanceInfoButton.connect("button_press_event", rundialog, self.distanceDialog)
        self.ind_atomsInfoButton.connect("button_press_event", rundialog, self.ind_atomsDialog)
        self.check_rotationInfoButton.connect("button_press_event", rundialog, self.check_rotationDialog)
        self.neighborInfoButton.connect("button_press_event", rundialog, self.neighborDialog)
        self.use_covalentInfoButton.connect("button_press_event", rundialog, self.use_covalentDialog)
        self.covalent_scaleInfoButton.connect("button_press_event", rundialog, self.covalent_scaleDialog)
        self.brute_neighborsInfoButton.connect("button_press_event", rundialog, self.brute_neighborsDialog)
        self.typeInfoButton.connect("button_press_event", rundialog, self.typeDialog)
        self.num_jobsInfoButton.connect("button_press_event", rundialog, self.num_jobsDialog)
        self.jobs_per_bundleInfoButton.connect("button_press_event", rundialog, self.jobs_per_bundleDialog)
        self.client_pathInfoButton.connect("button_press_event", rundialog, self.client_pathDialog)
        self.number_of_cpusInfoButton.connect("button_press_event", rundialog, self.number_of_cpusDialog)
        self.script_pathInfoButton.connect("button_press_event", rundialog, self.script_pathDialog)
        self.name_prefixInfoButton.connect("button_press_event", rundialog, self.name_prefixDialog)
        self.queued_jobsInfoButton.connect("button_press_event", rundialog, self.queued_jobsDialog)
        self.submit_jobInfoButton.connect("button_press_event", rundialog, self.submit_jobDialog)
        self.cancel_jobInfoButton.connect("button_press_event", rundialog, self.cancel_jobDialog)
        self.project_dirInfoButton.connect("button_press_event", rundialog, self.project_dirDialog)
        self.wu_template_pathInfoButton.connect("button_press_event", rundialog, self.wu_template_pathDialog)
        self.re_template_pathInfoButton.connect("button_press_event", rundialog, self.re_template_pathDialog)
        self.appnameInfoButton.connect("button_press_event", rundialog, self.appnameDialog)
        self.results_pathInfoButton.connect("button_press_event", rundialog, self.appnameDialog)
        
        #adds default config file if no config file exists
        self.config = ConfigParser.SafeConfigParser()
        self.config.read(os.path.join(pathfix.path, "default_config.ini"))
        try:
            self.config.read("./config.ini")
        except:
            print "No config.ini found in local directory, using default values."
        
        #stores default data from cfg file into objects
        temperature = self.config.get("Main", "temperature")
        random_seed = self.config.get("Main", "random_seed")
        try:
            job = self.config.get("Main", "job")
        except:
            job="akmc"
        try:
            potential = self.config.get("Main", "potential")
        except:
            potential = "eam_al"        

        commtype = self.config.get("Communicator", "type")
        num_jobs = self.config.get("Communicator", "num_jobs")
        jobs_per_bundle = self.config.get("Communicator", "jobs_per_bundle")
        client_path = self.config.get("Communicator", "client_path")
        num_cpus = self.config.get("Communicator", "number_of_CPUs")
        script_path = self.config.get("Communicator", "script_path")
        name_prefix = self.config.get("Communicator", "name_prefix")
        queued_jobs = self.config.get("Communicator", "queued_jobs")
        cancel_job = self.config.get("Communicator", "cancel_job")
        submit_job = self.config.get("Communicator", "submit_job")
        appname = self.config.get("Communicator", "boinc_appname")
        try:
            project_dir = self.config.get("Communicator", "boinc_project_dir")
        except:
            project_dir = ""
        try:
            wu_template_path = self.config.get("Communicator", "boinc_wu_template_path")
        except:
            wu_template_path = ""
        try:
            re_template_path = self.config.get("Communicator", "boinc_re_template_path")
        except:
            re_template_path = ""
            #results path gives error, default results path is currently set manuely
        try:
            results_path = self.config.get("Communicator", "boinc_results_path")
        except:
            results_path = "boinc_project_dir/sample_results"
        
        energy = self.config.get("Structure Comparison" , "energy_difference")
        distance = self.config.get("Structure Comparison", "distance_difference")
        neighbor = self.config.get("Structure Comparison", "neighbor_cutoff")
        covalent_scale = self.config.get("Structure Comparison", "covalent_scale")
        ind_atoms = self.config.get("Structure Comparison","indistinguishable_atoms")
        check_rotation = self.config.get("Structure Comparison", "check_rotation")
        use_covalent = self.config.get("Structure Comparison", "use_covalent")
        brute_neighbors = self.config.get("Structure Comparison", "brute_neighbors")
        
        #stores widgets from eoncfg.glade into objects
        self.tempEntry = self.builder.get_object("tempEntry")
        self.random_seedEntry = self.builder.get_object("random_seedEntry")
        self.energyEntry = self.builder.get_object("energyEntry")
        self.distanceEntry = self.builder.get_object("distanceEntry")
        self.neighborEntry = self.builder.get_object("neighborEntry")
        self.covalent_scaleEntry = self.builder.get_object("covalent_scaleEntry")
        self.num_jobsEntry = self.builder.get_object("num_jobsEntry")
        self.jobs_per_bundleEntry = self.builder.get_object("jobs_per_bundleEntry")
        self.client_pathEntry = self.builder.get_object("client_pathEntry")
        self.num_cpusEntry = self.builder.get_object("num_cpusEntry")
        self.script_pathEntry = self.builder.get_object("script_pathEntry")
        self.name_prefixEntry = self.builder.get_object("name_prefixEntry")
        self.queued_jobsEntry = self.builder.get_object("queued_jobsEntry")
        self.submit_jobEntry = self.builder.get_object("submit_jobEntry")
        self.cancel_jobEntry = self.builder.get_object("cancel_jobEntry")
        self.project_dirEntry = self.builder.get_object("project_dirEntry")
        self.wu_template_pathEntry = self.builder.get_object("wu_template_pathEntry")
        self.re_template_pathEntry = self.builder.get_object("re_template_pathEntry")
        self.appnameEntry = self.builder.get_object("appnameEntry")
        self.results_pathEntry = self.builder.get_object("results_pathEntry")
        self.jobButton = self.builder.get_object("jobButton")
        self.potentialButton = self.builder.get_object("potentialButton")        
        self.typeButton = self.builder.get_object("typeButton")
        self.ind_atomsTrue = self.builder.get_object("ind_atomsTrue")
        self.ind_atomsFalse = self.builder.get_object("ind_atomsFalse")
        self.check_rotationTrue = self.builder.get_object("check_rotationTrue")
        self.check_rotationFalse = self.builder.get_object("check_rotationFalse")
        self.use_covalentTrue = self.builder.get_object("use_covalentTrue")
        self.use_covalentFalse = self.builder.get_object("use_covalentFalse")
        self.brute_neighborsTrue = self.builder.get_object("brute_neighborsTrue")
        self.brute_neighborsFalse = self.builder.get_object("brute_neighborsFalse")
        

        # sets widgets to default to data from default cfg file
        self.tempEntry.set_text(temperature)
        self.random_seedEntry.set_text(random_seed)
        #.set_active can only take integers
        if commtype == "local":
            commtype = 0
        elif commtype == "cluster":
            commtype = 1
        elif commtype == "boinc":
            commtype = 2
        elif commtype == "arc":
            commtype = 3
        self.typeButton.set_active(commtype)
        if job == "akmc":
            job = 0
        elif job == "parallel_replica":
            job = 1
        elif job == "process_search":
            job = 2
        elif job == "saddle_search":
            job = 3
        elif job == "minimization":
            job = 4
        elif job == "hessian":
            job = 5
        elif job == "dimer_dr":
            job = 6
        elif job == "dimer_rotation":
            job = 7
        elif job == "displacement_sampling":
            job = 8
        else: 
            job = (-1)
        self.jobButton.set_active(job)
        if potential == "lj":
            potential = 0
        elif potential == "morse_pt":
            potential = 1
        elif potential == "emt":
            potential = 2
        elif potential == "epip":
            potential = 3
        elif potential == "vasp":
            potential = 4 
        elif potential == "tersoff_si":
            potential = 5 
        elif potential == "sw_si":
            potential = 6
        elif potential == "lenosky_si":
            potential = 7
        elif potential == "eam_al":
            potential = 8 
        elif potential == "qsc":
            potential = 9 
        elif potential == "zpice":
            potential = 10 
        elif potential == "tip4p":
            potential = 11 
        elif potential == "bopfox":
            potential = 12 
        else:
            potential = -1    
        self.potentialButton.set_active(potential)
        self.num_jobsEntry.set_text(num_jobs)
        self.jobs_per_bundleEntry.set_text(jobs_per_bundle)
        self.client_pathEntry.set_text(client_path)
        self.num_cpusEntry.set_text(num_cpus)
        self.script_pathEntry.set_text(script_path)
        self.name_prefixEntry.set_text(name_prefix)
        self.queued_jobsEntry.set_text(queued_jobs)
        self.cancel_jobEntry.set_text(cancel_job)
        self.submit_jobEntry.set_text(submit_job)
        self.appnameEntry.set_text(appname)
        self.results_pathEntry.set_text(results_path)
        self.project_dirEntry.set_text(project_dir)
        self.wu_template_pathEntry.set_text(wu_template_path)
        self.re_template_pathEntry.set_text(re_template_path)        
        self.energyEntry.set_text(energy)
        self.distanceEntry.set_text(distance)
        self.neighborEntry.set_text(neighbor)
        self.covalent_scaleEntry.set_text(covalent_scale)
        #if statements change radiobutton groups to display defaults
        if ind_atoms == "True":
            self.ind_atomsFalse.set_group(self.ind_atomsTrue)
        else:
            self.ind_atomsTrue.set_group(self.ind_atomsFalse)

        if check_rotation == "True":
            self.check_rotationFalse.set_group(self.check_rotationTrue)
        else:
            self.check_rotationTrue.set_group(self.check_rotationFalse)

        if use_covalent == "True":
            self.use_covalentFalse.set_group(self.use_covalentTrue)
        else:
            self.use_covalentTrue.set_group(self.use_covalentFalse)

        if brute_neighbors == "True":
            self.brute_neighborsFalse.set_group(self.brute_neighborsTrue)
        else:
            self.brute_neighborsTrue.set_group(self.brute_neighborsFalse)
        
        
        #adds save button from .glade and connects it to save function
        self.saveButton = self.builder.get_object("saveButton") 
        self.saveButton.connect("button_press_event", self.save)
        
        #shows window        
        self.window.show()

        

    #saves options to current directories config.ini file
    def save(self, widget, data=None):
        self.config.set('Main', 'temperature', self.tempEntry.get_text())
        self.config.set('Main', 'random_seed', self.random_seedEntry.get_text())
        self.config.set('Communicator', 'num_jobs', self.num_jobsEntry.get_text())
        self.config.set('Communicator', 'jobs_per_bundle', self.jobs_per_bundleEntry.get_text())
        self.config.set('Communicator', 'client_path', self.client_pathEntry.get_text())
        self.config.set('Communicator', 'number_of_CPUs', self.num_cpusEntry.get_text())
        self.config.set('Communicator', 'script_path', self.script_pathEntry.get_text())
        self.config.set('Communicator', 'name_prefix', self.name_prefixEntry.get_text())
        self.config.set('Communicator', 'queued_jobs', self.queued_jobsEntry.get_text())
        self.config.set('Communicator', 'cancel_job', self.cancel_jobEntry.get_text())
        self.config.set('Communicator', 'submit_job', self.submit_jobEntry.get_text())
        self.config.set('Communicator', 'boinc_appname', self.appnameEntry.get_text())
        self.config.set('Communicator', 'boinc_results_path', self.results_pathEntry.get_text())
        self.config.set('Communicator', 'boinc_project_dir', self.project_dirEntry.get_text())
        self.config.set('Communicator', 'boinc_wu_template_path', self.wu_template_pathEntry.get_text())
        self.config.set('Communicator', 'boinc_re_template_path', self.re_template_pathEntry.get_text())
        self.config.set('Structure Comparison', 'energy_difference', self.energyEntry.get_text())
        self.config.set('Structure Comparison', 'distance_difference', self.distanceEntry.get_text())
        self.config.set('Structure Comparison', 'neighbor_cutoff', self.neighborEntry.get_text())
        self.config.set('Structure Comparison', 'covalent_scale', self.covalent_scaleEntry.get_text())        
        #combo buttons are stored in ints and have to being changed back to strings
        newType = self.typeButton.get_active()
        if newType == 0:
            newType = "local"
        elif newType == 1:
            newType = "cluster"
        elif newType == 2:
            newType == "arc"
        self.config.set('Communicator', 'type', newType)
        newJob = self.jobButton.get_active()
        if newJob== 0:
            newJob ="akmc"
        elif newJob == 1:
            newJob ="parallel_replica"
        elif newJob == 2:
            newJob = "process_search"
        elif newJob == 3:
            newJob = "saddle_search"
        elif newJob == 4:
            newJob = "minimization"
        elif newJob == 5:
            newJob = "hessian"
        elif newJob == 6:
            newJob = "dimer_dr"
        elif newJob == 7:
            newJob = "dimer_rotation"
        elif newJob == 8:
            newJob = "displacement_sampling"
        else:
            newJob = ""
        self.config.set('Main', 'job', newJob)
        newPotential = self.potentialButton.get_active()
        if newPotential == 0:
            newPotential = "lj"
        elif newPotential == 1:
            newPotential = "morse_pt"
        elif newPotential == 2:
            newPotential = "emt"
        elif newPotential == 3:
            newPotential = "epip"
        elif newPotential == 4:
            newPotential = "vasp"
        elif newPotential == 5:
            newPotential = "tersoff_si" 
        elif newPotential == 6:
            newPotential = "sw_si"
        elif newPotential == 7:
            newPotential = "lenosky_si"
        elif newPotential == 8:
            newPotential = "eam_al"
        elif newPotential == 9:
            newPotential = "qsc" 
        elif newPotential == 10:
            newPotential = "zpice" 
        elif newPotential == 11:
            newPotential = "tip4p" 
        elif newPotential == 12:
            newPotential = "bopfox" 
        else:
            newPotential = ""
        self.config.set('Main', 'potential', newPotential)
        #if statements show which RB is active and returns an according string
        if self.ind_atomsTrue.get_active() == True:
            self.config.set('Structure Comparison', 'indistinguishable_atoms', 'True')
        else:
            self.config.set('Structure Comparison', 'indistinguishable_atoms', 'False')
        if self.check_rotationTrue.get_active() == True:
            self.config.set('Structure Comparison', 'check_rotation', 'True')
        else:
            self.config.set('Structure Comparison', 'check_rotation', 'False')
        if self.use_covalentTrue.get_active() == True:
            self.config.set('Structure Comparison', 'use_covalent', 'True')
        else:
            self.config.set('Structure Comparison', 'use_covalent', 'False')
        if self.brute_neighborsTrue.get_active() == True:
            self.config.set('Structure Comparison', 'brute_neighbors', 'True')           
        else:
            self.config.set('Structure Comparison', 'brute_neighbors', 'False')

        #changes label next to save button and allows config.set to work                   
        self.saveLabel = self.builder.get_object("saveLabel")
        self.saveLabel.set_text("Changes saved")
        f = open("config.ini", 'w')
        self.config.write(f)
        f.close()
        
            
        
      

        

            
    #allows program to exit when exit button is pressed
    def gtk_main_quit(self, widget):
        gtk.main_quit()
        
def main():
    GUI = eoncfg()
    gtk.main()
    
if __name__ == '__main__':
    pid = os.fork()
    if pid:
        os._exit(0)
    main()
