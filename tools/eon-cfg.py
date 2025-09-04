#!/usr/bin/env python

import os
import math
import time
import gtk.gdk as gdk
import gtk.glade as glade
import numpy as np
import pygtk
import gtk
import gtk.glade
import ConfigParser
import pathfix



class eoncfg(object):
    def __init__(self):
        #adds glade file, creates an object for main window
        gladetree = gtk.glade.XML(os.path.join(pathfix.path,"tools/eon-cfg.glade"))
        self.window = gladetree.get_widget("window1")
        #creates dialog windows
        self.jobDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, None)
        self.jobDialog.set_markup("""<b>Options:</b>
                                           <b>akmc:</b> Run an adaptive kinetic monte
                                                          carlo simulation.
                   <b>parallel_replica:</b> Calculate the rare-event dynamics
                                                          of the system by combining
                                                          transitions observed from multiple
                                                          trajectories run in parallel.
                    <b>process_search:</b> Combined saddle search,
                                                          minimizations, and prefactor
                                                          calculations. Used by the aKMC
                                                          method.
                      <b>saddle_search:</b> Do a saddle point search using a
                                                          minimum mode method.
                         <b>minimization:</b> Find the minimum from an initial
                                                          configuration.
                                      <b>hessian:</b> Calculate the Hessian matrix for
                                                          the specified configuration in a
                                                          process.
                                   <b>dimer_dr:</b> Rye is changing this.
                    <b>dimer_rotation:</b> Rye is changing this.
<b>displacement_sampling:</b> Job to sample different
                                                          displacement methods and
                                                          parameters to see
                                                          which are the most efficient.""")
        self.potentialDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, None)
        self.potentialDialog.set_markup("""<b>Options:</b>
                    <b>lj:</b> Lennard-Jones potential in reduced units ().
 <b>morse_pt:</b> Morse potential for platinum.
              <b>emt:</b> Effective medium theory, for metals.
              <b>edip:</b> Environment-Dependent Interatomic Potential,
                           for carbon.
             <b>vasp:</b> Vienna Ab-Initio Simulation Program (VASP)
                           interface.
 <b>tersoff_si:</b> Tersoff pair potential with angular terms, for
                           silicon.
            <b>sw_si:</b> Stillinger-Weber potential, for silicon.
<b>lenosky_si:</b> Lenosky potential, for silicon.
        <b>eam_al:</b> Embedded atom method parameterized for
                           aluminum.
                <b>qsc:</b> Quantum Sutton-Chen potential, for FCC metals.
            <b>zpice:</b> Water on platinum.
           <b>tip4p:</b> Point charge model for water.
       <b>bopfox:</b> Bond order potential, for metals.""")
        self.tempDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "The temperature of the simulation, in Kelvin.")
        self.random_seedDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Takes an integer number for the random seed. If this number is less than zero the current time is used as the random seed.")
        self.energyDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "How close in energy two configurations must be to be considered energetically equivalent.")
        self.distanceDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "The maximum distance two mapped atoms may be for two configurations to be considered equivalent.")
        self.ind_atomsDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Use an algorithm to compare structures that does not distinguish between atoms of the same element. That is to say the numbering of the atoms does not affect the structural comparison.")
        self.check_rotationDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Finds optimal overlap of structures via rotation before comparing them. Use this option in systems where structures can become rotated, such as nanoparticles.")
        self.neighborDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Atoms within this distance of each other are considered neighbors.")
        self.use_covalentDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Use the covalent radii of atoms to determine neighbors.")
        self.covalent_scaleDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT, gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Multiply covalent radii by this amount before determining neighbors.")
        self.brute_neighborsDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Determine neighbors by brute force (use this with nonorthogonal boxes).")
        self.typeDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, None)
        self.typeDialog.set_markup("""<b>Options:</b>
     <b>local:</b> The local communicator runs the calculations on the
                   same computer that the server is run on.
<b>cluster:</b> A job scheduler can be used to run jobs through
                   user supplied shell scripts. Examples are given for
                   SGE.
   <b>boinc:</b> Jobs can be submitted to a BOINC project.
        <b>arc:</b> Jobs can be submitted to the grid computing
                  software ARC.""")
        self.num_jobsDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, None)
        self.num_jobsDialog.set_markup("""     <b>local:</b> num_jobs' number of jobs will be run every time the
                   program is invoked.
<b>cluster:</b>  num_jobs is the desired sum of the queued and
                   running jobs. That is it should be set to the total
                   number of jobs that the user would like to run at once
                   and the script will submit jobs to the queue as
		   needed to try to achieve the target number.
   <b>boinc:</b>: one often wants to make enough work to keep all of
                   the clients busy. So instead of num_jobs being the
		   total number of jobs to run it sets the number of jobs
                   to keep in the queue. This way a buffer of
                   num_jobs/jobs_per_bundle
                   workunits are always kept in the BOINC queue.""")
        self.jobs_per_bundleDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "In eon a job is defined as task that the eon client executes, such as a process search or a parallel replica run. Sometimes it makes sense to run more than one of the same type of job at a time.")
        self.client_pathDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "Either the name or path to the eon client binary. If only a name and not a path is given then eon looks for the binary in same directory as config.ini failing to find it there it will search though the directories in the $PATH environment variable.")
        self.number_of_cpusDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "The number of jobs that will run simultaneously.")
        self.script_pathDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "The path to the user defined scripts for submitting jobs to the communicator.")
        self.name_prefixDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "When jobs are submitted to the scheduler they are given a unique internally used named. In order to make the jobs identifiable by the user the name_prefix can be set to a meaningful string that will always be prepended to the job names.")
        self.queued_jobsDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the name of the script that returns the job ids of all the running and queued jobs. It does not have to return the job ids of only eon related jobs.")
        self.submit_jobDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the name of the script that submits a single job to the queuing system. It takes two command line arguments. The first is the name of the job. This is not required for eon use, but is highly recommended so that users can identify which job is which. The second argument is the working directory. This is the path where the eon client should be executed. All of the needed client files will be placed in this directory. The script must return the job id of the submitted job. This is how eon internally keeps track of jobs.")
        self.cancel_jobDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the name of the script that cancels a job. It takes a single argument the job id.")
        self.project_dirDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the full path to the root of the BOINC project directory.")
        self.wu_template_pathDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the path, relative from the boinc_project_dir, to the boinc workunit template.")
        self.re_template_pathDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the path, relative from the boinc_project_dir, to the boinc result template.")
        self.appnameDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "boinc_appname Description")
        self.results_pathDialog = gtk.MessageDialog(self.window,gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,gtk.MESSAGE_INFO, gtk.BUTTONS_CLOSE, "This is the path where BOINC puts the final results. If you are using the sample_assimilator the results are stored in the project directory in a folder named sample_results.")


        #adds info buttons
        self.jobInfoButton = gladetree.get_widget("jobInfoButton")
        self.potentialInfoButton = gladetree.get_widget("potentialInfoButton")
        self.tempInfoButton = gladetree.get_widget("tempInfoButton")
        self.random_seedInfoButton = gladetree.get_widget("random_seedInfoButton")
        self.energyInfoButton = gladetree.get_widget("energyInfoButton")
        self.distanceInfoButton = gladetree.get_widget("distanceInfoButton")
        self.ind_atomsInfoButton = gladetree.get_widget("ind_atomsInfoButton")
        self.check_rotationInfoButton = gladetree.get_widget("check_rotationInfoButton")
        self.neighborInfoButton = gladetree.get_widget("neighborInfoButton")
        self.use_covalentInfoButton = gladetree.get_widget("use_covalentInfoButton")
        self.covalent_scaleInfoButton = gladetree.get_widget("covalent_scaleInfoButton")
        self.brute_neighborsInfoButton = gladetree.get_widget("brute_neighborsInfoButton")
        self.typeInfoButton = gladetree.get_widget("typeInfoButton")
        self.num_jobsInfoButton = gladetree.get_widget("num_jobsInfoButton")
        self.jobs_per_bundleInfoButton = gladetree.get_widget("jobs_per_bundleInfoButton")
        self.client_pathInfoButton = gladetree.get_widget("client_pathInfoButton")
        self.number_of_cpusInfoButton = gladetree.get_widget("number_of_cpusInfoButton")
        self.script_pathInfoButton = gladetree.get_widget("script_pathInfoButton")
        self.name_prefixInfoButton = gladetree.get_widget("name_prefixInfoButton")
        self.queued_jobsInfoButton = gladetree.get_widget("queued_jobsInfoButton")
        self.submit_jobInfoButton = gladetree.get_widget("submit_jobInfoButton")
        self.cancel_jobInfoButton = gladetree.get_widget("cancel_jobInfoButton")
        self.project_dirInfoButton = gladetree.get_widget("project_dirInfoButton")
        self.wu_template_pathInfoButton = gladetree.get_widget("wu_template_pathInfoButton")
        self.re_template_pathInfoButton = gladetree.get_widget("re_template_pathInfoButton")
        self.appnameInfoButton = gladetree.get_widget("appnameInfoButton")
        self.results_pathInfoButton = gladetree.get_widget("results_pathInfoButton")


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
        self.tempEntry = gladetree.get_widget("tempEntry")
        self.random_seedEntry = gladetree.get_widget("random_seedEntry")
        self.energyEntry = gladetree.get_widget("energyEntry")
        self.distanceEntry = gladetree.get_widget("distanceEntry")
        self.neighborEntry = gladetree.get_widget("neighborEntry")
        self.covalent_scaleEntry = gladetree.get_widget("covalent_scaleEntry")
        self.num_jobsEntry = gladetree.get_widget("num_jobsEntry")
        self.jobs_per_bundleEntry = gladetree.get_widget("jobs_per_bundleEntry")
        self.client_pathEntry = gladetree.get_widget("client_pathEntry")
        self.num_cpusEntry = gladetree.get_widget("num_cpusEntry")
        self.script_pathEntry = gladetree.get_widget("script_pathEntry")
        self.name_prefixEntry = gladetree.get_widget("name_prefixEntry")
        self.queued_jobsEntry = gladetree.get_widget("queued_jobsEntry")
        self.submit_jobEntry = gladetree.get_widget("submit_jobEntry")
        self.cancel_jobEntry = gladetree.get_widget("cancel_jobEntry")
        self.project_dirEntry = gladetree.get_widget("project_dirEntry")
        self.wu_template_pathEntry = gladetree.get_widget("wu_template_pathEntry")
        self.re_template_pathEntry = gladetree.get_widget("re_template_pathEntry")
        self.appnameEntry = gladetree.get_widget("appnameEntry")
        self.results_pathEntry = gladetree.get_widget("results_pathEntry")
        self.jobButton = gladetree.get_widget("jobButton")
        self.potentialButton = gladetree.get_widget("potentialButton")
        self.typeButton = gladetree.get_widget("typeButton")
        self.ind_atomsTrue = gladetree.get_widget("ind_atomsTrue")
        self.ind_atomsFalse = gladetree.get_widget("ind_atomsFalse")
        self.check_rotationTrue = gladetree.get_widget("check_rotationTrue")
        self.check_rotationFalse = gladetree.get_widget("check_rotationFalse")
        self.use_covalentTrue = gladetree.get_widget("use_covalentTrue")
        self.use_covalentFalse = gladetree.get_widget("use_covalentFalse")
        self.brute_neighborsTrue = gladetree.get_widget("brute_neighborsTrue")
        self.brute_neighborsFalse = gladetree.get_widget("brute_neighborsFalse")


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


        #adds save button & label from .glade and connects it to save function
        self.saveButton = gladetree.get_widget("saveButton")
        self.saveButton.connect("button_press_event", self.save)
        self.saveLabel = gladetree.get_widget("saveLabel")


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
        #changes save label
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
