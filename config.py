##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

# Authors: Ian Johnson, Rye Terrell

import ConfigParser
import numpy
import os.path
import sys
import string
import config
import yaml

config.init_done = False

config.format = []

class ConfigSection:
    def __init__(self, name, description):
        self.name = name
        self.description = description
        self.keys = []
class ConfigKey:
    def __init__(self, name, kind, description, default):
        self.name = name
        self.kind = kind
        self.description = description
        self.default = default
        self.values = []
class ConfigValue:
    def __init__(self, name, description):
        self.name = name
        self.description = description

y = yaml.load(open(os.path.join(sys.path[0], "config.yaml"), 'r'))

for sectionName in y:
    section = ConfigSection(sectionName, y[sectionName]['description'])
    for key in y[sectionName]['options']:
        kattr = y[sectionName]['options'][key]
        ck = ConfigKey(key, kattr['kind'], kattr['description'], kattr['default'])
        section.keys.append(ck)
        if 'values' in kattr:
            for value in kattr['values']:
                ck.values.append(ConfigValue(value, kattr['values'][value]))
    config.format.append(section)


def init(config_file = ""):
    if config.init_done:
        return None
    config.init_done = True

    parser = ConfigParser.SafeConfigParser()
    
    for i in range(len(config.format)):
        parser.add_section(config.format[i].name)
        for j in config.format[i].keys:
            parser.set(config.format[i].name, j.name,str (j.default))

    gave_config = True
    if config_file != "":
        if os.path.isfile(config_file):
            parser.read(config_file)
            config.config_path = os.path.abspath(config_file)
        else:
            print >> sys.stderr, "specified configuration file %s does not exist" % config_file
            sys.exit(2)
    elif os.path.isfile('config.ini'):
        parser.read('config.ini')
        config.config_path = os.path.abspath('config.ini')
        gave_config = False
    else:
        print >> sys.stderr, "You must provide a configuration file either by providing it as a command line argument or by placing a config.ini in the current directory."
        sys.exit(2)        
    sections = False
    options = False

    #print parser.items("Main")

    psections = parser.sections()
    fsections = [j.name for j in config.format]
    for a in range(len(psections)):
        psections[a] = psections[a].lower()

    for a in range(len(fsections)):
        fsections[a] = fsections[a].lower()
    
    config_error = False
    for s in psections:
        if s not in fsections:
            config_error = True
            sys.stderr.write("unknown section %s\n")

    for i in parser.sections():
        b = parser.options(i)
        for a in range(len(b)):
            b[a] = b[a].lower()
        for k in config.format:
            if k.name == i:
                foptions = [j.name for j in k.keys]
                for a in range(len(foptions)):
                    foptions[a] = foptions[a].lower()
                for o in b:
                    if o not in foptions:
                        config_error = True
                        sys.stderr.write("unknown option %s in section %s\n" % (o, k.name))
      
    for psection in parser.sections():
        poptions = parser.options(psection)
        for fsection in config.format:
            if psection == fsection.name:
                for k in fsection.keys:
                    for o in poptions:
                        if o == k.name:
                            if k.kind == "int":
                                try:
                                    x = int(parser.get(psection,k.name))
                                except:
                                    config_error = True
                                    sys.stderr.write("option %s in section %s should be an integer\n" % (o,psection))
                            elif k.kind == "float":
                                try:
                                    x = float(parser.get(psection, k.name))
                                except:
                                    config_error = True
                                    sys.stderr.write("option %s of section %s should be a float\n" %(o,psection))
                            elif k.kind == "boolean":
                                booleans = ['True', 'true', 'T', 't', '0', 'False', 'false', 'F', 'f', '1']
                                if parser.get(psection,k.name) not in booleans:
                                    config_error = True
                                    sys.stderr.write("option %s of section %s should be boolean\n" %(o,psection))
                            elif k.kind == "string" and len(k.values) !=0:   
                                values = [m.name for m in k.values]
                                if parser.get(psection,k.name) not in values:
                                    Vnames = ", ".join([v.name for v in k.values])
                                    config_error = True
                                    sys.stderr.write("option %s should be one of: %s\n" %(parser.get(psection,k.name),Vnames))

    if config_error:
        sys.stderr.write("aborting: could not parse config.ini\n")
        sys.exit(1)

    #Main options
    config.main_job = parser.get('Main', 'job')
    config.main_temperature = parser.getfloat('Main', 'temperature')

    try:
        config.main_random_seed = parser.getint('Main', 'random_seed')
        numpy.random.seed(config.main_random_seed)
    except:
        config.main_random_seed = None

    #Structure Comparison options
    config.comp_eps_e = parser.getfloat('Structure Comparison', 'energy_difference')
    config.comp_eps_r = parser.getfloat('Structure Comparison', 'distance_difference')
    config.comp_use_identical = parser.getboolean('Structure Comparison', 'indistinguishable_atoms')
    config.comp_brute_neighbors = parser.getboolean('Structure Comparison', 'brute_neighbors')
    config.comp_neighbor_cutoff = parser.getfloat('Structure Comparison', 'neighbor_cutoff')
    config.comp_use_covalent = parser.getboolean('Structure Comparison', 'use_covalent')
    config.comp_covalent_scale = parser.getfloat('Structure Comparison', 'covalent_scale')

    #AKMC options
    config.akmc_confidence              = parser.getfloat('AKMC', 'confidence')
    config.akmc_thermal_window          = parser.getfloat('AKMC', 'thermally_accessible_window')
    config.akmc_max_thermal_window      = parser.getfloat('AKMC', 'thermally_accessible_buffer')
    config.akmc_max_kmc_steps           = parser.getint('AKMC', 'max_kmc_steps')
    config.akmc_confidence_scheme       = parser.get('AKMC', 'confidence_scheme')
    config.akmc_confidence_correction   = parser.getboolean('AKMC', "confidence_correction")

    #Basin Hopping options
    config.bh_md_probability = parser.getfloat('Basin Hopping', 'md_probability')
    
    #path options
    config.path_root         = parser.get('Paths', 'main_directory')
    config.path_jobs_out     = parser.get('Paths', 'jobs_out')
    config.path_jobs_in      = parser.get('Paths', 'jobs_in')
    config.path_incomplete   = parser.get('Paths', 'incomplete')
    config.path_states       = parser.get('Paths', 'states')
    config.path_results      = parser.get('Paths', 'results')
    config.path_pot          = parser.get('Paths', 'potential_files')
    config.path_bh_minima    = parser.get('Paths', 'bh_minima')

    #Rye-requested check
    if not gave_config and not os.path.samefile(config.path_root, os.getcwd()):
        res = raw_input("The config.ini file in the current directory does not point to the current directory. Are you sure you want to continue? (y/N) ").lower()
        if len(res)>0 and res[0] == 'y':
            pass
        else:
            sys.exit(3)
        
    #communicator options
    config.comm_type = parser.get('Communicator', 'type')
    #print comm_type
    config.comm_job_bundle_size = parser.getint('Communicator', 'jobs_per_bundle')
    config.comm_job_buffer_size = parser.getint('Communicator', 'num_jobs')
    config.path_scratch = parser.get('Paths', 'scratch')
    #print path_scratch
    if config.comm_type == 'local':
        config.comm_local_client = parser.get('Communicator', 'client_path')
        config.comm_local_client = parser.get('Communicator', 'client_path')
        config.comm_local_ncpus = parser.getint('Communicator', 'number_of_CPUs')
    if config.comm_type == 'cluster':
        config.comm_script_path = parser.get('Communicator', 'script_path')
        config.comm_script_name_prefix = parser.get('Communicator', 'name_prefix')
        config.comm_script_queued_jobs_cmd = parser.get('Communicator', 'queued_jobs')
        config.comm_script_cancel_job_cmd = parser.get('Communicator', 'cancel_job')
        config.comm_script_submit_job_cmd = parser.get('Communicator', 'submit_job')
    if config.comm_type == 'mpi':
        pass
    if config.comm_type == 'boinc':
        config.comm_boinc_project_dir = parser.get('Communicator', 'boinc_project_dir')
        config.comm_boinc_wu_template_path = parser.get('Communicator', 'boinc_wu_template_path')
        config.comm_boinc_re_template_path = parser.get('Communicator', 'boinc_re_template_path')
        config.comm_boinc_appname = parser.get('Communicator', 'boinc_appname')
        config.comm_boinc_results_path = parser.get('Communicator', 'boinc_results_path')
        config.comm_boinc_priority = parser.getint('Communicator', 'boinc_priority')
        #print comm_boinc_project_dir
        #print comm_boinc_wu_template_path
        #print comm_boinc_re_template_path
        #print comm_boinc_appname
        #print comm_boinc_results_path
    if config.comm_type == 'arc':
        if parser.has_option('Communicator', 'client_path'):
            config.comm_client_path = parser.get('Communicator', 'client_path')
        else:
            config.comm_client_path = ""

        if parser.has_option('Communicator', 'blacklist'):
            config.comm_blacklist = [ string.strip(c) for c in parser.get('Communicator', 'blacklist').split(',') ]
        else:
            config.comm_blacklist = []

    #Saddle Search options
    config.disp_type = parser.get('Saddle Search', 'displace_type')
    if config.disp_type == 'water':
        config.stdev_translation = parser.getfloat('Saddle Search', 'stdev_translation') # undocumented
        config.stdev_rotation = parser.getfloat('Saddle Search', 'stdev_rotation') # undocumented
        config.molecule_list = eval(parser.get('Saddle Search', 'molecule_list')) # undocumented
        config.disp_at_random = parser.getint('Saddle Search', 'disp_at_random') # undocumented
    else:
        config.disp_magnitude= parser.getfloat('Saddle Search', 'displace_magnitude')
        config.disp_radius = parser.getfloat('Saddle Search', 'displace_radius')
        config.disp_min_norm = parser.getfloat('Saddle Search', 'displace_min_norm')
    if config.disp_type == 'under_coordinated':
        config.disp_max_coord = parser.getint('Saddle Search', 'displace_max_coordination')
    if config.disp_type == 'listed_atoms':
        config.disp_listed_atoms = [ int(string.strip(c)) for c in parser.get('Saddle Search', 'displace_atomlist').split(',') ]

    #KDB
    config.kdb_on = parser.getboolean('KDB', 'use_kdb')
    config.kdb_only = parser.getboolean('KDB', 'kdb_only')
    config.kdb_scratch_path = parser.get('Paths', 'kdb_scratch')
    config.kdb_path = parser.get('Paths', 'kdb')
    config.kdb_wait = parser.get('KDB', 'wait')
    config.kdb_addpath = parser.get('KDB', 'addpath')
    if config.kdb_addpath == "False":
        config.kdb_addpath = os.path.join(os.path.dirname(__file__), "kdb", "kdbinsert.py")
    config.kdb_querypath = parser.get('KDB', 'querypath')
    if config.kdb_querypath == "False":
        config.kdb_querypath = os.path.join(os.path.dirname(__file__), "kdb", "kdbquery.py")

    #Recycling
    config.recycling_on = parser.getboolean('Recycling', 'use_recycling')
    config.recycling_save_sugg = parser.getboolean('Recycling', 'save_suggestions')
    if not config.recycling_on:
        config.disp_moved_only = False
    else:
        config.disp_moved_only = parser.getboolean('Recycling', 'displace_moved_only')
    config.recycling_move_distance = parser.getfloat('Recycling', 'move_distance')
    config.sb_recycling_on = parser.getboolean('Recycling','use_sb_recycling')
    config.sb_recycling_path = None
    if config.sb_recycling_on:
        config.sb_recycling_path = parser.get('Paths', 'superbasin_recycling')

    #Coarse Graining
    config.sb_on = parser.getboolean('Coarse Graining', 'use_projective_dynamics')
    config.sb_state_file = parser.get('Coarse Graining', 'state_file') 
    config.sb_path = None
    if config.sb_on:
        config.sb_path = parser.get('Paths', 'superbasins')
        config.sb_scheme = parser.get('Coarse Graining', 'superbasin_scheme')
        if config.sb_scheme == 'transition_counting':
            config.sb_tc_ntrans = parser.getint('Coarse Graining', 'number_of_transitions')
        elif config.sb_scheme == 'energy_level':
            config.sb_el_energy_increment = parser.getfloat('Coarse Graining', 'energy_increment')

    config.askmc_on = parser.getboolean('Coarse Graining','use_askmc')
    if config.askmc_on:
        config.askmc_confidence = parser.getfloat('Coarse Graining','askmc_confidence')
        config.askmc_alpha = parser.getfloat('Coarse Graining','askmc_barrier_raise_param')
        config.askmc_gamma = parser.getfloat('Coarse Graining','askmc_high_barrier_def')
        config.askmc_barrier_test_on = parser.getboolean('Coarse Graining','askmc_barrier_test_on')
        config.askmc_connections_test_on = parser.getboolean('Coarse Graining','askmc_connections_test_on')

    #Debug options
    config.debug_interactive_shell = parser.getboolean('Debug', 'interactive_shell')
    if config.debug_interactive_shell:
        import signal, code
        signal.signal(signal.SIGQUIT, lambda signum, frame: code.interact(local=locals()))
    config.debug_keep_bad_saddles  = parser.getboolean('Debug', 'keep_bad_saddles')
    config.debug_keep_all_results  = parser.getboolean('Debug', 'keep_all_result_files')
    config.debug_results_path = parser.get('Debug', 'result_files_path')
    config.debug_register_extra_results = parser.getboolean('Debug', 'register_extra_results')
    config.debug_use_mean_time = parser.getboolean('Debug', 'use_mean_time')
    config.debug_target_trajectory = parser.get('Debug', 'target_trajectory')

    del parser
