##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##
##-----------------------------------------------------------------------------------
import ConfigParser
import logging
logger = logging.getLogger('config')
import StringIO
import numpy
import os.path
import sys
import string
import config

def init(config_file = ""):
    parser = ConfigParser.SafeConfigParser()

    parser.read(os.path.join(sys.path[0], 'default_config.ini'))

    gave_config = True
    if config_file != "":
        if os.path.isfile(config_file):
            parser.read(config_file)
        else:
            print >> sys.stderr, "specified configuration file %s does not exist" % config_file
            sys.exit(2)
    elif os.path.isfile('config.ini'):
        parser.read('config.ini')
        gave_config = False
    else:
        print >> sys.stderr, "You must provide a configuration file either by providing it as a command line argument or by placing a config.ini in the current directory."
        sys.exit(2)        

    #aKMC options
    config.akmc_temperature = parser.getfloat('aKMC', 'temperature')
    config.akmc_confidence  = parser.getfloat('aKMC', 'confidence')
    config.akmc_thermal_window = parser.getfloat('aKMC', 'thermally_accessible_window')
    config.akmc_max_thermal_window = parser.getfloat('aKMC', 'maximum_thermally_accessible_window') *akmc_thermal_window
    config.akmc_max_kmc_steps = parser.getint('aKMC', 'max_kmc_steps')

    #Debug Options
    config.debug_interactive_shell = parser.getboolean('Debug', 'interactive_shell')
    if debug_interactive_shell:
        import signal, code
        signal.signal(signal.SIGQUIT, lambda signum, frame: code.interact(local=locals()))
    config.debug_keep_bad_saddles  = parser.getboolean('Debug', 'keep_bad_saddles')
    config.debug_keep_all_results  = parser.getboolean('Debug', 'keep_all_result_files')
    try:
        config.debug_random_seed   = parser.getint('Debug', 'random_seed')
    except:
        config.debug_random_seed   = None
    if config.debug_random_seed:
        numpy.random.seed(config.debug_random_seed)
        logger.debug("Set random state from seed")
    config.debug_register_extra_results = parser.getboolean('Debug', 'register_extra_results')
    config.debug_list_search_results = parser.getboolean('Debug', 'list_search_results')
    config.debug_use_mean_time = parser.getboolean('Debug', 'use_mean_time')
    config.debug_target_trajectory = parser.get('Debug', 'target_trajectory')

    #path options
    config.path_root         = parser.get('Paths', 'main_directory')
    config.path_searches_out = parser.get('Paths', 'searches_out')
    config.path_searches_in  = parser.get('Paths', 'searches_in')
    config.path_states       = parser.get('Paths', 'states')
    config.path_results      = parser.get('Paths', 'results')
    config.path_pot          = parser.get('Paths', 'potential_files')

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
    config.comm_job_bundle_size = parser.getint('Communicator', 'job_bundle_size')
    config.comm_search_buffer_size = parser.getint('Communicator', 'search_buffer_size')
    config.path_scratch = parser.get('Paths', 'scratch')
    #print path_scratch
    if config.comm_type == 'local':
        config.comm_local_client = parser.get('Communicator', 'client_path')
        config.comm_local_ncpus = parser.getint('Communicator', 'number_of_CPUs')
    if config.comm_type == 'cluster':
        config.comm_script_path = parser.get('Communicator', 'script_path')
        config.comm_script_name_prefix = parser.get('Communicator', 'name_prefix')
        config.comm_script_queued_jobs_cmd = parser.get('Communicator', 'queued_jobs')
        config.comm_script_cancel_job_cmd = parser.get('Communicator', 'cancel_job')
        config.comm_script_submit_job_cmd = parser.get('Communicator', 'submit_job')
    if config.comm_type == 'mpi':
        config.comm_mpi_client = parser.get('Communicator', 'client_path')
        config.comm_mpi_mpicommand = parser.get('Communicator', 'mpi_command')
    if config.comm_type == 'boinc':
        config.comm_boinc_project_dir = parser.get('Communicator', 'boinc_project_dir')
        config.comm_boinc_wu_template_path = parser.get('Communicator', 'boinc_wu_template_path')
        config.comm_boinc_re_template_path = parser.get('Communicator', 'boinc_re_template_path')
        config.comm_boinc_appname = parser.get('Communicator', 'boinc_appname')
        config.comm_boinc_results_path = parser.get('Communicator', 'boinc_results_path')
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

    #
    #displacement options
    #

    #KDB
    config.kdb_on = parser.getboolean('KDB', 'use_kdb')
    if config.kdb_on:
        config.kdb_path = parser.get('Paths', 'kdb')
        config.kdb_addpath = parser.get('KDB', 'addpath')
        config.kdb_querypath = parser.get('KDB', 'querypath')
        config.kdb_wait = parser.get('KDB', 'wait')
        config.kdb_keep = parser.get('KDB', 'keep')
        config.kdb_rhsco = parser.getfloat('KDB', 'rhsco')


    #Recycling
    config.recycling_on = parser.getboolean('Recycling', 'use_recycling')
    config.recycling_save_sugg = parser.getboolean('Recycling', 'save_suggestions')
    if not config.recycling_on:
        config.disp_moved_only = False
    else:
        config.disp_moved_only = parser.getboolean('Recycling', 'displace_moved_only')
    config.recycling_move_distance = parser.getfloat('Recycling', 'move_distance')
    config.sb_recycling_on = parser.getboolean('Recycling','use_sb_recycling')
    if config.sb_recycling_on:
        config.sb_recycling_path = parser.get('Paths', 'superbasin_recycling')

    #Random Displacement
    config.disp_type = parser.get('Displacement', 'type')
    config.disp_brute_neighbors = parser.getboolean('Displacement', 'brute_neighbors')
    config.disp_cutoff = parser.getfloat('Displacement', 'cutoff')
    config.disp_use_covalent = parser.getboolean('Displacement', 'use_covalent')
    config.disp_covalent_scale = parser.getfloat('Displacement', 'covalent_scale')
    if config.disp_type == 'water':
        config.stdev_translation = parser.getfloat('Displacement', 'stdev_translation')
        config.stdev_rotation = parser.getfloat('Displacement', 'stdev_rotation')
    else:
        config.disp_size = parser.getfloat('Displacement', 'size')
        config.disp_radius = parser.getfloat('Displacement', 'radius')
    if config.disp_type == 'undercoordinated':
        config.disp_max_coord = parser.getint('Displacement', 'maximum_coordination')

    #Superbasins
    config.sb_on = parser.getboolean('Superbasins', 'use_superbasins')
    if config.sb_on:
        config.sb_path = parser.get('Paths', 'superbasins')
        config.sb_scheme = parser.get('Superbasins', 'scheme')
        if config.sb_scheme == 'transition_counting':
            config.sb_tc_ntrans = parser.getint('Superbasins', 'number_of_transitions')
        elif config.sb_scheme == 'energy_level':
            config.sb_el_energy_increment = parser.getfloat('Superbasins', 'energy_increment')

    config.askmc_on = parser.getboolean('Superbasins','use_askmc')
    if config.askmc_on:
        config.askmc_confidence = parser.getfloat('Superbasins','askmc_confidence')
        config.askmc_alpha = parser.getfloat('Superbasins','askmc_barrier_raise_param')
        config.askmc_gamma = parser.getfloat('Superbasins','askmc_high_barrier_def')
        config.askmc_barrier_test_on = parser.getboolean('Superbasins','askmc_barrier_test_on')
        config.askmc_connections_test_on = parser.getboolean('Superbasins','askmc_connections_test_on')

    #State comparison
    config.comp_eps_e = parser.getfloat('Structure Comparison', 'energy_difference')
    config.comp_eps_r = parser.getfloat('Structure Comparison', 'distance_difference')
    config.comp_use_identical = parser.getboolean('Structure Comparison', 'use_identical')

    del parser
