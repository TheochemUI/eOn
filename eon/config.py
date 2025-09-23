import configparser
import numpy
import os.path
import sys
import string
import yaml


class ConfigSection:
    def __init__(self, name):
        self.name = name
        self.keys = []

class ConfigKey:
    def __init__(self, name, kind, default):
        self.name = name
        self.kind = kind
        self.default = default
        self.values = []

class ConfigClass:
    def __init__(self):
        self.init_done = False
        self.format = []

        yaml_file = open(os.path.join(os.path.dirname(__file__), 'config.yaml'))
        y = yaml.load(yaml_file, Loader=yaml.BaseLoader)
#        y = yaml.load(yaml_file)
        yaml_file.close()

        for sectionName in y:
            section = ConfigSection(sectionName)
            for key in y[sectionName]['options']:
                kattr = y[sectionName]['options'][key]
                ck = ConfigKey(key, kattr['kind'], kattr['default'])
                section.keys.append(ck)
                if 'values' in kattr:
                    for value in kattr['values']:
                        ck.values.append(value)
            self.format.append(section)

    def init(self, config_file = ""):
        if self.init_done:
            return None
        self.init_done = True

        parser = configparser.ConfigParser()

        for i in range(len(self.format)):
            parser.add_section(self.format[i].name)
            for j in self.format[i].keys:
                parser.set(self.format[i].name, j.name,str (j.default))

        gave_config = True
        if config_file != "":
            if os.path.isfile(config_file):
                parser.read(config_file)
                self.config_path = os.path.abspath(config_file)
            else:
                print("Specified configuration file %s does not exist" % config_file, sys.stderr)
                sys.exit(2)
        elif os.path.isfile('config.ini'):
            parser.read('config.ini')
            self.config_path = os.path.abspath('config.ini')
            gave_config = False
        else:
            print("You must provide a configuration file either by providing its name as a command line argument or by placing a config.ini in the current directory", sys.stderr)
            sys.exit(2)
        sections = False
        options = False

        # makes sections all lowercase
        psections = parser.sections()
        fsections = [j.name for j in self.format]
        for a in range(len(psections)):
            psections[a] = psections[a].lower()
        for a in range(len(fsections)):
            fsections[a] = fsections[a].lower()

        # checks all sections in config.ini are in configparser
        config_error = False
        for s in psections:
            if s not in fsections:
                config_error = True
                sys.stderr.write('unknown section "%s"\n' % s)

        # checks all options in config.ini are in configparser
        for i in parser.sections():
            b = parser.options(i)
            for a in range(len(b)):
                b[a] = b[a].lower()
            for k in self.format:
                if k.name == i:
                    foptions = [j.name for j in k.keys]
                    for a in range(len(foptions)):
                        foptions[a] = foptions[a].lower()
                    for o in b:
                        if o not in foptions:
                            config_error = True
                            sys.stderr.write('unknown option "%s" in section "%s"\n' % (o, k.name))

        # checks all options are the right type and if parser option has values config.ini is using one of those values
        for psection in parser.sections():
            poptions = parser.options(psection)
            for fsection in self.format:
                if psection == fsection.name:
                    for k in fsection.keys:
                        for o in poptions:
                            if o == k.name:
                                if k.kind == "int":
                                    try:
                                        x = int(parser.get(psection,k.name))
                                    except:
                                        config_error = True
                                        sys.stderr.write('option "%s" in section "%s" should be an integer\n' % (o,psection))
                                elif k.kind == "float":
                                    try:
                                        x = float(parser.get(psection, k.name))
                                    except:
                                        config_error = True
                                        sys.stderr.write('option "%s" of section "%s" should be a float\n' %(o,psection))
                                elif k.kind == "boolean":
                                    booleans = ['True', 'true', 'T', 't', '0', 'False', 'false', 'F', 'f', '1']
                                    if parser.get(psection,k.name) not in booleans:
                                        config_error = True
                                        sys.stderr.write('option "%s" of section "%s" should be boolean\n' %(o,psection))
                                elif k.kind == "string" and len(k.values) !=0:
                                    values = k.values
                                    if parser.get(psection,k.name) not in values:
                                        Vnames = ", ".join(k.values)
                                        config_error = True
                                        sys.stderr.write('option "%s" should be one of: %s\n' %(parser.get(psection,k.name),Vnames))

        if config_error:
            sys.stderr.write("aborting: could not parse config.ini\n")
            sys.exit(1)

        # Main options
        self.main_job = parser.get('Main', 'job')
        self.main_temperature = parser.getfloat('Main', 'temperature')
        self.main_checkpoint = parser.getboolean('Main', 'checkpoint')

        self.main_random_seed = parser.getint('Main', 'random_seed')

#       try:
#           self.main_random_seed = parser.getint('Main', 'random_seed')
#           numpy.random.seed(self.main_random_seed)
#       except:
#           self.main_random_seed = None

        # Structure Comparison options
        self.comp_eps_e = parser.getfloat('Structure Comparison', 'energy_difference')
        self.comp_eps_r = parser.getfloat('Structure Comparison', 'distance_difference')
        self.comp_use_identical = parser.getboolean('Structure Comparison', 'indistinguishable_atoms')
        self.comp_check_rotation = parser.getboolean('Structure Comparison', 'check_rotation')
        self.comp_brute_neighbors = parser.getboolean('Structure Comparison', 'brute_neighbors')
        self.comp_neighbor_cutoff = parser.getfloat('Structure Comparison', 'neighbor_cutoff')
        self.comp_use_covalent = parser.getboolean('Structure Comparison', 'use_covalent')
        self.comp_covalent_scale = parser.getfloat('Structure Comparison', 'covalent_scale')
        self.comp_remove_translation = parser.getboolean('Structure Comparison', 'remove_translation')

        # AKMC options
        self.akmc_confidence                 = parser.getfloat('AKMC', 'confidence')
        self.akmc_server_side_process_search = parser.getboolean('AKMC', 'server_side_process_search')
        if self.akmc_server_side_process_search:
            if parser.getfloat('Prefactor', 'default_value') == 0:
                print("Error: you must provide a default prefactor when using server-side process search mode.")
                sys.exit()
            if parser.getint('Communicator', 'jobs_per_bundle') != 1:
                print("Error: you cannot use a bundle size other than 1 when using server-side process search mode.")
                sys.exit()
        self.akmc_thermal_window = parser.getfloat('AKMC', 'thermally_accessible_window')
        self.akmc_max_thermal_window = parser.getfloat('AKMC', 'thermally_accessible_buffer')
        self.akmc_max_kmc_steps = parser.getint('AKMC', 'max_kmc_steps')
        self.akmc_confidence_scheme = parser.get('AKMC', 'confidence_scheme')
        self.akmc_confidence_correction = parser.getboolean('AKMC', "confidence_correction")
        self.akmc_max_rate = parser.getfloat('AKMC', "max_rate")
        self.akmc_eq_rate = parser.getfloat('AKMC', "eq_rate")

        # Basin Hopping options
        self.bh_initial_random_structure_probability = parser.getfloat('Basin Hopping', 'initial_random_structure_probability')
        self.bh_initial_state_pool_size = parser.getint('Basin Hopping', 'initial_state_pool_size')

        # Path options
        self.path_root         = parser.get('Paths', 'main_directory')
        self.path_jobs_out     = parser.get('Paths', 'jobs_out')
        self.path_jobs_in      = parser.get('Paths', 'jobs_in')
        self.path_incomplete   = parser.get('Paths', 'incomplete')
        self.path_states       = parser.get('Paths', 'states')
        self.path_results      = parser.get('Paths', 'results')
        self.path_pot          = parser.get('Paths', 'potential_files')
        self.path_bh_minima    = parser.get('Paths', 'bh_minima')

        # Rye-requested check
        # Should we have some kind of sanity-check module/function somewhere?
        if not gave_config and not os.path.samefile(self.path_root, os.getcwd()):
            res = input("The config.ini file in the current directory does not point to the current directory. Are you sure you want to continue? (y/N) ").lower()
            if len(res)>0 and res[0] == 'y':
                pass
            else:
                sys.exit(3)

        if int(self.main_random_seed) >= 0:
            if os.path.isfile(os.path.join(self.path_root, 'prng.pkl')):
                from eon import fileio as io
                io.get_prng_state()
            else:
                numpy.random.seed(self.main_random_seed)
        else:
            self.main_random_seed = None

        # Communicator options
        self.comm_type = parser.get('Communicator', 'type')
        self.comm_job_bundle_size = parser.getint('Communicator', 'jobs_per_bundle')
        self.comm_job_buffer_size = parser.getint('Communicator', 'num_jobs')
        self.comm_job_max_size = parser.getint('Communicator', 'max_jobs')
        self.path_scratch = parser.get('Paths', 'scratch')
        if self.comm_type == 'local':
            self.comm_local_client = parser.get('Communicator', 'client_path')
            self.comm_local_ncpus = parser.getint('Communicator', 'number_of_CPUs')
        if self.comm_type == 'cluster':
            self.comm_script_path = parser.get('Communicator', 'script_path')
            self.comm_script_name_prefix = parser.get('Communicator', 'name_prefix')
            self.comm_script_queued_jobs_cmd = parser.get('Communicator', 'queued_jobs')
            self.comm_script_cancel_job_cmd = parser.get('Communicator', 'cancel_job')
            self.comm_script_submit_job_cmd = parser.get('Communicator', 'submit_job')
        if self.comm_type == 'mpi':
            def mpiexcepthook(type, value, traceback):
                sys.__excepthook__(type, value, traceback)
                from mpi4py import MPI
                sys.stderr.write("exception occured on rank %i\n" % MPI.COMM_WORLD.rank);
                MPI.COMM_WORLD.Abort()
            sys.excepthook = mpiexcepthook

        # Process Search options
        self.process_search_minimization_offset = parser.getfloat('Process Search', 'minimization_offset')
        self.process_search_default_prefactor = parser.getfloat('Prefactor', 'default_value')
        self.process_search_minimize_first = parser.getboolean('Process Search', 'minimize_first')

        # Saddle Search options
        self.saddle_method = parser.get('Saddle Search', 'method')
        self.saddle_search_max_iterations = parser.getint('Saddle Search', 'max_iterations')
        self.saddle_dynamics_temperature = parser.getfloat('Saddle Search', 'dynamics_temperature')
        self.displace_random_weight = parser.getfloat('Saddle Search', 'displace_random_weight')
        self.displace_listed_atom_weight = parser.getfloat('Saddle Search', 'displace_listed_atom_weight')
        self.displace_listed_type_weight = parser.getfloat('Saddle Search', 'displace_listed_type_weight')
        self.displace_all_listed = parser.getboolean('Saddle Search', 'displace_all_listed')
        self.displace_under_coordinated_weight = parser.getfloat('Saddle Search', 'displace_under_coordinated_weight')
        self.displace_least_coordinated_weight = parser.getfloat('Saddle Search', 'displace_least_coordinated_weight')
        self.displace_not_FCC_HCP_weight = parser.getfloat('Saddle Search', 'displace_not_FCC_HCP_weight')
        self.displace_not_TCP_BCC_weight = parser.getfloat('Saddle Search', 'displace_not_TCP_BCC_weight')
        self.displace_not_TCP_weight = parser.getfloat('Saddle Search', 'displace_not_TCP_weight')
        self.displace_softest_mode_weight = parser.getfloat('Saddle Search', 'displace_softest_mode_weight')
        self.displace_water_weight = parser.getfloat('Saddle Search', 'displace_water_weight') # undocumented
        self.stdev_translation = parser.getfloat('Saddle Search', 'stdev_translation') # undocumented
        self.stdev_rotation = parser.getfloat('Saddle Search', 'stdev_rotation') # undocumented
        self.molecule_list = eval(parser.get('Saddle Search', 'molecule_list')) # undocumented
        self.disp_at_random = parser.getint('Saddle Search', 'disp_at_random') # undocumented
        self.disp_magnitude= parser.getfloat('Saddle Search', 'displace_magnitude')
        self.disp_radius = parser.getfloat('Saddle Search', 'displace_radius')
        self.disp_min_norm = parser.getfloat('Saddle Search', 'displace_min_norm')
        self.void_bias_fraction = parser.getfloat('Saddle Search', 'void_bias_fraction')
        self.disp_max_coord = parser.getint('Saddle Search', 'displace_max_coordination')
        self.random_mode = parser.getboolean('Saddle Search', 'random_mode')
        if self.displace_listed_atom_weight != 0.0:
            self.disp_listed_atoms = [ int(c.lstrip()) for c in parser.get('Saddle Search', 'displace_atom_list').split(',') ]
            if self.disp_listed_atoms == ['None']:
                self.disp_listed_atoms = []
        if self.displace_listed_type_weight != 0.0:
            self.disp_listed_types = [ (c.lstrip()) for c in parser.get('Saddle Search', 'displace_type_list').split(',') ]
            if self.disp_listed_types == ['None']:
                self.disp_listed_types = []
        self.displace_1d = parser.getboolean('Saddle Search', 'displace_1d')
        self.dynamics_max_init_curvature = parser.getfloat('Saddle Search', 'dynamics_max_init_curvature')
        self.zero_mode_abort_curvature = parser.getfloat('Saddle Search', 'zero_mode_abort_curvature')

        # KDB
        self.kdb_on = parser.getboolean('KDB', 'use_kdb')
        self.kdb_only = parser.getboolean('KDB', 'kdb_only')
        self.kdb_scratch_path = parser.get('Paths', 'kdb_scratch')
        self.kdb_path = parser.get('Paths', 'kdb')
        self.kdb_nodupes = parser.getboolean('KDB', 'remove_duplicates')
        self.kdb_name = parser.get('KDB', 'kdb_name')
        self.kdb_nf = parser.get('KDB', 'kdb_nf')
        self.kdb_dc = parser.get('KDB', 'kdb_dc')
        self.kdb_mac = parser.get('KDB', 'kdb_mac')

        # Recycling
        self.recycling_on = parser.getboolean('Recycling', 'use_recycling')
        self.recycling_save_sugg = parser.getboolean('Recycling', 'save_suggestions')
        if not self.recycling_on:
            self.disp_moved_only = False
        else:
            self.disp_moved_only = parser.getboolean('Recycling', 'displace_moved_only')
        self.recycling_move_distance = parser.getfloat('Recycling', 'move_distance')
        self.recycling_active_region = parser.getfloat('Recycling', 'active_region')
        self.recycling_mass_weight_factor = parser.getfloat('Recycling', 'mass_weight_factor')
        self.sb_recycling_on = parser.getboolean('Recycling','use_sb_recycling')
        self.sb_recycling_path = None
        if self.sb_recycling_on:
            self.sb_recycling_path = parser.get('Paths', 'superbasin_recycling')

        # Coarse Graining
        self.sb_on = parser.getboolean('Coarse Graining', 'use_mcamc')
        self.sb_state_file = parser.get('Coarse Graining', 'state_file') 
        self.sb_path = None
        self.sb_path = parser.get('Paths', 'superbasins')
        self.sb_scheme = parser.get('Coarse Graining', 'superbasin_scheme')
        self.sb_max_size = parser.getint('Coarse Graining', 'max_size')
        if self.sb_scheme == 'transition_counting':
            self.sb_tc_ntrans = parser.getint('Coarse Graining', 'number_of_transitions')
        elif self.sb_scheme == 'energy_level':
            self.sb_el_energy_increment = parser.getfloat('Coarse Graining', 'energy_increment')
        elif self.sb_scheme == 'rate':
            self.sb_rt_rate_threshold = parser.getfloat('Coarse Graining', 'rate_threshold')
        self.sb_superbasin_confidence = parser.getboolean('Coarse Graining', 'superbasin_confidence')

        self.askmc_on = parser.getboolean('Coarse Graining','use_askmc')
        if self.askmc_on:
            self.askmc_confidence = parser.getfloat('Coarse Graining','askmc_confidence')
            self.askmc_alpha = parser.getfloat('Coarse Graining','askmc_barrier_raise_param')
            self.askmc_gamma = parser.getfloat('Coarse Graining','askmc_high_barrier_def')
            self.askmc_barrier_test_on = parser.getboolean('Coarse Graining','askmc_barrier_test_on')
            self.askmc_connections_test_on = parser.getboolean('Coarse Graining','askmc_connections_test_on')

        # Optimizers
        self.optimizers_max_iterations = parser.getint('Optimizer', 'max_iterations')

        # Debug options
        self.debug_interactive_shell = parser.getboolean('Debug', 'interactive_shell')
        if self.debug_interactive_shell:
            import signal, code
            signal.signal(signal.SIGQUIT, lambda signum, frame: code.interact(local=locals()))
        self.debug_keep_bad_saddles = parser.getboolean('Debug', 'keep_bad_saddles')
        self.debug_keep_all_results = parser.getboolean('Debug', 'keep_all_result_files')
        self.debug_results_path = parser.get('Debug', 'result_files_path')
        self.debug_register_extra_results = parser.getboolean('Debug', 'register_extra_results')
        self.debug_use_mean_time = parser.getboolean('Debug', 'use_mean_time')
        self.debug_target_trajectory = parser.get('Debug', 'target_trajectory')
        self.debug_stop_criterion = parser.getfloat('Debug', 'stop_criterion')

        self.mpi_poll_period = parser.getfloat('Potential', 'mpi_poll_period')

        del parser

# XXX(rg): Very leaky global state here because.........?
config = ConfigClass()
