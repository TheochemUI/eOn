import ConfigParser
import StringIO
import os.path
import sys


parser = ConfigParser.SafeConfigParser()

parser.read(os.path.join(sys.path[0], 'default_config.ini'))

if len(sys.argv)>1:
    if os.path.isfile(sys.argv[-1]):
        parser.read(sys.argv[-1])
    else:
        print >> sys.stderr, "specified configuration file %s does not exist" % sys.argv[-1]
        sys.exit(2)
elif os.path.isfile('config.ini'):
    parser.read('config.ini')
else:
    print >> sys.stderr, "You must provide a configuration file either by providing it as a command line argument or by placing a config.ini in the current directory."
    sys.exit(2)        

#aKMC options
akmc_temperature = parser.getfloat('aKMC', 'temperature')
akmc_confidence  = parser.getfloat('aKMC', 'confidence')
akmc_thermal_window = parser.getfloat('aKMC', 'thermally_accessible_window')
akmc_max_thermal_window = parser.getfloat('aKMC', 'maximum_thermally_accessible_window') *akmc_thermal_window

#Debug Options
debug_keep_bad_saddles  = parser.getboolean('Debug', 'keep_bad_saddles')
debug_keep_all_results  = parser.getboolean('Debug', 'keep_all_result_files')
try:
    debug_random_seed   = parser.getint('Debug', 'random_seed')
except:
    debug_random_seed   = None

#path options
path_root         = parser.get('Paths', 'main_directory')
path_searches_out = parser.get('Paths', 'searches_out')
path_searches_in  = parser.get('Paths', 'searches_in')
path_states       = parser.get('Paths', 'states')
path_results      = parser.get('Paths', 'results')

#communicator options
comm_type = parser.get('Communicator', 'type')
#print comm_type
comm_job_buffer_size = parser.getint('Communicator', 'job_buffer_size')
path_scratch = parser.get('Paths', 'scratch')
#print path_scratch
if comm_type == 'local':
    comm_local_client = parser.get('Communicator', 'client_path')
    comm_local_ncpus = parser.getint('Communicator', 'number_of_CPUs')
if comm_type == 'boinc':
    comm_boinc_project_dir = parser.get('Communicator', 'boinc_project_dir')
    comm_boinc_wu_template_path = parser.get('Communicator', 'boinc_wu_template_path')
    comm_boinc_re_template_path = parser.get('Communicator', 'boinc_re_template_path')
    comm_boinc_appname = parser.get('Communicator', 'boinc_appname')
    comm_boinc_results_path = parser.get('Communicator', 'boinc_results_path')
    #print comm_boinc_project_dir
    #print comm_boinc_wu_template_path
    #print comm_boinc_re_template_path
    #print comm_boinc_appname
    #print comm_boinc_results_path


#
#displacement options
#

#KDB
kdb_on = parser.getboolean('KDB', 'use_kdb')
if kdb_on:
    kdb_path = parser.get('Paths', 'kdb')
    kdb_addpath = parser.get('KDB', 'addpath')
    kdb_querypath = parser.get('KDB', 'querypath')


#Recycling
recycling_on = parser.getboolean('Recycling', 'use_recycling')
recycling_move_distance = parser.getfloat('Recycling', 'move_distance')

#Random Displacement
disp_type = parser.get('Displacement', 'type')
disp_size = parser.getfloat('Displacement', 'size')
disp_radius = parser.getfloat('Displacement', 'radius')
if disp_type == 'undercoordinated':
    disp_max_coord = parser.getint('Displacement', 'maximum_coordination')

#Superbasins
sb_on = parser.getboolean('Superbasins', 'use_superbasins')
sb_scheme = parser.get('Superbasins', 'scheme')
if sb_scheme == 'transition_counting':
    sb_tc_ntrans = parser.getint('Superbasins', 'number_of_transitions')
    sb_path = parser.get('Paths', 'superbasins')

#State comparison
comp_eps_e = parser.getfloat('Structure Comparison', 'energy_difference')
comp_eps_r = parser.getfloat('Structure Comparison', 'distance_difference')
comp_use_identical = parser.getboolean('Structure Comparison', 'use_identical')

del parser

if __name__=='__main__':
    for i in dir():
        print i, eval(i)
